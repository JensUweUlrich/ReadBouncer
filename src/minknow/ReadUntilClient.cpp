/*
 * ReadUntilClient.cpp
 *
 *  Created on: 28.10.2019
 *      Author: jens-uwe.ulrich
 */

#include "ReadUntilClient.hpp"


using json = nlohmann::json;



class Authentication : public grpc::MetadataCredentialsPlugin {
public:
    Authentication(const grpc::string& token) : m_token(token) {}

    grpc::Status GetMetadata(
        grpc::string_ref service_url, grpc::string_ref method_name,
        const grpc::AuthContext& channel_auth_context,
        std::multimap<grpc::string, grpc::string>* metadata
    ) override {
        metadata->insert(std::make_pair("local-auth", m_token));
        return grpc::Status::OK;
    }

private:
    grpc::string m_token;

	
};

namespace readuntil
{
	std::filesystem::path ReadUntilClientLog{};
	/**
	*	Establishes a connection to a given device/flowcell
	*	@device: Name of the flowcell to connect to
	*	@return: true if connection was sucessfully established, false otherwise
	*	@throws: DeviceServiceException, ReadUntilClientException
	*/
	bool ReadUntilClient::connect(std::string device, ConfigReader config)
	{
		ReadUntilClientLog /= "ReadUntilClientLog.txt";
		try
		{
			connection_logger = spdlog::rotating_logger_mt("RUClientLog", ReadUntilClientLog.string() , 1048576 * 5, 100);
		}
		catch(const spdlog::spdlog_ex& e)
		{
			std::cerr << "Log initialization failed: " << e.what() << std::endl;
		}
		
		connection_logger->flush_on(spdlog::level::info);
		connection_logger->set_level(spdlog::level::debug);

		channel_args = grpc::ChannelArguments();
		channel_args.SetMaxSendMessageSize(16 * 1024 * 1024);
		channel_args.SetMaxReceiveMessageSize(16 * 1024 * 1024);
		//channel_args.SetSslTargetNameOverride("localhost");
		uint32_t rpcPort{};

		std::stringstream info_str;
		bool secure_connect = false;
		
		if (mk_port == "9501")
		{
			info_str << "Connect to MinKNOW instance via unsecure connection to " << mk_host << " on port " << mk_port;
			connection_logger->info(info_str.str());
			connection_logger->flush();
			channel_creds = grpc::InsecureChannelCredentials();
			rpcPort = 8000;
		 }

		else // secure connection && mk_port = 9502;
		{
			info_str << "Connect to MinKNOW instance via secure SSL/TLS connection to " << mk_host << " on port " << mk_port;
			connection_logger->info(info_str.str());
			connection_logger->flush();
		
			//std::filesystem::path cert_file = "C:/ReadBouncer/src/rpc-certs/ca.crt";
			
			std::filesystem::path cert_file = ReadBouncerRoot;
			cert_file.append("rpc-certs");
			cert_file /= "ca.crt";
		
			if (!std::filesystem::exists(cert_file))
			{
				connection_logger->error("Could not find SSL/TLS certificate file : " + cert_file.string());
				connection_logger->flush();
				throw MissingCertificateException("Could not find SSL/TLS certificate file : " + cert_file.string());
			}

			std::ifstream ca(cert_file);
			std::string root_cert((std::istreambuf_iterator<char>(ca)),
			std::istreambuf_iterator<char>());
 
			grpc::SslCredentialsOptions opt;
			opt.pem_root_certs = root_cert;
			channel_creds = grpc::SslCredentials(opt); // for using with crt
			secure_connect = true;
			
			std::filesystem::path token_file;

			if (mk_host == "127.0.0.1"|| mk_host == "localhost"){

			std::stringstream s_1;
		 	s_1 << mk_host << ":" << mk_port;
			std::shared_ptr<::grpc::Channel> mgrCh_ = grpc::CreateCustomChannel(s_1.str(), channel_creds, channel_args);
			readuntil::Manager *mgr_ = new Manager(mgrCh_, secure_connect);
			token_file = mgr_->getTokenFilePath();
			}

			else{

				token_file = config.MinKNOW_Parsed.token_path;
				
			}

			if (!std::filesystem::exists(token_file))
			{
				connection_logger->error("Could not find token file : " + token_file.string());
				connection_logger->flush();
				throw MissingCertificateException("Could not find token file : " + token_file.string());
			}

			std::ifstream token_(token_file.string());
			std::string root_token((std::istreambuf_iterator<char>(token_)),
			std::istreambuf_iterator<char>());

			// read token from json file instead of substr()
			json options = json::parse(root_token);
			std::string token_auth = options["token"];


			std::cout << "Connecting to MinKNOW using the token: "<< token_auth << std::endl;

			Authentication *SSLConnection = new Authentication(token_auth);
			channel_creds = grpc::CompositeChannelCredentials(channel_creds,
		    	            grpc::MetadataCredentialsFromPlugin(std::unique_ptr<grpc::MetadataCredentialsPlugin>(SSLConnection)));

			std::stringstream s;
			s << mk_host << ":" << mk_port;

			std::shared_ptr<::grpc::Channel> mgrCh = grpc::CreateCustomChannel(s.str(), channel_creds, channel_args);
			readuntil::Manager* mgr = new Manager(mgrCh, secure_connect);

			// get RPC port for given device
			rpcPort = mgr->resolveRpcPort(device);

		}


		info_str.str("");
		info_str << "RPC port for " << device << " resolved";
		connection_logger->info(info_str.str());
		connection_logger->flush();
		std::stringstream connect_str;
		
		info_str.str("");
		connect_str << mk_host << ":" << rpcPort;
		info_str << "Trying to connect to Minknow on " << connect_str.str();
		connection_logger->info(info_str.str());
		connection_logger->flush();
		int retry_count = 5;

		for (int i = 1; i <= 5; ++i)
		{
			try
			{
				channel = grpc::CreateCustomChannel(connect_str.str(), channel_creds, channel_args);
				Instance* inst = (Instance*)getMinKnowService(MinKnowServiceType::INSTANCE);
				std::stringstream dm;
				 
				dm << "Sucessfully connected to minknow instance (version " << (*inst).get_version_info() << ")";
				connection_logger->info(dm.str());
				connected = true;
				break;
			}
			catch (ReadUntilClientException e)
			{
				std::stringstream em;
				em << "Failed to connect to minknow instance (retry " << i << "/" << retry_count << ") : " << e.what();
				connection_logger->error(em.str());
				connected = false;
				throw;
			}
			std::this_thread::sleep_for(std::chrono::seconds(1));
		}
		
		Device *dev = (readuntil::Device*) getMinKnowService(readuntil::MinKnowServiceType::DEVICE);
		try
		{
			std::stringstream devss;
			devss << "Detected " << (*dev).getDeviceType() << " Device with ID : " << (*dev).getDeviceId();
			connection_logger->info(devss.str());

		}
		catch (readuntil::DeviceServiceException e)
		{
			std::stringstream em;
			em << "Could not get device type/id : " << e.what();
			connection_logger->error(em.str());
			connected = false;
			throw;
		}

		return connected;
	}

	MinKnowService* ReadUntilClient::getMinKnowService(const MinKnowServiceType type)
	{
		
		MinKnowService *s;
		switch (type)
		{
			case ACQUISITION:
				return new Acquisition(channel);
			case ANALYSIS_CONFIGURATION:
				return new AnalysisConfiguration(channel);
			case INSTANCE:
				return new Instance(channel);
			case DATA:
				return new Data(channel);
			case DEVICE:
				return new Device(channel);
			default:
				return s;

		}
	}

}
