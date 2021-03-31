/*
 * ReadUntilClient.cpp
 *
 *  Created on: 28.10.2019
 *      Author: jens-uwe.ulrich
 */

#include "ReadUntilClient.hpp"

namespace readuntil
{
	/**
	*	Establishes a connection to a given device/flowcell
	*	@device: Name of the flowcell to connect to
	*	@return: true if connection was sucessfully established, false otherwise
	*	@throws: DeviceServiceException, ReadUntilClientException
	*/
	bool ReadUntilClient::connect(std::string device)
	{
		try
		{
			connection_logger = spdlog::rotating_logger_mt("RUClientLog", "logs/ReadUntilClientLog.txt", 1048576 * 5, 100);
		}
		catch(const spdlog::spdlog_ex& e)
		{
			std::cerr << "Log initialization failed: " << e.what() << std::endl;
		}
		


		// connect with MinKNOW Manager to get FlowCell Connection
		std::stringstream s;
		s << mk_host << ":" << mk_port;
		std::shared_ptr<::grpc::Channel> mgrCh = grpc::CreateChannel(s.str(), grpc::InsecureChannelCredentials());
		readuntil::Manager *mgr = new Manager(mgrCh);
		// get RPC port for given device
		uint32_t rpcPort = mgr->resolveRpcPort(device);


		std::stringstream connect_str;
		std::stringstream info_str;
		connect_str << mk_host << ":" << rpcPort;
		info_str << "Trying to connect to Minknow on " << connect_str.str();
		int retry_count = 5;
		connection_logger->flush_on(spdlog::level::info);
		connection_logger->set_level(spdlog::level::debug);
		connection_logger->info(info_str.str());

		for (int i = 1; i <= 5; ++i)
		{
			channel = grpc::CreateChannel(connect_str.str(), grpc::InsecureChannelCredentials());
			Instance *inst = (Instance*) getMinKnowService(MinKnowServiceType::INSTANCE);
			try
			{
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
