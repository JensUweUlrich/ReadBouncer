/*
 * ReadUntilClient.cpp
 *
 *  Created on: 28.10.2019
 *      Author: jens-uwe.ulrich
 */

#include "ReadUntilClient.hpp"

namespace readuntil
{
	bool ReadUntilClient::connect()
	{
		try
		{
			connection_logger = spdlog::rotating_logger_mt("RUClientLog", "logs/ReadUntilClientLog.txt", 1048576 * 5, 100);
		}
		catch(const spdlog::spdlog_ex& e)
		{
			std::cerr << "Log initialization failed: " << e.what() << std::endl;
		}
		
		bool connected = false;
		std::stringstream connect_str;
		std::stringstream info_str;
		connect_str << mk_host << ":" << mk_port;
		info_str << "Trying to connect to Minknow on " << connect_str.str();
		int retry_count = 5;
		connection_logger->flush_on(spdlog::level::err);
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
		}

		readuntil::Manager *mgr = (readuntil::Manager*) getMinKnowService(readuntil::MinKnowServiceType::MANAGER);

		try
		{
			std::stringstream gupss;
			gupss << "guppy version : " << (*mgr).getGuppyVersion();
			connection_logger->info(gupss.str());
		}
		catch (readuntil::ReadUntilClientException e)
		{
			std::stringstream em;
			em << "Could not get guppy version : " << e.what();
			connection_logger->error(em.str());
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
			case INSTANCE:
				return new Instance(channel);
			case DATA:
				return new Data(channel);
			case DEVICE:
				return new Device(channel);
			case MANAGER:
			{
				std::stringstream s;
				s << mk_host << ":" << "9501";
				std::shared_ptr<::grpc::Channel> mgrCh = grpc::CreateChannel("localhost:9501", grpc::InsecureChannelCredentials());
				return new Manager(mgrCh);
			}
			default:
				return s;

		}
	}

}
