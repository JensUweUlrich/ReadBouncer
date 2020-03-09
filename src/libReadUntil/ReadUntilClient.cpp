/*
 * ReadUntilClient.cpp
 *
 *  Created on: 28.10.2019
 *      Author: jens-uwe.ulrich
 */

#include "ReadUntilClient.hpp"

namespace readuntil
{
	void ReadUntilClient::connect()
	{
		auto connection_logger = std::make_shared<spdlog::logger>("ClientConnection", daily_sink);
		spdlog::register_logger(connection_logger);
		std::stringstream s;
		s << "Trying to connect to " << mk_host << ":" << mk_port;
		int retry_count = 5;
		connection_logger->set_level(spdlog::level::debug);
		connection_logger->info(s.str());

		for (int i = 1; i <= 5; ++i)
		{
			channel = grpc::CreateChannel(s.str(), grpc::InsecureChannelCredentials());
			Instance *inst = (Instance*) getMinKnowService(MinKnowServiceType::INSTANCE);
			try
			{
				std::stringstream dm;
				dm << "Sucessfully connected to minknow instance (version " << (*inst).get_version_info() << ")";
				connection_logger->info(dm.str());
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
	}

	MinKnowService* ReadUntilClient::getMinKnowService(const MinKnowServiceType type)
	{
		MinKnowService *s;
		switch (type)
		{
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
