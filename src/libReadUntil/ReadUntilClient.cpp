/*
 * ReadUntilClient.cpp
 *
 *  Created on: 28.10.2019
 *      Author: jens
 */

#include "ReadUntilClient.hpp"

namespace readuntil
{
	void ReadUntilClient::connect()
	{
		std::stringstream s;
		s << mk_host << ":" << mk_port;
		int retry_count = 5;
		std::cout << "Trying to connect to " << s.str() << std::endl;

		for (int i = 1; i <= 5; ++i)
		{
			channel = grpc::CreateChannel(s.str(), grpc::InsecureChannelCredentials());
			Instance *inst = (Instance*) getMinKnowService(MinKnowServiceType::INSTANCE);
			try
			{
				std::stringstream dm;
				dm << "Sucessfully connected to minknow instance (version " << (*inst).get_version_info() << ")";
				std::cout << dm.str() << std::endl;
				break;
			}
			catch (ReadUntilClientException e)
			{
				std::stringstream em;
				em << "Failed to connect to minknow instance (retry " << i << "/" << retry_count << ") : " << e.what();
				std::cerr << em.str() << std::endl;
			}
			std::this_thread::sleep_for(std::chrono::seconds(1));
		}

		Device *dev = (readuntil::Device*) getMinKnowService(readuntil::MinKnowServiceType::DEVICE);
		try
		{

			std::cout << "Detected " << (*dev).getDeviceType() << " Device with ID : " << (*dev).getDeviceId() << std::endl;

		}
		catch (readuntil::DeviceServiceException e)
		{
			std::cerr << "Could not get device type/id : " << e.what() << std::endl;
		}

		readuntil::Manager *mgr = (readuntil::Manager*) getMinKnowService(readuntil::MinKnowServiceType::MANAGER);

		try
		{

			std::cout << "guppy version : " << (*mgr).getGuppyVersion() << std::endl;

		}
		catch (readuntil::ReadUntilClientException e)
		{
			std::cerr << "Could not get guppy version : " << e.what() << std::endl;
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
