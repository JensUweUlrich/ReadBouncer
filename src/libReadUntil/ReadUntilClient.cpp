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
		DEBUGMESSAGE(std::cout, "Trying to connect to " << s.str());


		for (int i = 1; i <= 5 ; ++i)
		{
			channel = grpc::CreateChannel(s.str(), grpc::InsecureChannelCredentials());
			Instance *inst = (Instance*) getMinKnowService(MinKnowServiceType::INSTANCE);
			try
			{
				std::stringstream dm;
				dm << "Sucessfully connected to minknow instance (version " << (*inst).get_version_info() << ")";
				break;
			}
			catch (ReadUntilClientException e)
			{
				std::stringstream em;
				em << "Failed to connect to minknow instance (retry " << i << "/" << retry_count << ") : " << e.what();
				std::cerr <<  em.str() << std::endl;
			}
			std::this_thread::sleep_for(std::chrono::seconds(1));
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
				return new Manager(channel);
			default:
				return s;

		}
	}

}
