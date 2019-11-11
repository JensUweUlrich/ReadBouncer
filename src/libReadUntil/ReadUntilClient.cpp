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
		uint8_t retry_count = 5;
		for (uint8_t i = 1; i <= 5 ; ++i)
		{
			channel = grpc::CreateChannel(s.str(), grpc::InsecureChannelCredentials());

			ReadUntilClient &client = ReadUntilClient::getClient();
			Instance *inst = (Instance*) client.getMinKnowService(MinKnowServiceType::INSTANCE);
			try
			{
				std::stringstream dm;
				dm << "Sucessfully connected to minknow instance (version " << (*inst).get_version_info() << ")" << std::endl;
				std::cout << dm.str() << std::endl;
				break;
			}
			catch (ReadUntilClientException e)
			{
				std::stringstream em;
				em << "Failed to connect to minknow instance (retry " << i << "/" << retry_count << " : " << e.what();
				std::cerr <<  em.str() << std::endl;
			}
			std::this_thread::sleep_for(std::chrono::seconds(1));
		}
	}

	MinKnowService* ReadUntilClient::getMinKnowService(const MinKnowServiceType type)
	{
		switch (type)
		{
			case INSTANCE:
				MinKnowService *s = new Instance(channel);
				return s;
		}
	}

}
