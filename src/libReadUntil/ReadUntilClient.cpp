/*
 * ReadUntilClient.cpp
 *
 *  Created on: 28.10.2019
 *      Author: jens
 */

#include "ReadUntilClient.hpp"

ReadUntilClient::ReadUntilClient()
{
	std::stringstream s;
	s << mk_host << ":" << mk_port;
	auto channel = grpc::CreateChannel(s.str(), grpc::InsecureChannelCredentials());

}

ReadUntilClient::~ReadUntilClient()
{
	// TODO Auto-generated destructor stub
}

