/*
 * Data.cpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
 */

#include "Data.hpp"

namespace readuntil
{
	Data::Data(std::shared_ptr<::grpc::Channel> channel)
	{
		stub = DataService::NewStub(channel);
	}



}


