/*
 * Device.cpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
 */

#include "Device.hpp"

namespace readuntil
{
	Device::Device(std::shared_ptr<::grpc::Channel> channel)
	{
		stub = DeviceService::NewStub(channel);
	}



}


