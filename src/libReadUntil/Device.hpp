/*
 * Device.hpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
 */

#include <string>
#include <grpcpp/grpcpp.h>
#include <minknow/rpc/device.grpc.pb.h>

#include "MinKnowService.hpp"
#include "ReadUntilClientException.hpp"

using namespace ::ont::rpc::device;

#ifndef LIBREADUNTIL_DEVICE_HPP_
#define LIBREADUNTIL_DEVICE_HPP_

namespace readuntil
{

	class Device: public MinKnowService
	{
		private:
			std::unique_ptr<DeviceService::Stub> stub;
		public:
			Device(std::shared_ptr<::grpc::Channel> channel);
			~Device();

	};

} //namespace

#endif /* LIBREADUNTIL_DEVICE_HPP_ */
