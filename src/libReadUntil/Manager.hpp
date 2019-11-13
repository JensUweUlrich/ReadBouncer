/*
 * Manager.hpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
 */

#include <string>
#include <vector>
#include <grpcpp/grpcpp.h>
#include <minknow/rpc/manager.grpc.pb.h>
#include <google/protobuf/repeated_field.h>

#include "MinKnowService.hpp"
#include "ReadUntilClientException.hpp"

using namespace ::ont::rpc::manager;

#ifndef LIBREADUNTIL_MANAGER_HPP_
#define LIBREADUNTIL_MANAGER_HPP_

namespace readuntil
{
	class Manager: public MinKnowService
	{
		private:
			std::unique_ptr<ManagerService::Stub> stub;
		public:
			Manager(std::shared_ptr<::grpc::Channel> channel);
			~Manager();

			std::vector<ListDevicesResponse::ActiveDevice> getActiveDevices();
			std::string getDeviceName(ListDevicesResponse::ActiveDevice &dev);
			uint32_t getRpcPort(ListDevicesResponse::ActiveDevice &dev);
			std::vector<std::string> getPendingDevices();
			std::vector<std::string> getInactiveDevices();

	};
} // namespace readuntil



#endif /* LIBREADUNTIL_MANAGER_HPP_ */
