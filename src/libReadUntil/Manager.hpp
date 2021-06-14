/*
 * Manager.hpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
 */

#include <string>
#include <vector>
#include <grpcpp/grpcpp.h>
#include <minknow_api/manager.grpc.pb.h>
#include <google/protobuf/repeated_field.h>

#include "MinKnowService.hpp"
#include "ReadUntilClientException.hpp"

using namespace ::minknow_api::manager;

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

			std::vector<FlowCellPosition> getFlowCells();
			std::string getFlowCellName(FlowCellPosition &dev);
			uint32_t getRpcPort(FlowCellPosition &dev);
			std::string getGuppyVersion();
			uint32_t resolveRpcPort(std::string &deviceName);

	};
} // namespace readuntil



#endif /* LIBREADUNTIL_MANAGER_HPP_ */
