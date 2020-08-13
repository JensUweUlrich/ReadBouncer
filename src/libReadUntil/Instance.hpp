/*
 * Instance.hpp
 *
 * Header for Minknow API Instance
 *
 *  Created on: 28.10.2019
 *      Author: jens
 */
#include <string>
#include <grpcpp/grpcpp.h>
#include <minknow_api/instance.grpc.pb.h>

#include "MinKnowService.hpp"
#include "ReadUntilClientException.hpp"
#include "../debug_messages.hpp"

using namespace ::minknow_api::instance;

#ifndef LIBREADUNTIL_INSTANCE_HPP_
#define LIBREADUNTIL_INSTANCE_HPP_

namespace readuntil
{
	class Instance: public MinKnowService
	{
		private:
			std::unique_ptr<InstanceService::Stub> stub;
		public:
			Instance(std::shared_ptr<::grpc::Channel> channel);
			~Instance();

			std::string get_version_info();

	};
} // namespace readuntil

#endif // LIBREADUNTIL_INSTANCE_HPP_
