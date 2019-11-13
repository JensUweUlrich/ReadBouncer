/*
 * Data.hpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
 */

#include <string>
#include <grpcpp/grpcpp.h>
#include <minknow/rpc/data.grpc.pb.h>

#include "MinKnowService.hpp"
#include "ReadUntilClientException.hpp"

using namespace ::ont::rpc::data;

#ifndef LIBREADUNTIL_DATA_HPP_
#define LIBREADUNTIL_DATA_HPP_

namespace readuntil
{

	class Data: public MinKnowService
	{
		private:
			std::unique_ptr<DataService::Stub> stub;
		public:
			Data(std::shared_ptr<::grpc::Channel> channel);
			~Data();
	};

} //namespace
#endif /* LIBREADUNTIL_DATA_HPP_ */
