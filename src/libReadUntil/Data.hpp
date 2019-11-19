/*
 * Data.hpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
 */

#include <string>
#include <grpcpp/grpcpp.h>
#include <google/protobuf/map.h>
#include <minknow/rpc/data.grpc.pb.h>
#include <minknow/rpc/data.pb.h>

#include "../debug_messages.hpp"

#include "MinKnowService.hpp"
#include "DataServiceException.hpp"

using namespace ::ont::rpc::data;
using namespace ::google::protobuf;

#ifndef LIBREADUNTIL_DATA_HPP_
#define LIBREADUNTIL_DATA_HPP_

namespace readuntil
{

	class Data: public MinKnowService
	{
		private:
			std::unique_ptr<DataService::Stub> stub;
			GetLiveReadsRequest_Actions actionList;
		public:
			Data(std::shared_ptr<::grpc::Channel> channel);
			~Data();

			void getLiveReads();
			void addAction();
	};

} //namespace
#endif /* LIBREADUNTIL_DATA_HPP_ */
