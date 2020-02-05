/*
 * Data.hpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
 */
#include <chrono>
#include <thread>
#include <string>
#include <sstream>
#include <queue>
#include <set>
#include <mutex>
#include <unistd.h>
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
    struct ReadCache
    {
            uint32 channelNr{};
            uint32 readNr{};
    };

    class Data: public MinKnowService
    {
        private:
            std::unique_ptr<DataService::Stub> stub;
            std::unique_ptr<grpc::ClientReaderWriter<GetLiveReadsRequest, GetLiveReadsResponse>> stream;
            std::queue<ReadCache> reads;
            std::mutex mutex;
            bool runs = false;
            void createSetupMessage();
            void getLiveSignals();
            void addActions();
        public:
            Data() = default;
            Data(std::shared_ptr<::grpc::Channel> channel);

            ~Data()
            {
                stub.release();
                //actionList.Clear();
            }

            Data& operator=(const Data &other)
            {
                if (this != &other)
                {
                    stub.reset(other.stub.get());
                    //actionList.CopyFrom(other.actionList);
                }
                return *this;
            }

	        grpc::ClientContext* getContext()
	        {
		        return &context;
	        }

            void getLiveReads();
            bool isRunning();
    };

} //namespace
#endif /* LIBREADUNTIL_DATA_HPP_ */
