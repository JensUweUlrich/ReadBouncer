/*
 * Data.hpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
 */
#include <chrono>
#include <thread>
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
            bool runs = false;
            void createSetupMessage();
            void getLiveSignals();
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

            void getLiveReads();
            void addAction();
            bool isRunning();
    };

} //namespace
#endif /* LIBREADUNTIL_DATA_HPP_ */
