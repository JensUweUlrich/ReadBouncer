/*
 * Data.hpp
 *
 *  Created on: 12.11.2019
 *      Author: jens-uwe.ulrich
 */
#include <algorithm>
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
            string id{};
    };

    class Data: public MinKnowService
    {
        private:
            std::unique_ptr<DataService::Stub> stub;
            std::unique_ptr<grpc::ClientReaderWriter<GetLiveReadsRequest, GetLiveReadsResponse>> stream;
            std::queue<ReadCache> reads;
            std::vector<std::string> uniqueReadIds;
            std::mutex mutex;
            bool runs = false;
            uint8_t unblockChannels;
            uint8_t unblockReads;
            uint8_t actionBatchSize;
            void createSetupMessage();
            void getLiveSignals();
            void addActions();
            void addUnblockAction(GetLiveReadsRequest_Actions *actionList, ReadCache &read);
            void addStopReceivingDataAction(GetLiveReadsRequest_Actions *actionList, ReadCache &read);
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
            
            inline void setUnblockChannels(const uint8_t &unblock)
            {
                unblockChannels = unblock;
            }

            inline void setUnblockReads(const uint8_t &unblock)
            {
                unblockReads = unblock;
            }

            inline void setActionBatchSize(const uint8_t &size)
            {
                actionBatchSize = size;
            }
    };

} //namespace
#endif /* LIBREADUNTIL_DATA_HPP_ */
