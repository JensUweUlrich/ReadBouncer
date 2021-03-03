/*
 * Data.hpp
 *
 *  Created on: 12.11.2019
 *      Author: jens-uwe.ulrich
 */

#pragma once

#include <algorithm>
#include <chrono>
#include <thread>
#include <string>
#include <sstream>
#include <fstream>
#include <queue>
#include <set>
#include <mutex>
#include <io.h>
#include <grpcpp/grpcpp.h>
#include <google/protobuf/map.h>
#include <minknow_api/data.grpc.pb.h>
#include <minknow_api/data.pb.h>

#include "spdlog/spdlog.h"

#include "SafeQueue.hpp"
#include "StopClock.hpp"

#include "Acquisition.hpp"
#include "Analysis_Configuration.hpp"
#include "MinKnowService.hpp"
#include "DataServiceException.hpp"

using namespace ::minknow_api::data;
using namespace ::google::protobuf;


namespace readuntil
{
    struct SignalRead
    {
            uint32 channelNr{};
            uint32 readNr{};
            string id{};
            std::vector<float> raw_signals{};

    };

    struct ActionResponse
    {
        uint32 channelNr{};
        uint32 readNr{};
        string id{};
        bool response;
    };


    class Data: public MinKnowService
    {
        private:
            std::unique_ptr<DataService::Stub> stub;
            std::unique_ptr<grpc::ClientReaderWriter<GetLiveReadsRequest, GetLiveReadsResponse>> stream;
            readuntil::Acquisition* acq;
            readuntil::AnalysisConfiguration* conf;
            std::set<int32> filterClasses;
            std::shared_ptr<spdlog::logger> data_logger;
            bool runs = false;
            bool unblock_all = false;
            uint8_t actionBatchSize = 50;
            
            void addUnblockAction(GetLiveReadsRequest_Actions* actionList, ActionResponse& response, const double unblock_duration);
            void addStopReceivingDataAction(GetLiveReadsRequest_Actions* actionList, ActionResponse& response);
            void adaptActionBatchSize(const int queue_size);
            void resolveFilterClasses();
            std::vector<float> string_to_float(std::string const & s);
        public:
            Data() = default;
            Data(std::shared_ptr<::grpc::Channel> channel);

            ~Data()
            {
                stub.release();
                delete acq;
                delete conf;
                //actionList.Clear();
            }

            Data& operator=(const Data &other)
            {
                if (this != &other)
                {
                    stub.reset(other.stub.get());
                }
                return *this;
            }

	        grpc::ClientContext* getContext()
	        {
		        return &context;
	        }

            void getLiveSignals(SafeQueue<SignalRead>& basecall_queue);
            void sendActions(SafeQueue<readuntil::ActionResponse>& action_queue);
            void startLiveStream();
            void stopLiveStream();

            bool isRunning();
            
            inline void setActionBatchSize(const uint8_t &size)
            {
                actionBatchSize = size;
            }

            inline void setUnblockAll(const bool &unblock)
            {
                unblock_all = unblock;
            }
    };

} //namespace
