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

#include "DeepNano2.h"

#include "Acquisition.hpp"
#include "Analysis_Configuration.hpp"
#include "MinKnowService.hpp"
#include "DataServiceException.hpp"

using namespace ::minknow_api::data;
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
            std::string raw_data{};
    };

    struct ReadResponse
    {
        uint32 channelNr{};
        uint32 readNr{};
        string id{};
        uint8_t response{0};
        std::string sequence;
        uint64 samples_since_start{};
        double seconds_since_start{};
        uint64 start_sample{};
        uint64 chunk_start_sample{};
        uint64 chunk_length{};
    };


    class Data: public MinKnowService
    {
        private:
            std::unique_ptr<DataService::Stub> stub;
            std::unique_ptr<grpc::ClientReaderWriter<GetLiveReadsRequest, GetLiveReadsResponse>> stream;
            readuntil::Acquisition *acq;
            readuntil::AnalysisConfiguration *conf;
            std::queue<ReadCache> reads;
            std::queue<GetLiveReadsResponse_ActionResponse> responseQueue;
            std::map<string, ReadResponse> responseCache;
            std::vector<std::string> uniqueReadIds;
            std::set<int32> filterClasses;
            std::mutex readMutex;
            std::mutex responseMutex;
            std::mutex respQueueMutex;
            Caller *caller;
            std::shared_ptr<spdlog::logger> data_logger;
            bool runs = false;
            bool unblock_all = false;
            uint8_t unblockChannels;
            uint8_t unblockReads;
            uint8_t actionBatchSize = 50;
            void createSetupMessage();
            void getLiveSignals();
            void addActions();
            void addUnblockAction(GetLiveReadsRequest_Actions *actionList, ReadCache &read, const double unblock_duration);
            void addStopReceivingDataAction(GetLiveReadsRequest_Actions *actionList, ReadCache &read);
            void printResponseData();
            void adaptActionBatchSize();
            void resolveFilterClasses();
            std::vector<float> string_to_float(std::string const & s);
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

            void getLiveReads(std::string &weights);
            bool isRunning();
            void processSignals();
            
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

            inline void setUnblockAll(const bool &unblock)
            {
                unblock_all = unblock;
            }
    };

} //namespace
#endif /* LIBREADUNTIL_DATA_HPP_ */
