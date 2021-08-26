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

#if defined(_WIN32)
    #include <io.h>
#else
    #include <inttypes.h>
    #include <unistd.h>
#endif

#include <grpcpp/grpcpp.h>
#include <google/protobuf/map.h>
#include <minknow_api/data.grpc.pb.h>
#include <minknow_api/data.pb.h>

#include "spdlog/spdlog.h"

#include <uuid.h>

#include "SafeQueue.hpp"
#include "SafeMap.hpp"
#include "SafeSet.hpp"
#include "StopClock.hpp"
#include <CSVfile.hpp>

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
            TimeMeasures processingTimes{};
    };

    struct ActionResponse
    {
        uint32 channelNr{};
        uint32 readNr{};
        string id{};
        uint32_t readlen;
        TimeMeasures processingTimes{};
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
	        uuids::uuid_random_generator uuid_generator{};
            bool runs = false;
            bool unblock_all = false;
            uint8_t actionBatchSize = 50;
            
            void addUnblockAction(GetLiveReadsRequest_Actions* actionList, ActionResponse& response, const double unblock_duration);
            void addStopReceivingDataAction(GetLiveReadsRequest_Actions* actionList, ActionResponse& response);
            void adaptActionBatchSize(const int queue_size);
            void resolveFilterClasses();

            /**
            *   converts a given string of signals into a vector of float signals
            *   @signalString   : string of nanopore signals
            *   @return         : vector of float nanopore signals
            */
            inline std::vector<float> string_to_float(std::string const& signalString)
            {
                //assert(signalString.size() % sizeof(float) == 0);

                std::vector<float> result(signalString.size() / sizeof(float));

                if (!result.empty())
                {
                    std::copy(signalString.data(), signalString.data() + signalString.size(),
                        reinterpret_cast<char*>(&result.front()));
                }

                return result;
            }

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
            void sendActions(SafeQueue<readuntil::ActionResponse>& action_queue, SafeQueue<Durations>& duration_queue);
            void controlResponses(SafeQueue<readuntil::ActionResponse>& action_queue);
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
