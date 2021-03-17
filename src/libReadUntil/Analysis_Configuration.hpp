/*
 * Analysis_Configuration.hpp
 *
 *  Created on: 14.08.2020
 *      Author: jens-uwe.ulrich
 */
#include <string>
#include <sstream>
#include <grpcpp/grpcpp.h>
#include <google/protobuf/map.h>
#include <minknow_api/analysis_configuration.grpc.pb.h>
#include <minknow_api/analysis_configuration.pb.h>

//#include "spdlog/spdlog.h"

#include "MinKnowService.hpp"
#include "ReadUntilClientException.hpp"

using namespace ::minknow_api::analysis_configuration;
using namespace ::google::protobuf;

#ifndef LIBREADUNTIL_ANALYSIS_CONFIGURATION_HPP_
#define LIBREADUNTIL_ANALYSIS_CONFIGURATION_HPP_

namespace readuntil
{

    class AnalysisConfiguration: public MinKnowService
    {
        private:
            std::unique_ptr<AnalysisConfigurationService::Stub> stub;
//            std::shared_ptr<spdlog::logger> Analysis_Configuration_logger;
            
        public:
            AnalysisConfiguration() = default;
            AnalysisConfiguration(std::shared_ptr<::grpc::Channel> channel);

            ~AnalysisConfiguration()
            {
                stub.release();
            }

            AnalysisConfiguration& operator=(const AnalysisConfiguration &other)
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
            // return map of read classifications mapped to id (used internally in analysis)
            Map<int32, string> getReadClassifications();
                
    };

} //namespace
#endif /* LIBREADUNTIL_ANALYSIS_CONFIGURATION_HPP_ */
