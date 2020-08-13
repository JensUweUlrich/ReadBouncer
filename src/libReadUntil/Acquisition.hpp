/*
 * Acquisition.hpp
 *
 *  Created on: 07.04.2020
 *      Author: jens-uwe.ulrich
 */
#include <string>
#include <sstream>
#include <grpcpp/grpcpp.h>
#include <google/protobuf/map.h>
#include <minknow_api/acquisition.grpc.pb.h>
#include <minknow_api/acquisition.pb.h>

#include "spdlog/spdlog.h"

#include "MinKnowService.hpp"
#include "ReadUntilClientException.hpp"

using namespace ::minknow_api::acquisition;
using namespace ::google::protobuf;

#ifndef LIBREADUNTIL_ACQUISITION_HPP_
#define LIBREADUNTIL_ACQUISITION_HPP_

namespace readuntil
{

    class Acquisition: public MinKnowService
    {
        private:
            std::unique_ptr<AcquisitionService::Stub> stub;
            std::unique_ptr<grpc::ClientReaderWriter<WatchForStatusChangeRequest, WatchForStatusChangeResponse>> stream;
            std::shared_ptr<spdlog::logger> acquisition_logger;
            
        public:
            Acquisition() = default;
            Acquisition(std::shared_ptr<::grpc::Channel> channel);

            ~Acquisition()
            {
                stub.release();
            }

            Acquisition& operator=(const Acquisition &other)
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
            // return true immediately after MinKNOW status has changed to PROCESSING
            bool hasStarted();
            // return true immediately after MinKNOW status has changed to FINISHING
            bool isFinished();
            
    };

} //namespace
#endif /* LIBREADUNTIL_ACQUISITION_HPP_ */
