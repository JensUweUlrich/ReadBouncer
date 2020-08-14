/*
 * Analysis_Configuration.cpp
 *
 *  Created on: 07.04.2020
 *      Author: jens-uwe.ulrich
 */

#include "Analysis_Configuration.hpp"

namespace readuntil
{

    AnalysisConfiguration::AnalysisConfiguration(std::shared_ptr<::grpc::Channel> channel)
    {
        stub = AnalysisConfigurationService::NewStub(channel);
        Analysis_Configuration_logger = spdlog::get("RUClientLog");
    }

    Map<int32, string> AnalysisConfiguration::getReadClassifications()
    {
        GetReadClassificationsRequest request;
		GetReadClassificationsResponse response;
		::grpc::ClientContext context;
		::grpc::Status status = stub->get_read_classifications(&context, request, &response);
		if (status.ok())
		{
            return response.read_classifications();
		}
		else
		{
			throw ReadUntilClientException(status.error_message());
		}
    }

}

