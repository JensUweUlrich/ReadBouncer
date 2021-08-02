/*
 * Analysis_Configuration.cpp
 *
 *  Created on	: 14.08.2020
 *	Last Change : 16.04.2021 
 *      Author	: jens-uwe.ulrich
 */

#include "Analysis_Configuration.hpp"

namespace readuntil
{
	/**
	*   Constructor of Data class
	*   @channel:   shared pointer to already opened GRPC channel for communication with MinKNOW
	*/
    AnalysisConfiguration::AnalysisConfiguration(std::shared_ptr<::grpc::Channel> channel)
    {
        stub = AnalysisConfigurationService::NewStub(channel);
        Analysis_Configuration_logger = spdlog::get("RUClientLog");
    }

	/**
	*	set the seconds after which MinKNOW breaks the read signals and streams them to the ReadUntil client
	*	one chunk corresponds to those seconds
	*	@seconds	: seconds after which to break the reads
	*/
	void AnalysisConfiguration::set_break_reads_after_seconds(double seconds)
	{
		GetAnalysisConfigurationRequest request;
		SetAnalysisConfigurationResponse response;
		minknow_api::analysis_configuration::AnalysisConfiguration conf;
		
		::grpc::Status status = stub->get_analysis_configuration(&context, request, &conf);

		if (!status.ok())
            throw ReadUntilClientException(status.error_message());
        //ReadUntilClientException(status.error_message());
		
		ReadDetectionParams current_rd = conf.read_detection();
		std::stringstream sstr;

		ReadDetectionParams* new_rd = conf.mutable_read_detection();
		DoubleValue* s = new_rd->mutable_break_reads_after_seconds();
		sstr.str("");
		sstr << "Applications default value for break_reads_after_seconds : " << s->value();
		Analysis_Configuration_logger->info(sstr.str());

		s->set_value(seconds);
		::grpc::ClientContext context2;
		status = stub->set_analysis_configuration(&context2, conf, &response);
		if (!status.ok())
            throw ReadUntilClientException(status.error_message());
        //ReadUntilClientException(status.error_message());

		::grpc::ClientContext context3;
		status = stub->get_analysis_configuration(&context3, request, &conf);
		if (!status.ok())
            //ReadUntilClientException(status.error_message());
            throw ReadUntilClientException(status.error_message());

		sstr.str("");
		sstr << "Set new value for break_reads_after_seconds to : " << conf.read_detection().break_reads_after_seconds().value();
		Analysis_Configuration_logger->info(sstr.str());
		Analysis_Configuration_logger->flush();
	}

	/**
	*	get all possible classifications for reads
	*/
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
            //ReadUntilClientException(status.error_message());
            throw ReadUntilClientException(status.error_message());
		}
    }

}

