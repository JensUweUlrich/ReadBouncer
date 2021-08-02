/*
 * Instance.cpp
 *
 *  Created on: 08.11.2019
 *      Author: jens-uwe ulrich
 */

#include "Instance.hpp"

namespace readuntil
{
	Instance::Instance(std::shared_ptr<::grpc::Channel> channel)
	{
		stub = InstanceService::NewStub(channel);
	}

	std::string Instance::get_version_info()
	{
		GetVersionInfoRequest request;
		GetVersionInfoResponse response;
		::grpc::ClientContext context;
		::grpc::Status status = stub->get_version_info(&context, request, &response);
		if (status.ok())
		{
			return response.minknow().full();
		}
		else
		{
            throw ReadUntilClientException(status.error_message());
            //ReadUntilClientException(status.error_message());
		}
	}
}
