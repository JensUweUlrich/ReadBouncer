/*
 * Manager.cpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
 */

#include "Manager.hpp"

namespace readuntil
{
	Manager::Manager(std::shared_ptr<::grpc::Channel> channel)
	{
		stub = ManagerService::NewStub(channel);
	}

	std::vector<ListDevicesResponse::ActiveDevice> Manager::getActiveDevices()
	{
		ListDevicesRequest request;
		ListDevicesResponse response;
		::grpc::Status status = stub->list_devices(&context, request, &response);
		std::vector < ListDevicesResponse::ActiveDevice > actDev
		{ };
		if (status.ok())
		{
			actDev.assign(response.active().begin(), response.active().end());

		}
		else
		{
			throw ReadUntilClientException(status.error_message());
		}
		return actDev;
	}

	std::string Manager::getDeviceName(ListDevicesResponse::ActiveDevice &dev)
	{
		return dev.name();
	}

	uint32_t Manager::getRpcPort(ListDevicesResponse::ActiveDevice &dev)
	{
		return dev.ports().insecure_grpc();
	}

	std::vector<std::string> Manager::getPendingDevices()
	{
		ListDevicesRequest request;
		ListDevicesResponse response;
		::grpc::Status status = stub->list_devices(&context, request, &response);
		std::vector < ::std::string > pendDev
		{ };
		if (status.ok())
		{
			pendDev.assign(response.pending().begin(), response.pending().end());

		}
		else
		{
			throw ReadUntilClientException(status.error_message());
		}
		return pendDev;
	}

	std::vector<std::string> Manager::getInactiveDevices()
	{
		ListDevicesRequest request;
		ListDevicesResponse response;
		::grpc::Status status = stub->list_devices(&context, request, &response);
		std::vector < ::std::string > inactDev
		{ };
		if (status.ok())
		{
			inactDev.assign(response.inactive().begin(), response.inactive().end());

		}
		else
		{
			throw ReadUntilClientException(status.error_message());
		}
		return inactDev;
	}

}
