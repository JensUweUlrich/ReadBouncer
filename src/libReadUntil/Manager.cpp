/*
 * Manager.cpp
 *
 * Created on: 12.11.2019
 * Author: Jens-Uwe Ulrich
 *
 * Fetching information about devices connected to MinKNOW via Remote Procedure Calls
 *
 * [update] 14.11.2019
 *		seems to use another port (9502?) and tls secure connection with rpc certificates from the MinKNOW conf directory
 *
 *
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
		std::vector<ListDevicesResponse::ActiveDevice> actDev
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
		::ont::rpc::manager::ListDevicesRequest request;
		::ont::rpc::manager::ListDevicesResponse response;
		::grpc::ClientContext context;
		::grpc::Status status = stub->list_devices(&context, request, &response);
		std::vector<::std::string> pendDev
		{ };
		if (status.ok())
		{
			if (response.pending_size() > 0)
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
		std::vector<::std::string> inactDev
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

	std::string Manager::getGuppyVersion()
	{
		::ont::rpc::manager::GetVersionInfoRequest request;
		::ont::rpc::manager::GetVersionInfoResponse response;
		::grpc::Status status = stub->get_version_info(&context, request, &response);
		if (status.ok())
		{
			return response.guppy_connected_version();
		}
		else
		{
			throw ReadUntilClientException(status.error_message());
		}
	}

}
