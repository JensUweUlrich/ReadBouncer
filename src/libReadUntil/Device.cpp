/*
 * Device.cpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
 */

#include "Device.hpp"

namespace readuntil
{

	Device::Device(std::shared_ptr<::grpc::Channel> channel)
	{
		stub = DeviceService::NewStub(channel);
		deviceInfo = new GetDeviceInfoResponse();
		deviceState = new GetDeviceStateResponse();
		flowCellInfo = new GetFlowCellInfoResponse();
	}

	std::string Device::getDeviceId()
	{
		if (!deviceInfoLoaded)
		{
			getDeviceInfo();
		}
		return deviceInfo->device_id();

	}

	std::string Device::getDeviceType()
	{
		if (!deviceInfoLoaded)
		{
			getDeviceInfo();
		}

		switch (deviceInfo->device_type())
		{
			case GetDeviceInfoResponse::DeviceType::GetDeviceInfoResponse_DeviceType_MINION:
				return "MinION";
			case GetDeviceInfoResponse::DeviceType::GetDeviceInfoResponse_DeviceType_MINION_MK1C:
				return "MinION Mk1C";
			case GetDeviceInfoResponse::DeviceType::GetDeviceInfoResponse_DeviceType_GRIDION:
				return "GridION";
			case GetDeviceInfoResponse::DeviceType::GetDeviceInfoResponse_DeviceType_PROMETHION:
				return "PromethION";
			default:
				return "Unknown";
		}

	}

	bool Device::isReady()
	{
		if (!deviceStateLoaded)
		{
			getDeviceState();
		}

		if (deviceState->device_state() == GetDeviceStateResponse::DeviceState::GetDeviceStateResponse_DeviceState_DEVICE_READY)
			return true;
		else
			return false;

	}

	bool Device::isDisconnected()
	{
		if (!deviceStateLoaded)
		{
			getDeviceState();
		}
		if (deviceState->device_state() == GetDeviceStateResponse::DeviceState::GetDeviceStateResponse_DeviceState_DEVICE_DISCONNECTED)
			return true;
		else
			return false;

	}

	bool Device::hasFlongleAdapter()
	{
		if (!flowCellInfoLoaded)
		{
			getFlowcellInfo();
		}
		return flowCellInfo->has_adapter();

	}

	bool Device::hasFlowCell()
	{
		if (!flowCellInfoLoaded)
		{
			getFlowcellInfo();
		}
		return flowCellInfo->has_flow_cell();
	}

	std::unique_ptr<grpc::ClientReader<GetDeviceStateResponse>> Device::streamDeviceState()
	{
		StreamDeviceStateRequest request;
		std::unique_ptr<grpc::ClientReader<GetDeviceStateResponse>> reader(stub->stream_device_state(&context,request));
		return reader;
	}

}

