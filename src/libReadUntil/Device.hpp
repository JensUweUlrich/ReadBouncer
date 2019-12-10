/*
 * Device.hpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
 */

#include <string>
#include <grpcpp/grpcpp.h>
#include <minknow/rpc/device.grpc.pb.h>

#include "MinKnowService.hpp"
#include "DeviceServiceException.hpp"
#include "../debug_messages.hpp"

using namespace ::ont::rpc::device;

#ifndef LIBREADUNTIL_DEVICE_HPP_
#define LIBREADUNTIL_DEVICE_HPP_

namespace readuntil
{

	class Device: public MinKnowService
	{
		private:
			std::unique_ptr<DeviceService::Stub> stub;
			GetFlowCellInfoResponse *flowCellInfo;
			GetDeviceInfoResponse *deviceInfo;
			GetDeviceStateResponse *deviceState;
			bool deviceInfoLoaded = false;
			bool deviceStateLoaded = false;
			bool flowCellInfoLoaded = false;

			inline void getDeviceInfo()
			{
				GetDeviceInfoRequest request;
				::grpc::Status status = stub->get_device_info(&context, request, *&deviceInfo);
				if (!status.ok())
				{
					throw DeviceServiceException("Could not get Device information : " + status.error_message());
				}
				deviceInfoLoaded = true;
			}

			inline void getDeviceState()
			{
				GetDeviceStateRequest request;
				::grpc::Status status = stub->get_device_state(&context, request, *&deviceState);
				if (!status.ok())
				{
					throw DeviceServiceException("Could not get Device state : " + status.error_message());
				}
				deviceStateLoaded = true;
			}

			inline void getFlowcellInfo()
			{
				GetFlowCellInfoRequest request;
				::grpc::Status status = stub->get_flow_cell_info(&context, request, *&flowCellInfo);
				if (!status.ok())
				{
					throw DeviceServiceException("Could not get Flowcell information : " + status.error_message());
				}
				flowCellInfoLoaded = true;
			}

		public:
			Device(std::shared_ptr<::grpc::Channel> channel);
			~Device()
			{
				stub.release();
			}
			std::string getDeviceId();
			std::string getDeviceType();
			bool isReady();
			bool isDisconnected();
			bool hasFlongleAdapter();
			bool hasFlowCell();
			std::unique_ptr<grpc::ClientReader<GetDeviceStateResponse>> streamDeviceState();
	}
	;

} //namespace

#endif /* LIBREADUNTIL_DEVICE_HPP_ */
