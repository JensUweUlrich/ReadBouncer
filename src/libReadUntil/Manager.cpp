/*
 * Manager.cpp
 *
 * Created on: 12.11.2019
 * Author: Jens-Uwe Ulrich
 *
 * Fetching information about devices connected to MinKNOW via Remote Procedure Calls
 *
 * [update] 27.08.2020 resolveRpcPort for a given device name
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

	std::vector<FlowCellPosition> Manager::getFlowCells()
	{
		FlowCellPositionsRequest request;
		FlowCellPositionsResponse response;
		std::unique_ptr<grpc::ClientReader<FlowCellPositionsResponse>> reader(stub->flow_cell_positions(&context, request));
		std::vector<FlowCellPosition> cells{ };
		while (reader->Read(&response))
		{
			cells.assign(response.positions().begin(), response.positions().end());

		}
		grpc::Status status = reader->Finish();
		if(!status.ok())
		{
            //throw ReadUntilClientException(status.error_message());
            ReadUntilClientException(status.error_message());
		}
		return cells;
	}

	std::string Manager::getFlowCellName(FlowCellPosition &dev)
	{
		return dev.name();
	}

	uint32_t Manager::getRpcPort(FlowCellPosition &dev)
	{
		return dev.rpc_ports().insecure();
	}

	uint32_t Manager::resolveRpcPort(std::string &deviceName)
	{
		for(FlowCellPosition fp : getFlowCells())
		{
			if(deviceName == getFlowCellName(fp))
			{
				return getRpcPort(fp);
			}
		}
		return 0;
	}

	std::string Manager::getGuppyVersion()
	{
		GetVersionInfoRequest request;
		GetVersionInfoResponse response;
		::grpc::Status status = stub->get_version_info(&context, request, &response);
		if (status.ok())
		{
			return response.guppy_connected_version();
		}
		else
		{
            ReadUntilClientException(status.error_message());
            //throw ReadUntilClientException(status.error_message());
		}
	}

}
