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
	Manager::Manager(std::shared_ptr<::grpc::Channel> channel, bool secure_connect)
	{
		stub = ManagerService::NewStub(channel);
		this->secure_connect = secure_connect;
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
		
		//std::cout << cells.size() << 'n';

		grpc::Status status = reader->Finish();

		
		if(!status.ok())
		{
			
			throw ReadUntilClientException(status.error_message());
		}
		return cells;
	}

	std::string Manager::getFlowCellName(FlowCellPosition &dev)
	{
		return dev.name();
	}

	uint32_t Manager::getRpcPort(FlowCellPosition &dev)
	{
		if (this->secure_connect)
			return dev.rpc_ports().secure();
		//else // return no secure as minknow 5.0.0 has no insecure connection! 
			//return dev.rpc_ports().insecure();
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
		//GetVersionInfoResponse response; // minknow_api 4.5.0
		minknow_api::instance::GetVersionInfoResponse response;
		grpc::ClientContext c;
		::grpc::Status status = stub->get_version_info(&c, request, &response);
		if (status.ok())
		{
			return response.guppy_connected_version();
		}
		else
		{
			throw ReadUntilClientException(status.error_message());
		}
	}

	std::string Manager::getTokenFilePath()
	{
		minknow_api::manager::LocalAuthenticationTokenPathRequest request;
		minknow_api::manager::LocalAuthenticationTokenPathResponse response;
		grpc::ClientContext c;
		::grpc::Status status = stub->local_authentication_token_path(&c, request, &response);

		return response.path();
	}

}
