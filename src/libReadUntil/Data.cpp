/*
 * Data.cpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
 */

#include "Data.hpp"

namespace readuntil
{
	Data::Data(std::shared_ptr<::grpc::Channel> channel)
	{
		stub = DataService::NewStub(channel);
	}

	void Data::addAction()
	{
		// If unblock actions or stop further data
		GetLiveReadsRequest_Action *action = actionList.add_actions();
		GetLiveReadsRequest_StopFurtherData data
		{ };
		action->set_channel(1);
		action->set_allocated_stop_further_data(&data);
		action->set_number(1);
		action->set_action_id("stop_further_1_1");
	}

	void Data::getLiveReads()
	{

		// setup parameters of the live reads stream
		GetLiveReadsRequest_StreamSetup setup;
		setup.set_first_channel(1);
		setup.set_last_channel(1);
		setup.set_raw_data_type(GetLiveReadsRequest_RawDataType_CALIBRATED);
		setup.set_sample_minimum_chunk_size(0);

		//
		GetLiveReadsRequest setupRequest;
		setupRequest.set_allocated_setup(&setup);
		std::shared_ptr<grpc::ClientReaderWriter<GetLiveReadsRequest, GetLiveReadsResponse>> stream(stub->get_live_reads(&context));
		// first write setup the stream
		if (!stream->Write(setupRequest))
		{
			throw DataServiceException("Unable to setup a live read stream!");
		}

		// If unblock actions or stop further data

		if (actionList.GetCachedSize() > 0)
		{
			GetLiveReadsRequest actionRequest;
			actionRequest.set_allocated_actions(&actionList);

			// write actions to the stream
			if (!stream->Write(actionRequest))
			{
				throw DataServiceException("Unable to write action request to the live read stream!");
			}
		}

		stream->WritesDone();


		GetLiveReadsResponse response;
		while(stream -> Read(&response))
		{
			for (GetLiveReadsResponse_ActionResponse actResponse : response.action_reponses())
			{
				switch (actResponse.response())
				{
					case GetLiveReadsResponse_ActionResponse_Response_SUCCESS:
						DEBUGMESSAGE(std::cout, "Action " + actResponse.action_id() + " was successful");
					case GetLiveReadsResponse_ActionResponse_Response_FAILED_READ_FINISHED:
						DEBUGMESSAGE(std::cout, "Action " + actResponse.action_id() + " failed.");
				}
			}

			Map<uint32, GetLiveReadsResponse_ReadData> readData = response.channels();
			for (MapPair<uint32, GetLiveReadsResponse_ReadData> entry : readData)
			{
				DEBUGMESSAGE(std::cout, "Channel Number: " + entry.first);
				DEBUGMESSAGE(std::cout, "Read Number: " + entry.second.number());
				DEBUGMESSAGE(std::cout, "Chunk Length: " + entry.second.chunk_length());
			}

		}

		grpc::Status status = stream->Finish();
	}

}

