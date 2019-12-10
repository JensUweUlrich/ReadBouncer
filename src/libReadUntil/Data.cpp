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
        DEBUGMESSAGE(std::cout, "start action request thread");
        // If unblock actions or stop further data
        GetLiveReadsRequest actionRequest;
        GetLiveReadsRequest_Actions actionList = actionRequest.actions();
        GetLiveReadsRequest_Action *action = actionList.add_actions();
        GetLiveReadsRequest_StopFurtherData data
        { };
        DEBUGMESSAGE(std::cout, "stop further");
        action->set_channel(1);
        action->set_allocated_stop_further_data(&data);
        action->set_number(1);
        action->set_action_id("stop_further_1_1");
        if (!stream->Write(actionRequest))
        {
            throw DataServiceException("Unable to add action to a live read stream!");
        }
        DEBUGMESSAGE(std::cout, "action written");
        //stream->WritesDone();
        DEBUGMESSAGE(std::cout, "done");
    }

    bool Data::isRunning()
    {
        return runs;
    }

    void Data::createSetupMessage()
    {
        DEBUGMESSAGE(std::cout, "start setup message thread");

        GetLiveReadsRequest setupRequest{};
        GetLiveReadsRequest_StreamSetup *setup = setupRequest.mutable_setup();
        setup->set_first_channel(1);
        setup->set_last_channel(512);
        setup->set_raw_data_type(GetLiveReadsRequest_RawDataType_CALIBRATED);
        setup->set_sample_minimum_chunk_size(50);

        DEBUGVAR(std::cout, setupRequest.has_setup());
        if (!stream->Write(setupRequest))
        {
            throw DataServiceException("Unable to setup a live read stream!");
        }

        DEBUGMESSAGE(std::cout, "leaving setup message thread");
    }

    void Data::getLiveSignals()
    {
        DEBUGMESSAGE(std::cout, "start getting signals thread");

        GetLiveReadsResponse response;
        uint32 channel{};
        uint32 readNr{};
        stream->Read(&response);
        //     while (stream->Read(&response))
        //     {
        DEBUGVAR(std::cout, response.seconds_since_start());
        /*          for (GetLiveReadsResponse_ActionResponse actResponse : response.action_reponses())
        {
            switch (actResponse.response())
            {
                case GetLiveReadsResponse_ActionResponse_Response_SUCCESS:
                    DEBUGMESSAGE(std::cout, "Action " + actResponse.action_id() + " was successful");
                case GetLiveReadsResponse_ActionResponse_Response_FAILED_READ_FINISHED:
                    DEBUGMESSAGE(std::cout, "Action " + actResponse.action_id() + " failed.");
            }
        }
*/
        Map<uint32, GetLiveReadsResponse_ReadData> readData = response.channels();
        for (MapPair<uint32, GetLiveReadsResponse_ReadData> entry : readData)
        {
            channel = entry.first;
            readNr = entry.second.number();
            break;
        }

        //  break;
        //  }
        DEBUGVAR(std::cout, channel);
        DEBUGVAR(std::cout, readNr);

        // grpc::Status status = stream->Finish();
    }

    void Data::getLiveReads()
    {

        grpc::Status status;

        stream = stub->get_live_reads(&context);
        // first write setup the stream

        std::thread setupThread(&Data::createSetupMessage, this);

        setupThread.join();
        DEBUGMESSAGE(std::cout, "written");


        //grpc::ClientContext context2;
        //stream = stub->get_live_reads(&context);

        //std::thread actionThread(&Data::addAction, this);
        //actionThread.join();

        std::thread readerThread(&Data::getLiveSignals, this);


        readerThread.join();

        /*stream->WritesDone();
        status = stream->Finish();

        DEBUGVAR(std::cout, status.error_code());
        DEBUGVAR(std::cout, status.error_message());
        DEBUGVAR(std::cout, status.error_details());
*/
    }

}

