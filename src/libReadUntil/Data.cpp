/*
 * Data.cpp
 *
 *  Created on: 12.11.2019
 *      Author: jens-uwe.ulrich
 */

#include "Data.hpp"

namespace readuntil
{

    void processSignals()
    {
		
    }

    Data::Data(std::shared_ptr<::grpc::Channel> channel)
    {
        stub = DataService::NewStub(channel);
    }

    void Data::addUnblockAction(GetLiveReadsRequest_Actions *actionList, ReadCache &read)
    {
        GetLiveReadsRequest_Action *action = actionList->add_actions();
        action->set_channel(read.channelNr);
        GetLiveReadsRequest_UnblockAction *data = action->mutable_unblock();
        *data = action->unblock();
        data->set_duration(0.1);
        action->set_number(read.readNr);
        std::stringstream buf;
        buf << "unblock_" << read.channelNr << "_" << read.readNr;
        action->set_action_id(buf.str());
    }

    void Data::addStopReceivingDataAction(GetLiveReadsRequest_Actions *actionList, ReadCache &read)
    {
        GetLiveReadsRequest_Action *action = actionList->add_actions();
        action->set_channel(read.channelNr);
        GetLiveReadsRequest_StopFurtherData *data = action->mutable_stop_further_data();
        *data = action->stop_further_data();
        action->set_number(read.readNr);
        std::stringstream buf;
        buf << "stop_receiving_" << read.channelNr << "_" << read.readNr;
        action->set_action_id(buf.str());
    }

    void Data::addActions()
    {
        DEBUGMESSAGE(std::cout, "start action request thread");
        // as long as signals are received from MinKnow
        // iterate over received data and stop further data allocation for every odd read on every odd channel
        while (isRunning())
        {
            std::chrono::time_point<std::chrono::system_clock> start, end;
            start = std::chrono::system_clock::now();
            //DEBUGMESSAGE(std::cout, "setup new action request");
            GetLiveReadsRequest actionRequest{};
            GetLiveReadsRequest_Actions *actionList = actionRequest.mutable_actions();
            
            int counter = 0;
            while (counter < actionBatchSize)
            {
                ReadCache read{};
                bool hasElement = false;
                mutex.lock();
                // take read out of the queue if it's not empty
                if (!reads.empty())
                {
                    read.channelNr = reads.front().channelNr;
                    read.readNr = reads.front().readNr;
                    hasElement = true;
                    reads.pop();
                }
                mutex.unlock();

                if (hasElement)
                {
                    if (read.channelNr % unblockChannels != 0 || read.readNr % unblockReads != 0)
                    {
                        addUnblockAction(actionList, read);
                        addStopReceivingDataAction(actionList, read);
                        counter+=2;
                    }
                }
            }

            //DEBUGMESSAGE(std::cout, "try to send action request");
            if (!stream->Write(actionRequest))
            {
                throw DataServiceException("Unable to add action to a live read stream!");
            }
            //DEBUGMESSAGE(std::cout, "action request sent");

            end = std::chrono::system_clock::now();
            int elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            if(elapsed_milliseconds < 100)
            {
                sleep(100 - elapsed_milliseconds);
            }
        }
         DEBUGMESSAGE(std::cout, "leaving action request thread");
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
        setup->set_raw_data_type(GetLiveReadsRequest_RawDataType_UNCALIBRATED);
        setup->set_sample_minimum_chunk_size(0);

        DEBUGVAR(std::cout, setupRequest.has_setup());
        if (!stream->Write(setupRequest))
        {
            throw DataServiceException("Unable to setup a live read stream!");
        }

        DEBUGMESSAGE(std::cout, "leaving setup message thread");
    }

	/*Data::getSignalType()
	{
		GetDataTypesRequest request;
		GetDataTypesResponse response;
		::grpc::ClientContext context;
		::grpc::Status status = stub->get_data_types(&context, request, &response);
		if (status.ok())
		{
			GetDataTypesResponse_DataType dataType = response.calibrated_signal();
			
		}
		else
		{
			throw DataServiceException("Unable to resolve signal data type!");
		}
	}
*/
    void Data::getLiveSignals()
    {
        DEBUGMESSAGE(std::cout, "start getting signals thread");
        GetLiveReadsResponse response;
        int actNr = 0;
        int success = 0;
        int failed = 0;
        
        while (stream->Read(&response))
        {
            runs = true;
            
        	for (GetLiveReadsResponse_ActionResponse actResponse : response.action_reponses())
            {
                actNr++;
                switch (actResponse.response())
                {
                    case GetLiveReadsResponse_ActionResponse_Response_SUCCESS:
                        //DEBUGMESSAGE(std::cout, "Action " + actResponse.action_id() + " succeeded.");
                        success++;
                        break;
                    case GetLiveReadsResponse_ActionResponse_Response_FAILED_READ_FINISHED:
                        //DEBUGMESSAGE(std::cout, "Action " + actResponse.action_id() + " failed.");
                        failed++;
                        break;
                }
            }
            std::stringstream ss;
            ss << "Success/Failed rate = " << success <<"/" << failed; 
            DEBUGMESSAGE(std::cout, ss.str());
            
       		Map<uint32, GetLiveReadsResponse_ReadData> readData = response.channels();
       		for (MapPair<uint32, GetLiveReadsResponse_ReadData> entry : readData)
        	{
                // only add reads we did not already see to processing queue
                if (std::find(uniqueReadIds.begin(), uniqueReadIds.end(), entry.second.id()) == uniqueReadIds.end())
                {
           		    uint32 channel = entry.first;
           		    uint32 readNr = entry.second.number();
                    ReadCache r{};
                    r.channelNr = channel;
                    r.readNr = readNr;
                    mutex.lock();
                    reads.push(r);
                    mutex.unlock();
                    uniqueReadIds.push_back(entry.second.id());
                }
       		}
            std::stringstream ss2;
            ss2 << "ReadCacheSize : " << reads.size();
            DEBUGMESSAGE(std::cout, ss2.str());  
       }
       DEBUGMESSAGE(std::cout, "leaving signals thread");
       runs = false;
    }

    void Data::getLiveReads()
    {
        runs = true;
        grpc::Status status;

        stream = stub->get_live_reads(&context);
        // first write setup the stream

        std::thread setupThread(&Data::createSetupMessage, this);

        setupThread.join();
        DEBUGMESSAGE(std::cout, "written");

        std::thread readerThread(&Data::getLiveSignals, this);
        std::thread actionThread(&Data::addActions, this);


        readerThread.join();
        actionThread.join();
	
        stream->WritesDone();
	    context.TryCancel();
        status = stream->Finish();

        DEBUGVAR(std::cout, status.error_code());
        DEBUGVAR(std::cout, status.error_message());
        DEBUGVAR(std::cout, status.error_details());

    }

}

