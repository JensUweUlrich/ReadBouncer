/*
 * Data.cpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
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

    void Data::addActions()
    {
        DEBUGMESSAGE(std::cout, "start action request thread");
        // as long as signals are received from MinKnow
        // iterate over received data and stop further data allocation for every odd read on every odd channel
        while (isRunning())
        {
            DEBUGMESSAGE(std::cout, "setup new action request");
            GetLiveReadsRequest actionRequest{};
            GetLiveReadsRequest_Actions *actionList = actionRequest.mutable_actions();
            
            int counter = 0;
            while (counter < 2)
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
                    if (read.channelNr % 2 == 1 || read.readNr % 2 == 1)
                    {
                        GetLiveReadsRequest_Action *action = actionList->add_actions();
                        action->set_channel(read.channelNr);
                        GetLiveReadsRequest_StopFurtherData *data = action->mutable_stop_further_data();
                        *data = action->stop_further_data();
                        //DEBUGVAR(std::cout, action->has_stop_further_data());
                        action->set_number(read.readNr);
                        std::stringstream buf;
                        buf << "stop_further_" << read.channelNr << "_" << read.readNr;
                        action->set_action_id(buf.str());  
                        //DEBUGMESSAGE(std::cout, buf.str());
                        counter++;
                    }
                }
            }

            DEBUGMESSAGE(std::cout, "try to send action request");
            if (!stream->Write(actionRequest))
            {
                throw DataServiceException("Unable to add action to a live read stream!");
            }
            DEBUGMESSAGE(std::cout, "action request sent");

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
        setup->set_raw_data_type(GetLiveReadsRequest_RawDataType_CALIBRATED);
        setup->set_sample_minimum_chunk_size(50);

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
        std::set<std::string> unique_reads;
        while (stream->Read(&response))
        {
            runs = true;
            
        	for (GetLiveReadsResponse_ActionResponse actResponse : response.action_reponses())
            {
                actNr++;
                switch (actResponse.response())
                {
                    case GetLiveReadsResponse_ActionResponse_Response_SUCCESS:
                        DEBUGMESSAGE(std::cout, "Action " + actResponse.action_id() + " succeeded.");
                        //success++;
                        break;
                    case GetLiveReadsResponse_ActionResponse_Response_FAILED_READ_FINISHED:
                        DEBUGMESSAGE(std::cout, "Action " + actResponse.action_id() + " failed.");
                        break;
                }
            }
            
            
       		Map<uint32, GetLiveReadsResponse_ReadData> readData = response.channels();
       		for (MapPair<uint32, GetLiveReadsResponse_ReadData> entry : readData)
        	{
                auto ret = unique_reads.emplace(entry.second.id());
                // only add reads we did not already see to processing queue
                if (ret.second)
                {
           		    uint32 channel = entry.first;
           		    uint32 readNr = entry.second.number();
                    ReadCache r{};
                    r.channelNr = channel;
                    r.readNr = readNr;
                    mutex.lock();
                    reads.push(r);
                    mutex.unlock();
                }
       		}  
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
        sleep(1);
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

