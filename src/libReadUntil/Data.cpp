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
        acq = new readuntil::Acquisition(channel);
        data_logger = spdlog::get("RUClientLog");
    }

    void Data::addUnblockAction(GetLiveReadsRequest_Actions *actionList, ReadCache &read, const double unblock_duration)
    {
        GetLiveReadsRequest_Action *action = actionList->add_actions();
        action->set_channel(read.channelNr);
        GetLiveReadsRequest_UnblockAction *data = action->mutable_unblock();
        *data = action->unblock();
        data->set_duration(unblock_duration);
        action->set_number(read.readNr);
        std::stringstream buf;
        buf << read.channelNr << "_" << read.readNr;
        responseMutex.lock();
        std::map<string, ReadResponse>::iterator it = responseCache.find(buf.str());
        if (it != responseCache.end())
        {
            (*it).second.unblock_duration = unblock_duration;
        }
        responseMutex.unlock();
        buf << "unblock_" << buf.str();
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
        data_logger->debug("start action request thread");
        // as long as signals are received from MinKnow
        // iterate over received data and stop further data allocation for every odd read on every odd channel
        while (isRunning())
        {
            std::chrono::time_point<std::chrono::system_clock> start, end;
            start = std::chrono::system_clock::now();
            GetLiveReadsRequest actionRequest{};
            GetLiveReadsRequest_Actions *actionList = actionRequest.mutable_actions();
            
            int counter = 0;
            while (counter < actionBatchSize)
            {
                ReadCache read{};
                bool hasElement = false;
                readMutex.lock();
                // take read out of the queue if it's not empty
                if (!reads.empty())
                {
                    read.channelNr = reads.front().channelNr;
                    read.readNr = reads.front().readNr;
                    hasElement = true;
                    reads.pop();
                }
                readMutex.unlock();

                if (hasElement)
                {
                    // 
                    switch(read.channelNr % 4)
                    {
                        case 0:
                        {
                            //channel numbers with mod 4 == 0 are sequenced as usual
                            std::stringstream buf;
                            buf << read.channelNr << "_" << read.readNr;
                            responseMutex.lock();
                            std::map<string, ReadResponse>::iterator it = responseCache.find(buf.str());
                            if (it != responseCache.end())
                            {
                                (*it).second.unblock_duration = -1.0;
                            }
                            responseMutex.unlock();
                            break;
                        }
                        case 1:
                        {
                            // unblock odd numbered reads with duration 1 sec 
                            if (read.readNr % 2 == 1)
                            {
                                addUnblockAction(actionList, read, 1.0);
                                counter++;
                            }
                            break;
                        }
                        case 2:
                        {
                            // unblock odd numbered reads with duration 0.1 sec
                            if (read.readNr % 2 == 1)
                            {
                                addUnblockAction(actionList, read, 0.1);
                                counter++;
                            }
                            break;
                        }
                        case 3:
                        {
                            // unblock every fourth read with duration 0.1 sec
                            if (read.readNr % 4 == 0)
                            {
                                addUnblockAction(actionList, read, 0.1);
                                counter++;
                            }
                            break;
                        }
                        default:
                            // do nothing
                            break;
                    }
                }
            }

            if (!stream->Write(actionRequest))
            {
                throw DataServiceException("Unable to add action to a live read stream!");
            }

            end = std::chrono::system_clock::now();
            int elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            if(elapsed_milliseconds < 100)
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(100 - elapsed_milliseconds));
            }
        }
        data_logger->debug("leaving action request thread");
    }

    bool Data::isRunning()
    {
        return runs;
    }

    void Data::createSetupMessage()
    {
        data_logger->debug("start setup message thread");

        GetLiveReadsRequest setupRequest{};
        GetLiveReadsRequest_StreamSetup *setup = setupRequest.mutable_setup();
        setup->set_first_channel(1);
        setup->set_last_channel(512);
        setup->set_raw_data_type(GetLiveReadsRequest_RawDataType_UNCALIBRATED);
        setup->set_sample_minimum_chunk_size(4);

        if (!stream->Write(setupRequest))
        {
            throw DataServiceException("Unable to setup a live read stream!");
        }

        data_logger->debug("leaving setup message thread");
    }

    void Data::printResponseData()
    {
        std::ofstream ofs("readMetaData.csv", std::ofstream::out);
        ofs << "read_id\tchannel_nr\tread_nr\tresponse\tunblock_duration\tsamples_since_start\tseconds_since_start\tstart_sample\tchunk_start_sample\tchunk_length\n"; 
        while (isRunning())
        {
            ReadResponse r{};
            responseMutex.lock();
            if (!responseCache.empty())
            {
                std::map<string,ReadResponse>::iterator it = responseCache.begin();
                if((*it).second.response > 0 || (*it).second.unblock_duration < 0)
                {
                    ofs << (*it).first << "\t" << (*it).second.channelNr << "\t" << (*it).second.readNr << "\t" << (*it).second.response << "\t" << (*it).second.unblock_duration << "\t" << (*it).second.samples_since_start << "\t" << (*it).second.start_sample << "\t" << (*it).second.chunk_start_sample << "\t" << (*it).second.chunk_length << "\n";
                    responseCache.erase(it);
                }
            }
            responseMutex.unlock();

        }

        responseMutex.lock();
        for (std::map<string,ReadResponse>::iterator it = responseCache.begin(); it != responseCache.end(); ++it)
        {
            ofs << (*it).first << "\t" << (*it).second.channelNr << "\t" << (*it).second.readNr << "\t" << (*it).second.response << "\t" << (*it).second.unblock_duration << "\t" << (*it).second.samples_since_start << "\t" << (*it).second.start_sample << "\t" << (*it).second.chunk_start_sample << "\t" << (*it).second.chunk_length << "\n";
        }
        responseMutex.unlock();
        ofs.close();
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
        data_logger->debug("start getting signals thread");
        GetLiveReadsResponse response;
        int actNr = 0;
        int success = 0;
        int failed = 0;

        while (stream->Read(&response))
        {
            if (acq->isFinished())
            {
                runs = false;
                break;
            }

            runs = true;
            
        	for (GetLiveReadsResponse_ActionResponse actResponse : response.action_reponses())
            {
                // add action success/failed information to responseCache entry
                std::string id(actResponse.action_id());
                if (id.rfind("unblock_",0) == 0)
                {
                    id = id.substr(8);
                }
                responseMutex.lock();
                std::map<string, ReadResponse>::iterator it = responseCache.find(id);
                if (it != responseCache.end())
                {
                    switch(actResponse.response())
                    {
                        case GetLiveReadsResponse_ActionResponse_Response_SUCCESS:
                            (*it).second.response = 1;
                            break;
                        case GetLiveReadsResponse_ActionResponse_Response_FAILED_READ_FINISHED:
                            (*it).second.response = 2;
                            break;
                    }
                }
                responseMutex.unlock();


                actNr++;
                switch (actResponse.response())
                {
                    case GetLiveReadsResponse_ActionResponse_Response_SUCCESS:
                        success++;
                        break;
                    case GetLiveReadsResponse_ActionResponse_Response_FAILED_READ_FINISHED:
                        failed++;
                        break;
                }
            }
            std::stringstream ss;
            ss << "Success/Failed rate = " << success <<"/" << failed; 
            data_logger->info(ss.str());
            
       		Map<uint32, GetLiveReadsResponse_ReadData> readData = response.channels();
       		for (MapPair<uint32, GetLiveReadsResponse_ReadData> entry : readData)
        	{
                // only add reads we did not already see to processing queue
                if (std::find(uniqueReadIds.begin(), uniqueReadIds.end(), entry.second.id()) == uniqueReadIds.end())
                {

                    // store read data for csv printing
                    ReadResponse respData{};
                    respData.channelNr = entry.first;
                    respData.readNr = entry.second.number();
                    respData.samples_since_start = response.samples_since_start();
                    respData.seconds_since_start = response.seconds_since_start();
                    respData.start_sample = entry.second.start_sample();
                    respData.chunk_start_sample = entry.second.chunk_start_sample();
                    respData.chunk_length = entry.second.chunk_length();
                    std::stringstream buf;
                    buf << respData.channelNr << "_" << respData.readNr;
                    responseMutex.lock();
                    responseCache.emplace(buf.str(), respData);
                    responseMutex.unlock();

                    // add read to read cache for action request processing
           		    uint32 channel = entry.first;
           		    uint32 readNr = entry.second.number();
                    ReadCache r{};
                    r.channelNr = channel;
                    r.readNr = readNr;
                    readMutex.lock();
                    reads.push(r);
                    readMutex.unlock();
                    uniqueReadIds.push_back(entry.second.id());
                }
       		}
            std::stringstream ss2;
            ss2 << "ReadCacheSize : " << reads.size();
            data_logger->debug(ss2.str());  
       }
       data_logger->debug("leaving signals thread");
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

        std::thread readerThread(&Data::getLiveSignals, this);
        std::thread actionThread(&Data::addActions, this);
        std::thread printThread(&Data::printResponseData, this);


        readerThread.join();
        actionThread.join();
	
        stream->WritesDone();
	    context.TryCancel();
        status = stream->Finish();

        printThread.join();

        data_logger->debug(status.error_code());
        data_logger->debug(status.error_message());
        data_logger->debug(status.error_details());

    }

}

