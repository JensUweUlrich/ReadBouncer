/*
 * Data.cpp
 *
 *  Created on: 12.11.2019
 *      Author: jens-uwe.ulrich
 */

#include "Data.hpp"

namespace readuntil
{

    /**
    *   Constructor of Data class
    *   @channel:   shared pointer to already opened GRPC channel for communication with MinKNOW
    */
    Data::Data(std::shared_ptr<::grpc::Channel> channel)
    {
        stub = DataService::NewStub(channel);
        acq = new readuntil::Acquisition(channel);
        conf = new readuntil::AnalysisConfiguration(channel);
        // TODO: better usage of logging
        data_logger = spdlog::get("RUClientLog");
        resolveFilterClasses();
        
    }

    /**
        load all possible read classifications
        only reads classified as strand or adapter should be used
        adds corresponding integer code of the two classifications to the filterClasses set
        Possible read classifications by MinKNOW:
        83: "strand",
        67: "strand1",
        77: "multiple",
        90: "zero",
        65: "adapter",
        66: "mux_uncertain",
        70: "user2",
        68: "user1",
        69: "event",
        80: "pore",
        85: "unavailable",
        84: "transition",
        78: "unclassed",
    */
    void Data::resolveFilterClasses()
    {
        Map<int32, string> classes = conf->getReadClassifications();
		for (MapPair<int32, string> p : classes)
		{
            // we are currently only interested in reads classified as 'strand' or 'adapter'
			if (p.second=="strand" || p.second=="adapter")
            {
                filterClasses.insert(p.first);
            }
		}
    }

    /**
    *   creates an unblock action for a given read and adds the action to the action list
    *   @actionList         : List of action messages
    *   @response           : action response object storing infos like read number and channel number
    *   @unblock_duration   : unblock duration in seconds
    */
    void Data::addUnblockAction(GetLiveReadsRequest_Actions* actionList, ActionResponse& response, const double unblock_duration)
    {
        // create an action message and add it to the action list
        GetLiveReadsRequest_Action* action = actionList->add_actions();
        // FlowCell channel number on which we call the action
        action->set_channel(response.channelNr);

        // define action as unblock action
        GetLiveReadsRequest_UnblockAction* data = action->mutable_unblock();
        *data = action->unblock();

        // duration in seconds for which we reverse the current on that channel to unblock the pore
        data->set_duration(unblock_duration);

        // specify read number that shall be unblocked -> avoids unblock of consecutive read
        action->set_number(response.readNr);

        // create unique action id by using a timestamp
        std::stringstream buf;
        std::chrono::milliseconds ms = std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch());
        buf << "unblock_" << response.channelNr << "_" << response.readNr << "_" << ms.count();
        action->set_action_id(buf.str());
    }

    /**
    *   creates an action message for not receiving more data for tha read
    *   adds the action to the action list
    *   @actionList : List of action messages
    *   @response   : action response object storing infos like read number and channel number
    */
    void Data::addStopReceivingDataAction(GetLiveReadsRequest_Actions* actionList, ActionResponse& response)
    {
        // create an action message and add it to the action list
        GetLiveReadsRequest_Action* action = actionList->add_actions();
        // FlowCell channel number on which we call the action
        action->set_channel(response.channelNr);

        // define action as stop_further_data action
        GetLiveReadsRequest_StopFurtherData* data = action->mutable_stop_further_data();
        *data = action->stop_further_data();

        // specify read number for which we don't want to receive more data
        // but only for this read, for the next one we want to get live signals
        action->set_number(response.readNr);

        // create unique action id by using a timestamp
        std::stringstream buf;
        buf << "stop_receiving_" << response.channelNr << "_" << response.readNr;
        action->set_action_id(buf.str());
    }

    /**
    *   converts a given string of signals into a vector of float signals
    *   @signalString   : string of nanopore signals
        @return         : vector of float nanopore signals
    */
    std::vector<float> Data::string_to_float(std::string const & signalString)
    {
        assert(signalString.size() % sizeof(float) == 0);

        std::vector<float> result(signalString.size() / sizeof(float));

        if (!result.empty())
        {
            std::copy(signalString.data(), signalString.data() + signalString.size(),
                  reinterpret_cast<char *>(&result.front()));
        }

        return result;
    }

    /**
    *   take read action responses from the queue and add unblock and/or stop_receiving_data
    *   action to the action list and write action request to the bidirectional stream
    *   @action_queue   : safe queue with reads for which action messages shall be sent to MinKNOW
    */
    void Data::sendActions(SafeQueue<readuntil::ActionResponse>& action_queue, SafeQueue<Durations>& duration_queue)
    {
        data_logger->debug("start action request thread");
        // as long as signals are received from MinKnow
        // iterate over received data and stop further data allocation for every odd read on every odd channel
        
        // only send action messages while sequencing is still ongoing
        while (isRunning())
        {
            // we want take the time between creation of first and last action message
            std::chrono::time_point<std::chrono::system_clock> start, end;
            start = std::chrono::system_clock::now();

            // create an action request that will be written on the stream later
            GetLiveReadsRequest actionRequest{};
            // list of actions is attribute of an action request object
            GetLiveReadsRequest_Actions* actionList = actionRequest.mutable_actions();
            
            int counter = 0;
            // only actionBatchSize is the number of action messages that are send in one request
            while (counter < actionBatchSize)
            {
                if (!action_queue.empty())
                {
                    // take read from the queue and add unblock message for that read to the queue
                    // if read is classified as host read
                    ActionResponse readResponse = action_queue.pop();
                    if (readResponse.response)
                    {
                        addUnblockAction(actionList, readResponse, 0.1);
                        counter++;
                        
                    }

                    // add stop_receiving_data message for every read in the queue
                    addStopReceivingDataAction(actionList, readResponse);
                    counter++;
                    readResponse.processingTimes.timeCompleteRead.stop();
                    if (readResponse.response)
                    {
                        duration_queue.push(Durations{
                            readResponse.processingTimes.timeCompleteRead.elapsed(),
                            -1,
                            readResponse.processingTimes.timeBasecallRead.elapsed(),
                            readResponse.processingTimes.timeClassifyRead.elapsed(),
                            });
                    }
                    else
                    {
                        duration_queue.push(Durations{
                            -1,
                            readResponse.processingTimes.timeCompleteRead.elapsed(),
                            readResponse.processingTimes.timeBasecallRead.elapsed(),
                            readResponse.processingTimes.timeClassifyRead.elapsed(),
                            });
                    }
                    
                }
                else
                    break;
            }

            // increase number of action messages in one request if queue is too full
            // decrease otherwise
            adaptActionBatchSize(action_queue.size());

            // write action request to the stream -> send message to MinKNOW
            // throw exception if that was nut successful
            if (!stream->Write(actionRequest))
            {
                // TODO : throw more specific exception
                throw DataServiceException("Unable to add action to a live read stream!");
            }

            // there should be at least 100 ms between two separate action requests
            // otherwise MinKNOW action request queue is overflooded and drops requests
            end = std::chrono::system_clock::now();
            int elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            if(elapsed_milliseconds < 100)
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(100 - elapsed_milliseconds));
            }
            
        }
        data_logger->debug("leaving action request thread");
    }
    
    /**
    *   automatically adapt the number of action messages sent in one action request via the stream
    *   @queue_size : size of the action response queue
    * 
    *   TODO: evaluate if porportion of increase and decrease in size is ok
    */
    void Data::adaptActionBatchSize(const int queue_size)
    {
        if(queue_size > 0)
        {
            actionBatchSize += queue_size;
        }
        else
        {
            actionBatchSize *= 0.8; 
        }
    }

    /**
    *   checks whether there are still read signals arriving from MinKNOW
    *   @return: true, if sequencig is ongoing and false otherwise
    */
    bool Data::isRunning()
    {
        return runs;
    }

    /**
    *   starts the bidirectional stream with MinKNOW
    *   sends setup message to MinKNOW, which is needed initially receiving data
    *   @throws : DataServiceException
    */
    void Data::startLiveStream()
    {
        data_logger->debug("start setup message thread");

        runs = true;

        // start streaming live nanopore signals
        stream = stub->get_live_reads(&context);

        // create setup messae
        GetLiveReadsRequest setupRequest{};
        GetLiveReadsRequest_StreamSetup* setup = setupRequest.mutable_setup();

        // we want to receive signals from all 512 channels of the MinION
        // has to be changed in case Flongle or PromethION is used
        // TODO: set last channel based on device type
        setup->set_first_channel(1);
        setup->set_last_channel(512);

        // we only want to receive calibrated data
        setup->set_raw_data_type(GetLiveReadsRequest_RawDataType_CALIBRATED);
        // minimum number of signals for part of a read to be sent via the stream
        // zero means no limitation
        // TODO: Try to find out if other parameter like e.g. 4 is better
        setup->set_sample_minimum_chunk_size(0);

        // write setup message to stream and throw exception if that was not successful
        if (!stream->Write(setupRequest))
        {
            // TODO: create more specific exception
            throw DataServiceException("Unable to setup a live read stream!");
        }

        data_logger->debug("leaving setup message thread");
    }

	/**
    *   pull live nanopore signals from the stream and add the reads to the basecalling queue
    *   @basecall_queue : safe queue for storing reads ready for basecalling
    *   
    *   TODO: analyze action responses and log success and failed actions
    */
    void Data::getLiveSignals(SafeQueue<SignalRead>& basecall_queue)
    {
        data_logger->debug("start getting signals thread");
        GetLiveReadsResponse response;

        // as long as there is incoming data on the stream we read it and store
        // temporary data in response variable
        while (stream->Read(&response))
        {
            // stop trying to read from the stream if sequencing has been finished
            if (acq->isFinished())
            {
                runs = false;
                break;
            }

            // iterate over the read data from each channel
            Map<uint32, GetLiveReadsResponse_ReadData> readData = response.channels();
       		for (MapPair<uint32, GetLiveReadsResponse_ReadData> entry : readData)
        	{
                TimeMeasures times;
                times.timeCompleteRead.start();
                // only process read chunks from filter classes (e.g. strand or adapter)
                bool filtered = true;
                for (int32 cl : entry.second.chunk_classifications())
                {
                    if (filterClasses.find(cl) != filterClasses.end())
                    {
                        filtered = false;
                        break;
                    }
                }

                if (filtered)
                {
                    continue;
                }
                // add read to basecalling queue
                // nanopore signal string is converted into a vector of floating point numbers
                basecall_queue.push(SignalRead{ entry.first,
                                                    entry.second.number(),
                                                    entry.second.id(),
                                                    string_to_float(entry.second.raw_data()),
                                                    times});
       		}

       }
       runs = false;
    }

    /**
    *   stop streaming data with a clean finish
    */
    void Data::stopLiveStream()
    {
        runs = false;
        stream->WritesDone();
	    context.TryCancel();
        grpc::Status status = stream->Finish();

    }
    
}

