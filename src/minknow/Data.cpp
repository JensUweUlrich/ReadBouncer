/*
 * Data.cpp
 *
 *  Created on: 12.11.2019
 *      Author: jens-uwe.ulrich
 */

#include "Data.hpp"

namespace readuntil
{
    std::filesystem::path CSVFile{};

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

	std::random_device rd;
	auto seed_data = std::array<int, std::mt19937::state_size> {};
	std::generate(std::begin(seed_data), std::end(seed_data), std::ref(rd));
	std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
	std::mt19937 generator(seq);
	uuid_generator = uuids::uuid_random_generator(generator);
        
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
    void Data::addUnblockAction(GetLiveReadsRequest_Actions* actionList, uint32_t channelNr, uint32_t readNr , const double unblock_duration)
    {
        // create an action message and add it to the action list
        GetLiveReadsRequest_Action* action = actionList->add_actions();
        // FlowCell channel number on which we call the action
        action->set_channel(channelNr);

        // define action as unblock action
        GetLiveReadsRequest_UnblockAction* data = action->mutable_unblock();
        *data = action->unblock();

        // duration in seconds for which we reverse the current on that channel to unblock the pore
        data->set_duration(unblock_duration);

        // specify read number that shall be unblocked -> avoids unblock of consecutive read
        action->set_number(readNr);

        // create uuid for action
        
        uuids::uuid const id = uuid_generator();
        
        action->set_action_id(uuids::to_string(id));
        //std::cout << uuids::to_string(id) << std::endl;
        /*std::stringstream buf;
        std::chrono::milliseconds ms = std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch());
        buf << "unblock_" << response.channelNr << "_" << response.readNr << "_" << ms.count();
        action->set_action_id(buf.str());
        */

    }

    /**
    *   creates an action message for not receiving more data for tha read
    *   adds the action to the action list
    *   @actionList : List of action messages
    *   @response   : action response object storing infos like read number and channel number
    */
    void Data::addStopReceivingDataAction(GetLiveReadsRequest_Actions* actionList, uint32_t channelNr, uint32_t readNr)
    {
        // create an action message and add it to the action list
        GetLiveReadsRequest_Action* action = actionList->add_actions();
        // FlowCell channel number on which we call the action
        action->set_channel(channelNr);

        // define action as stop_further_data action
        GetLiveReadsRequest_StopFurtherData* data = action->mutable_stop_further_data();
        *data = action->stop_further_data();

        // specify read number for which we don't want to receive more data
        // but only for this read, for the next one we want to get live signals
        action->set_number(readNr);

        // create unique action id by using a timestamp
        std::stringstream buf;
        buf << "stop_receiving_" << channelNr << "_" << readNr;
        action->set_action_id(buf.str());
    }

    /**
    *   take read action responses from the queue and add unblock and/or stop_receiving_data
    *   action to the action list and write action request to the bidirectional stream
    *   @action_queue   : safe queue with reads for which action messages shall be sent to MinKNOW
    */
    void Data::sendActions(SafeQueue<RTPair>& action_queue, SafeQueue<Durations>& duration_queue)
    {
        data_logger->info("Start sending unblock actions to MinKNOW");
        data_logger->flush();
        // as long as signals are received from MinKnow
        // iterate over received data and stop further data allocation for every odd read on every odd channel
        CSVFile /= "read_until_decision_stats.csv";
        csvfile csv(CSVFile.string());
        // header
        csv << " " << "read_id" << "channel_nr" << "read_nr" << "sequence_length" << "decision" << "duration" << endrow;

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
                    RTPair rp = std::move(action_queue.pop());
                    if (rp.first.unblock)
                    {
                        addUnblockAction(actionList, rp.first.channelNr, rp.first.readNr, 0.1);
                        counter++;
                        rp.second.timeCompleteRead.stop();
                        csv << std::to_string(counter) << rp.first.id << rp.first.channelNr
                            << rp.first.readNr << rp.first.sequence.size() << "unblock" 
                            << std::to_string(rp.second.timeCompleteRead.elapsed()) << endrow;
                    }
                    else
                    {
                        // add stop_receiving_data message for every read in the queue
                        addStopReceivingDataAction(actionList, rp.first.channelNr, rp.first.readNr);
                        counter++;
                        rp.second.timeCompleteRead.stop();
                        csv << std::to_string(counter) << rp.first.id << rp.first.channelNr
                            << rp.first.readNr << rp.first.sequence.size() << "stop_receiving"
                            << std::to_string(rp.second.timeCompleteRead.elapsed()) << endrow;
                    }
                    
                    if (rp.first.unblock)
                    {
                        duration_queue.push(Durations{
                            rp.second.timeCompleteRead.elapsed(),
                            -1,
                            rp.second.timeBasecallRead.elapsed(),
                            rp.second.timeClassifyRead.elapsed(),
                            });
                        if (rp.second.timeCompleteRead.elapsed() > 100.0)
                        {
                            std::cerr << rp.second.timeCompleteRead.begin() << "\t" << rp.second.timeCompleteRead.end() << std::endl;
                        }
                    }
                    else
                    {
                        duration_queue.push(Durations{
                            -1,
                            rp.second.timeCompleteRead.elapsed(),
                            rp.second.timeBasecallRead.elapsed(),
                            rp.second.timeClassifyRead.elapsed(),
                            });

                        if (rp.second.timeCompleteRead.elapsed() > 100.0)
                        {
                            std::cerr << rp.second.timeCompleteRead.begin() << "\t" << rp.second.timeCompleteRead.end() << std::endl;
                        }
                    }
                    
                }
                else
                    break;
            }

            // increase number of action messages in one request if queue is too full
            // decrease otherwise
            adaptActionBatchSize(action_queue.size());

            // write action request to the stream -> send message to MinKNOW
            // try 5 times to send message
            // throw exception if that was not successfull
            for (uint8_t i = 1; i <= 5; ++i)
            {
                if (stream->Write(actionRequest))
                    break;
                if (i == 5)
                {
                    data_logger->error("Failed sending action request to MinKNOW");
                    data_logger->flush();
                    throw FailedActionRequestException("Could not send unblock actions to MinKNOW!");
                }
                data_logger->warn("Failed sending action request number " + i);
                data_logger->flush();
                // wait for 0.4 seconds before trying to send the request again
                std::this_thread::sleep_for(std::chrono::milliseconds(400));
            }

            // there should be at least 100 ms between two separate action requests
            // otherwise MinKNOW action request queue is overflooded and drops requests
            end = std::chrono::system_clock::now();
            int elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            if(elapsed_milliseconds < 400)
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(400 - elapsed_milliseconds));
            }
            
        }
        data_logger->info("Finished sending unblock actions to MinKNOW");
        data_logger->flush();
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
    *   @throws : FailedSetupMessageException
    */
    void Data::startLiveStream()
    {
        data_logger->info("Trying to send setup message to MinKNOW");
        data_logger->flush();

        runs = true;

        // start streaming live nanopore signals
        stream = stub->get_live_reads(&context);

        // create setup messae
        GetLiveReadsRequest setupRequest{};
        GetLiveReadsRequest_StreamSetup* setup = setupRequest.mutable_setup();

        // we want to receive signals from all 512 channels of the MinION
        // has to be changed in case Flongle or PromethION is used
        // TODO: set last channel based on device type

        setup->set_first_channel((int)minChannel);
        setup->set_last_channel((int)maxChannel);

        // we only want to receive calibrated data
        setup->set_raw_data_type(GetLiveReadsRequest_RawDataType_CALIBRATED);
        // minimum number of signals for part of a read to be sent via the stream
        // zero means no limitation
        // TODO: Try to find out if other parameter like e.g. 4 is better
        setup->set_sample_minimum_chunk_size(0);
        setup->set_max_unblock_read_length_samples(0);

        // write setup message to stream and throw exception if that was not successful
        if (!stream->Write(setupRequest))
        {
            data_logger->error("Failed starting live stream : Could not send setup message!");
            data_logger->flush();
            throw FailedSetupMessageException("Failed starting live stream : Could not send setup message!");
        }

        data_logger->info("Setup message successfully send to MinKNOW.");
        data_logger->flush();
    }

	/**
    *   pull live nanopore signals from the stream and add the reads to the basecalling queue
    *   @basecall_queue : safe queue for storing reads ready for basecalling
    *   
    *   TODO: analyze action responses and log success and failed actions
    */
    void Data::getLiveSignals(SafeQueue<RTPair>& basecall_queue)
    {
        data_logger->info("Live signal receiving thread started");
        data_logger->flush();
        GetLiveReadsResponse response;

        uint32_t success = 0;
        uint32_t finished = 0;
        uint32_t too_long = 0;
        uint32_t readTag = 0;

        StopClock::TimePoint begin = StopClock::Clock::now();

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

            for (GetLiveReadsResponse_ActionResponse actResp : response.action_responses())
            {
                if (actResp.response() == GetLiveReadsResponse_ActionResponse_Response_SUCCESS)
                    success++;
                else  if (actResp.response() == GetLiveReadsResponse_ActionResponse_Response_FAILED_READ_FINISHED)
                    finished++;
                else
                    too_long++;
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
                ONTRead read{};
                read.channelNr = entry.first;
                read.readNr = entry.second.number();
                read.id = entry.second.id();
                read.raw_signals = string_to_float(entry.second.raw_data());
                read.readTag = ++readTag;
                // add read to basecalling queue
                // nanopore signal string is converted into a vector of floating point numbers
                basecall_queue.push(std::make_pair(std::move(read), times));
       		}

            StopClock::TimePoint end = StopClock::Clock::now();
            std::chrono::duration< StopClock::Seconds > elapsed = end - begin;
            if (elapsed.count() > 60.0)
            {
                data_logger->info("----------------------------- Intermediate Results -------------------------------------------------------");
                std::stringstream sstr;
                sstr << "Number of successfully unblocked reads    : " << success;
                data_logger->info(sstr.str());
                sstr.str("");
                sstr << "Number of failed finished reads           : " << finished;
                data_logger->info(sstr.str());
                sstr.str("");
                sstr << "Number of failed too long reads           : " << too_long;
                data_logger->info(sstr.str());
                data_logger->info("----------------------------------------------------------------------------------------------------------");
                data_logger->flush();
                begin = end;
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

