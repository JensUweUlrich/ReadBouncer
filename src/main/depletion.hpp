/*
 * live-depletion.hpp
 *
 *  Created on: 31.03.2021
 *      Author: jens-uwe.ulrich
 */

 // global variables
#include<iostream>
#include <vector>
#include <string>

std::filesystem::path NanoLiveRoot;

readuntil::Data* data;
SafeMap<std::string, std::pair<uint8_t, interleave::Read> > pending{};
double avgDurationCompleteClassifiedRead = 0;
double avgDurationCompleteUnClassifiedRead = 0;
double avgDurationBasecallRead = 0;
double avgDurationClassifyRead = 0;

// too short reads that need more data
SafeQueue<interleave::Read> classifiedReads{};
SafeQueue<interleave::Read> unclassifiedReads{};

//-------------------------------------------------------------------



//--------------------------------------------------------------------
// functions referenced as asynchronous tasks
//----------------------------------------------------------------------------------------------------------

/**
*	take read from the basecalling queue, perform basecalling and push that read on the classification queue
*	@basecall_queue			: safe queue with reads ready for basecalling
*	@classification_queue	: safe queue with basecalled reads
*	@weights				: weights file path needed by DeepNano to perform basecalling
*	@acq					: Acquisition service checking if sequencing run is already finished
*/
void basecall_live_reads(SafeQueue<readuntil::SignalRead>& basecall_queue,
    SafeQueue<interleave::Read>& classification_queue,
    std::string& weights,
    std::filesystem::path& weights_file,
    readuntil::Acquisition* acq)
{
    std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
    std::string f = weights_file.string();

    // create DeepNano2 caller object
    Caller* caller = create_caller(weights.c_str(), f.c_str(), 5, 0.01);

    // stop clock for looging basecall queue size regularly => only for debugging
    StopClock::TimePoint begin = StopClock::Clock::now();

    while (true)
    {
        if (!basecall_queue.empty())
        {
            readuntil::SignalRead read = basecall_queue.pop();

            // if we see data from this read for the first time
            // basecall signals and push to classification queue
            if (!pending.contains(read.id))
            {
                read.processingTimes.timeBasecallRead.start();
                char* sequence = call_raw_signal(caller, read.raw_signals.data(), read.raw_signals.size());
                read.processingTimes.timeBasecallRead.stop();
                classification_queue.push(interleave::Read(read.id, sequence, read.channelNr, read.readNr, read.processingTimes));
            }
            else
            {
                // if read was not classified after first and second try
                // concatenate sequence data from actual signals with former basecalled sequence data
                // use stop clock from former tries again
                TimeMeasures m = pending[read.id].second.getProcessingTimes();
                m.timeBasecallRead.start();
                char* sequence = call_raw_signal(caller, read.raw_signals.data(), read.raw_signals.size());
                m.timeBasecallRead.stop();
                std::stringstream sstr;
                sstr << pending[read.id].second.getSeq() << sequence;
                // push prolonged sequence to classification queue
                // longer read may be classified better
                classification_queue.push(interleave::Read(read.id, sstr.str(), read.channelNr, read.readNr, m));

            }
            StopClock::TimePoint end = StopClock::Clock::now();
            std::chrono::duration< StopClock::Seconds > elapsed = end - begin;
            if (elapsed.count() > 60.0)
            {
                std::stringstream sstr;
                sstr << "Size of Basecall Queue       :	" << basecall_queue.size();
                nanolive_logger->debug(sstr.str());
                nanolive_logger->flush();
                begin = end;
            }
        }

        if (acq->isFinished())
            break;

    }
}

/**
*	take basecalled reads from classification queue, try to find read in host IBF
*	push read on action queue if classified as host read with unblock label
*	if 3 chunks were not classified as non-host,
* 	push read on action queue with stop_receiving_data label
*	@classifcation_queue	: safe queue with basecalled reads ready for classification
*	@action_queue			: safe queue with reads for which action messages shall be sent to MinKNOW
*	@filters				: vector of host IBFs
*	@acq					: Acquisition service checking if sequencing run is already finished
*
*/
void classify_live_reads(SafeQueue<interleave::Read>& classification_queue,
    SafeQueue<readuntil::ActionResponse>& action_queue,
    std::vector<interleave::TIbf>& DepletionFilters,
    std::vector<interleave::TIbf>& TargetFilters,
    const double significance,
    const double error_rate,
    readuntil::Acquisition* acq)
{
    std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
    interleave::ClassifyConfig conf{};
    conf.strata_filter = -1;
    conf.significance = significance;
    conf.error_rate = error_rate;
    uint16_t found = 0;
    uint16_t failed = 0;
    double avgReadLen = 0.0;
    uint64_t rCounter = 0;
    std::set<std::string> onceSeen{};
    bool withTarget = !(TargetFilters.empty());

    StopClock::TimePoint begin = StopClock::Clock::now();
    while (true)
    {
        if (!classification_queue.empty())
        {
            interleave::Read read = classification_queue.pop();

            TimeMeasures m = read.getProcessingTimes();
            string readID = seqan::toCString(read.getID());
            // Average Read length for classifcation

            try
            {

                // only classify reads with length >= 300
                // reads with length < 300 will wait for the next read chunks and combine them
                // longer read can be better classified
                if (read.getSeqLength() < 300)
                {
                    // add readid and read entry to pending map
                    pending.insert({ readID , std::make_pair(0, read) });
                }
                else
                {
                    m.timeClassifyRead.start();
                    bool classified = false;
                    // if additional target filter is given
                    // read is classified for depletion iff it was found in depletion filters but not in target filters
                    if (withTarget)
                        classified = read.classify(DepletionFilters, conf) && !read.classify(TargetFilters, conf);
                    else
                        classified = read.classify(DepletionFilters, conf);
                    m.timeClassifyRead.stop();

                    if (classified)
                    {
                        // store all read data and time measures for classified read
                        action_queue.push(readuntil::ActionResponse{ read.getChannelNr(), read.getReadNr(),
                                            seqan::toCString(read.getID()), m, true });
                        classifiedReads.push(read);
                        pending.erase(readID);
                        avgReadLen += ((double)read.getSeqLength() - avgReadLen) / (double) ++rCounter;
                    }
                    else
                    {
                        // check if we already marked read as unclassified
                        // if read unclassified for the second time => stop receiving further data from this read
                        std::set<std::string>::iterator it = onceSeen.find(readID);
                        if (it != onceSeen.end())
                        {
                            action_queue.push(readuntil::ActionResponse{ read.getChannelNr(), read.getReadNr(),
                                                    seqan::toCString(read.getID()), pending[readID].second.getProcessingTimes() , false });
                            unclassifiedReads.push(read);
                            pending.erase(readID);
                            onceSeen.erase(it);
                            avgReadLen += ((double)read.getSeqLength() - avgReadLen) / (double) ++rCounter;
                        }
                        else
                        {
                            // read unclassified for the first time => insert in Ste of already seen reads
                            // but erase from pendin if first chunks were too short for classification
                            onceSeen.insert(readID);
                            pending.erase(readID);
                        }
                    }
                }
            }
            catch (std::exception& e)
            {
                std::stringstream estr;
                estr << "Error classifying Read : " << read.getID() << "(Len=" << read.getSeqLength() << ")";
                nanolive_logger->error(estr.str());
                estr.str("");
                estr << "Error message          : " << e.what();
                nanolive_logger->error(estr.str());
                nanolive_logger->flush();
            }

            StopClock::TimePoint end = StopClock::Clock::now();
            std::chrono::duration< StopClock::Seconds > elapsed = end - begin;
            if (elapsed.count() > 60.0)
            {
                std::stringstream sstr;
                sstr << "Size of Classification Queue       :	" << classification_queue.size();
                nanolive_logger->debug(sstr.str());
                sstr.str("");
                sstr << "Average Read Length                :	" << avgReadLen;
                nanolive_logger->debug(sstr.str());
                nanolive_logger->flush();
                begin = end;
            }
        }

        if (acq->isFinished())
            break;

    }
}


/**
*	compute average duration for complete processing time, basecalling time and classification time per read
*	using online average computation
*	@duration_queue	:	safe queue elapsed time for reads that finished sequencing
*	@acq			:	Acquisition service checking if sequencing run is already finished
*/
void compute_average_durations(SafeQueue<Durations>& duration_queue, readuntil::Acquisition* acq)
{
    std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
    uint64_t classifiedReadCounter = 0;
    uint64_t unclassifiedReadCounter = 0;
    StopClock::TimePoint begin = StopClock::Clock::now();
    while (true)
    {
        if (!duration_queue.empty())
        {
            Durations dur = duration_queue.pop();
            if (dur.completeClassified > -1)
            {
                classifiedReadCounter++;
                avgDurationCompleteClassifiedRead += (dur.completeClassified - avgDurationCompleteClassifiedRead) / classifiedReadCounter;
            }
            else
            {
                unclassifiedReadCounter++;
                avgDurationCompleteUnClassifiedRead += (dur.completeUnclassified - avgDurationCompleteUnClassifiedRead) / unclassifiedReadCounter;
            }

            avgDurationBasecallRead += (dur.basecalling - avgDurationBasecallRead) / (classifiedReadCounter + unclassifiedReadCounter);
            avgDurationClassifyRead += (dur.classification - avgDurationClassifyRead) / (classifiedReadCounter + unclassifiedReadCounter);

            StopClock::TimePoint end = StopClock::Clock::now();
            std::chrono::duration< StopClock::Seconds > elapsed = end - begin;
            if (elapsed.count() > 60.0)
            {
                nanolive_logger->info("----------------------------- Intermediate Results -------------------------------------------------------");
                std::stringstream sstr;
                sstr << "Number of classified reads                        :	" << classifiedReadCounter;
                nanolive_logger->info(sstr.str());
                sstr.str("");
                sstr << "Number of unclassified reads                      :	" << unclassifiedReadCounter;
                nanolive_logger->info(sstr.str());
                sstr.str("");
                sstr << "Average Processing Time for classified Reads      :	" << avgDurationCompleteClassifiedRead;
                nanolive_logger->info(sstr.str());
                sstr.str("");
                sstr << "Average Processing Time for unclassified Reads    :	" << avgDurationCompleteUnClassifiedRead;
                nanolive_logger->info(sstr.str());
                sstr.str("");
                sstr << "Average Processing Time Read Basecalling          :	" << avgDurationBasecallRead;
                nanolive_logger->info(sstr.str());
                sstr.str("");
                sstr << "Average Processing Time Read Classification       :	" << avgDurationClassifyRead;
                nanolive_logger->info(sstr.str());
                nanolive_logger->info("----------------------------------------------------------------------------------------------------------");
                nanolive_logger->flush();
                begin = end;
            }

        }

        if (acq->isFinished() && duration_queue.empty())
            break;
    }

    // print average duration times
    std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Number of classified reads						:	" << classifiedReadCounter << std::endl;
    std::cout << "Number of unclassified reads						:	" << unclassifiedReadCounter << std::endl;
    std::cout << "Average Processing Time for classified Reads		:	" << avgDurationCompleteClassifiedRead << std::endl;
    std::cout << "Average Processing Time for unclassified Reads	:	" << avgDurationCompleteUnClassifiedRead << std::endl;
    std::cout << "Average Processing Time Read Basecalling			:	" << avgDurationBasecallRead << std::endl;
    std::cout << "Average Processing Time Read Classification		:	" << avgDurationClassifyRead << std::endl;

}

void writeReads(SafeQueue<interleave::Read>& read_queue,
    readuntil::Acquisition* acq,
    const std::string output_file)
{
    seqan::SeqFileOut SeqOut;

    if (!seqan::open(SeqOut, seqan::toCString(output_file)))
    {
        std::cerr << "ERROR: Unable to open the file: " << output_file << std::endl;
        return;
    }

    while (true)
    {
        if (!read_queue.empty())
        {
            interleave::Read read = read_queue.pop();
            try
            {
                std::stringstream sstr;
                sstr << read.getID() << " channelNr=" << read.getChannelNr() << " readNr=" << read.getReadNr();
                seqan::writeRecord(SeqOut, sstr.str(), read.getSeq());

            }
            catch (seqan::Exception const& e)
            {
                std::cerr << "ERROR: " << e.what() << " [@" << read.getID() << "]" << std::endl;
                break;
            }
        }

        if (acq->isFinished() && read_queue.empty())
            break;
    }

    seqan::close(SeqOut);
}

/**
 *	core method for live read depletion
 *	@parser: input from the command line
 *  @throws: IBFBuildException
 */
/*void live_read_depletion(live_depletion_parser& parser)
{
    std::cout<<"Hi"<<std::endl;
   std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
    bool withTarget = false;
    // first check if basecalling file exists
    std::filesystem::path weights_file = NanoLiveRoot;
    weights_file.append("data");
    weights_file /= "rnn" + parser.weights + ".txt";
    if (!std::filesystem::exists(weights_file))
    {
        nanolive_logger->error("Could not find DeepNano weights file : " + weights_file.string());
        nanolive_logger->flush();
        //throw;
    }

    std::vector<interleave::TIbf> DepletionFilters{};
    std::vector<interleave::TIbf> TargetFilters{};
    // first load IBFs of host reference sequence
    if (parser.verbose)
        std::cout << "Loading Interleaved Bloom Filter(s) for depletion!" << ::std::endl;


    try
    {
        interleave::IBFConfig config{};
        interleave::IBF filter{};
        config.input_filter_file = parser.ibf_deplete_file;
        interleave::FilterStats stats = filter.load_filter(config);
        if (parser.verbose)
            interleave::print_stats(stats);
        DepletionFilters.emplace_back(filter.getFilter());
    }
    catch (interleave::IBFBuildException& e)
    {
        nanolive_logger->error("Could not load IBF File : " + std::string(e.what()));
        nanolive_logger->flush();
        //throw;
    }

    // parse target IBF if given as parameter
    if (parser.ibf_target_file.length() > 0)
    {
        try
        {
            interleave::IBFConfig config{};
            interleave::IBF filter{};
            config.input_filter_file = parser.ibf_target_file;
            interleave::FilterStats stats = filter.load_filter(config);
            TargetFilters.emplace_back(filter.getFilter());
            interleave::print_stats(stats);
            withTarget = true;
        }
        catch (interleave::ParseIBFFileException& e)
        {
            nanolive_logger->error("Error parsing target IBF using the following parameters");
            nanolive_logger->error("Target IBF file                : " + parser.ibf_target_file);
            nanolive_logger->error("Error message : " + std::string(e.what()));
            nanolive_logger->flush();
            //throw;
        }
    }


    if (parser.verbose)
    {
        std::cout << "Successfully loaded Interleaved Bloom Filter(s)!" << ::std::endl;
        std::cout << "Trying to connect to MinKNOW" << std::endl;
        std::cout << "Host : " << parser.host << std::endl;
        std::cout << "Port : " << parser.port << std::endl;
    }

    nanolive_logger->info("Successfully loaded Interleaved Bloom Filter(s)!");
    nanolive_logger->info("Trying to connect to MinKNOW");
    nanolive_logger->info("Host : " + parser.host);
    std::stringstream sstr;
    sstr << "Port : " << parser.port;
    nanolive_logger->info(sstr.str());
    nanolive_logger->flush();



    // create ReadUntilClient object and connect to specified device
    readuntil::ReadUntilClient& client = readuntil::ReadUntilClient::getClient();
    client.setHost(parser.host);
    client.setPort(parser.port);

    // TODO: throw exception if connection could not be established
    if (client.connect(parser.device))
    {
        if (parser.verbose)
            std::cout << "Connection successfully established!" << ::std::endl;
        else
        {
            nanolive_logger->info("Connection successfully established!");
            nanolive_logger->flush();
        }
    }
    else
    {
        std::cerr << "Could not establish connection to MinKNOW or MinION device" << std::endl;
        nanolive_logger->error("Could not establish connection to MinKNOW or MinION device (" + parser.device + ")");
        nanolive_logger->flush();
    }

    // wait until sequencing run has been started
    if (parser.verbose)
        std::cout << "Waiting for device to start sequencing!" << ::std::endl;

    std::cout << "Please start the sequencing run now!" << ::std::endl;

    readuntil::Acquisition* acq = (readuntil::Acquisition*)client.getMinKnowService(readuntil::MinKnowServiceType::ACQUISITION);
    if (acq->hasStarted())
    {
        if (parser.verbose)
            std::cout << "Sequencing has begun. Starting live signal processing!" << ::std::endl;

        nanolive_logger->info("Sequencing has begun. Starting live signal processing!");
        nanolive_logger->flush();

    }

    // create Data Service object
    // used for streaming live nanopore signals from MinKNOW and sending action messages back
    data = (readuntil::Data*)client.getMinKnowService(readuntil::MinKnowServiceType::DATA);

    // start live streaming of data
    try
    {
        data->startLiveStream();
    }
    catch (readuntil::DataServiceException& e)
    {
        nanolive_logger->error("Could not start streaming signals from device (" + parser.device + ")");
        nanolive_logger->error("Error message : " + std::string(e.what()));
        nanolive_logger->flush();
        throw;
    }

    // thread safe queue storing reads ready for basecalling
    SafeQueue<readuntil::SignalRead> basecall_queue{};
    // thread safe queue storing basecalled reads ready for classification
    SafeQueue<interleave::Read> classification_queue{};
    // thread safe queue storing classified reads ready for action creation
    SafeQueue<readuntil::ActionResponse> action_queue{};
    // thread safe queue storing for every read the duration for the different tasks to complete
    SafeQueue<Durations> duration_queue{};

    // start live signal streaming from ONT MinKNOW
    std::vector< std::future< void > > tasks;

    if (parser.verbose)
    {
        std::cout << "Start receiving live signals thread" << std::endl;
        std::cout << "Start basecalling thread" << std::endl;
        std::cout << "Start read classification thread" << std::endl;
        std::cout << "Start sending unblock messages thread" << std::endl;
    }


    // create thread for receiving signals from MinKNOW
    tasks.emplace_back(std::async(std::launch::async, &readuntil::Data::getLiveSignals, data, std::ref(basecall_queue)));

    // create thread for live basecalling
    tasks.emplace_back(std::async(std::launch::async, &basecall_live_reads, std::ref(basecall_queue),
        std::ref(classification_queue), std::ref(parser.weights), std::ref(weights_file), acq));

    // create thread/task for classification
    tasks.emplace_back(std::async(std::launch::async, &classify_live_reads, std::ref(classification_queue),
        std::ref(action_queue), std::ref(DepletionFilters), std::ref(TargetFilters),
        parser.kmer_significance, parser.error_rate, acq));

    // create thread/task for sending action messages back to MinKNOW
    tasks.emplace_back(std::async(std::launch::async, &readuntil::Data::sendActions, data, std::ref(action_queue), std::ref(duration_queue)));

    // create task for calculating average times needed to complete the different tasks
    tasks.emplace_back(std::async(std::launch::async, &compute_average_durations, std::ref(duration_queue), acq));

    // create task for writing classified and unclassified reads to output fasta files
    tasks.emplace_back(std::async(std::launch::async, &writeReads, std::ref(classifiedReads), acq, "classifiedReads.fasta"));
    tasks.emplace_back(std::async(std::launch::async, &writeReads, std::ref(unclassifiedReads), acq, "unclassifiedReads.fasta"));


    for (auto& task : tasks)
    {
        task.get();
    }

    data->stopLiveStream();


}*/

void live_read_depletion(live_depletion_parser& parser)
{
   std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
    bool withTarget = false;
    // first check if basecalling file exists
    std::filesystem::path weights_file = NanoLiveRoot;
    //weights_file.append("data");
    //weights_file.append("/weights/");
    weights_file = "rnn" + parser.weights + ".txt";
    if (std::filesystem::exists(weights_file))
    {
        std::cout<<"Found the weights file: "<<weights_file.string()<<std::endl;
        nanolive_logger->info("Found the weights file: " + weights_file.string());
        nanolive_logger->flush();
        //throw;
    }
    if (!std::filesystem::exists(weights_file))
    {
        std::cerr<<"Could not find DeepNano weights file: "<<weights_file.string()<<std::endl;
        nanolive_logger->error("Could not find DeepNano weights file : " + weights_file.string());
        nanolive_logger->flush();
        //throw;
    }

    std::vector<interleave::TIbf> DepletionFilters{};
    std::vector<interleave::TIbf> TargetFilters{};
    // first load IBFs of host reference sequence
    if (parser.verbose)
        std::cout << "Loading Interleaved Bloom Filter(s) for depletion!" << ::std::endl;


    try
    {
        interleave::IBFConfig config{};
        interleave::IBF filter{};
        config.input_filter_file = parser.ibf_deplete_file;
        interleave::FilterStats stats = filter.load_filter(config);
        if (parser.verbose)
            interleave::print_stats(stats);
        DepletionFilters.emplace_back(filter.getFilter());
    }
    catch (interleave::IBFBuildException& e)
    {
        nanolive_logger->error("Could not load IBF File : " + std::string(e.what()));
        nanolive_logger->flush();
        std::cerr<<"Could not load IBF File: "<<std::string(e.what())<<std::endl;
        //throw;
    }

    // parse target IBF if given as parameter
    if (parser.ibf_target_file.length() > 0)
    {
        try
        {
            interleave::IBFConfig config{};
            interleave::IBF filter{};
            config.input_filter_file = parser.ibf_target_file;
            interleave::FilterStats stats = filter.load_filter(config);
            TargetFilters.emplace_back(filter.getFilter());
            interleave::print_stats(stats);
            withTarget = true;
        }
        catch (interleave::ParseIBFFileException& e)
        {
            std::cerr<<"Error parsing target IBF using the following parameters"<<std::endl;
            std::cerr<<"Target IBF file: " + parser.ibf_target_file<<std::endl;
            std::cerr<<"Error message : " + std::string(e.what())<<std::endl;
            nanolive_logger->error("Error parsing target IBF using the following parameters");
            nanolive_logger->error("Target IBF file                : " + parser.ibf_target_file);
            nanolive_logger->error("Error message : " + std::string(e.what()));
            nanolive_logger->flush();
            //throw;
        }
    }


    if (parser.verbose)
    {
        std::cout << "Successfully loaded Interleaved Bloom Filter(s)!" << ::std::endl;
        std::cout << "Trying to connect to MinKNOW" << std::endl;
        std::cout << "Host : " << parser.host << std::endl;
        std::cout << "Port : " << parser.port << std::endl;
    }

    nanolive_logger->info("Successfully loaded Interleaved Bloom Filter(s)!");
    nanolive_logger->info("Trying to connect to MinKNOW");
    nanolive_logger->info("Host : " + parser.host);
    std::stringstream sstr;
    sstr << "Port : " << parser.port;
    nanolive_logger->info(sstr.str());
    nanolive_logger->flush();



    // create ReadUntilClient object and connect to specified device
    readuntil::ReadUntilClient& client = readuntil::ReadUntilClient::getClient();
    client.setHost(parser.host);
    client.setPort(parser.port);

    // TODO: throw exception if connection could not be established
    if (client.connect(parser.device))
    {
        if (parser.verbose)
            std::cout << "Connection successfully established!" << ::std::endl;
        else
        {
            std::cout << "Connection successfully established!" << ::std::endl;
            nanolive_logger->info("Connection successfully established!");
            nanolive_logger->flush();
        }
    }
    else
    {
        std::cerr << "Could not establish connection to MinKNOW or MinION device ("<<parser.device << ")"<< std::endl;
        nanolive_logger->error("Could not establish connection to MinKNOW or MinION device (" + parser.device + ")");
        nanolive_logger->flush();
    }

    // wait until sequencing run has been started
    if (parser.verbose)
        std::cout << "Waiting for device to start sequencing!" << ::std::endl;

    std::cout << "Please start the sequencing run now!" << ::std::endl;

    readuntil::Acquisition* acq = (readuntil::Acquisition*)client.getMinKnowService(readuntil::MinKnowServiceType::ACQUISITION);
    if (acq->hasStarted())
    {
        if (parser.verbose)
            std::cout << "Sequencing has begun. Starting live signal processing!" << ::std::endl;

        nanolive_logger->info("Sequencing has begun. Starting live signal processing!");
        nanolive_logger->flush();

    }

    // create Data Service object
    // used for streaming live nanopore signals from MinKNOW and sending action messages back
    data = (readuntil::Data*)client.getMinKnowService(readuntil::MinKnowServiceType::DATA);

    // start live streaming of data
    try
    {
        data->startLiveStream();
    }
    catch (readuntil::DataServiceException& e)
    {
        std::cerr<<"Could not start streaming signals from device ("<<parser.device<<")"<<std::endl;
        std::cerr<<"Error message: "<<std::string(e.what())<<std::endl;
        nanolive_logger->error("Could not start streaming signals from device (" + parser.device + ")");
        nanolive_logger->error("Error message : " + std::string(e.what()));
        nanolive_logger->flush();
        //throw;
    }

    // thread safe queue storing reads ready for basecalling
    SafeQueue<readuntil::SignalRead> basecall_queue{};
    // thread safe queue storing basecalled reads ready for classification
    SafeQueue<interleave::Read> classification_queue{};
    // thread safe queue storing classified reads ready for action creation
    SafeQueue<readuntil::ActionResponse> action_queue{};
    // thread safe queue storing for every read the duration for the different tasks to complete
    SafeQueue<Durations> duration_queue{};

    // start live signal streaming from ONT MinKNOW
    std::vector< std::future< void > > tasks;

    if (parser.verbose)
    {
        std::cout << "Start receiving live signals thread" << std::endl;
        std::cout << "Start basecalling thread" << std::endl;
        std::cout << "Start read classification thread" << std::endl;
        std::cout << "Start sending unblock messages thread" << std::endl;
    }


    // create thread for receiving signals from MinKNOW
    tasks.emplace_back(std::async(std::launch::async, &readuntil::Data::getLiveSignals, data, std::ref(basecall_queue)));

    // create thread for live basecalling
    tasks.emplace_back(std::async(std::launch::async, &basecall_live_reads, std::ref(basecall_queue),
        std::ref(classification_queue), std::ref(parser.weights), std::ref(weights_file), acq));

    // create thread/task for classification
    tasks.emplace_back(std::async(std::launch::async, &classify_live_reads, std::ref(classification_queue),
        std::ref(action_queue), std::ref(DepletionFilters), std::ref(TargetFilters),
        parser.kmer_significance, parser.error_rate, acq));

    // create thread/task for sending action messages back to MinKNOW
    tasks.emplace_back(std::async(std::launch::async, &readuntil::Data::sendActions, data, std::ref(action_queue), std::ref(duration_queue)));

    // create task for calculating average times needed to complete the different tasks
    tasks.emplace_back(std::async(std::launch::async, &compute_average_durations, std::ref(duration_queue), acq));

    // create task for writing classified and unclassified reads to output fasta files
    tasks.emplace_back(std::async(std::launch::async, &writeReads, std::ref(classifiedReads), acq, "classifiedReads.fasta"));
    tasks.emplace_back(std::async(std::launch::async, &writeReads, std::ref(unclassifiedReads), acq, "unclassifiedReads.fasta"));


    for (auto& task : tasks)
    {
        task.get();
    }

    data->stopLiveStream();


}


//-----------Qt

void live_read_depletion_qt(std::string ibf_deplete_file_name, std::string ibf_target_file_name, std::string host_name, int port, std::string device_name, std::string weights_name,double kmer_significance, double error_rate)
{
    std::cout<<"Starting...."<<std::endl;
    std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
    bool withTarget = false;
    // first check if basecalling file exists
    std::filesystem::path weights_file = NanoLiveRoot;
    weights_file.append("data");
    weights_file /= "rnn" + weights_name + ".txt";
    if (!std::filesystem::exists(weights_file))
    {
        nanolive_logger->error("Could not find DeepNano weights file : " + weights_file.string());
        //std::cerr<<"Could not find DeepNano weights file: "<<weights_file.string()<<std::endl;

       nanolive_logger->flush();
        //throw;
    }

    std::vector<interleave::TIbf> DepletionFilters{};
    std::vector<interleave::TIbf> TargetFilters{};
    // first load IBFs of host reference sequence
    //if (parser.verbose)
      std::cout << "Loading Interleaved Bloom Filter(s) for depletion!" << ::std::endl;


    try
    {
        interleave::IBFConfig config{};
        interleave::IBF filter{};
        config.input_filter_file = ibf_deplete_file_name;
        interleave::FilterStats stats = filter.load_filter(config);
        //if (parser.verbose)
            interleave::print_stats(stats);
        DepletionFilters.emplace_back(filter.getFilter());
    }
    catch (interleave::IBFBuildException& e)
    {
        //std::cerr<<"Could not load IBF File: "<<std::string(e.what())<<std::endl;

       nanolive_logger->error("Could not load IBF File : " + std::string(e.what()));
       nanolive_logger->flush();
        //throw;
    }

    // parse target IBF if given as parameter
    if (ibf_target_file_name.length() > 0)
    {
        try
        {
            interleave::IBFConfig config{};
            interleave::IBF filter{};
            config.input_filter_file = ibf_target_file_name;
            interleave::FilterStats stats = filter.load_filter(config);
            TargetFilters.emplace_back(filter.getFilter());
            interleave::print_stats(stats);
            withTarget = true;
        }
        catch (interleave::ParseIBFFileException& e)
        {
            //std::cerr<<"Error parsing target IBF using the following parameters"<<std::endl;
            nanolive_logger->error("Error parsing target IBF using the following parameters");
            //std::cerr<<"Target IBF file                : " + ibf_target_file_name<<std::endl;
            nanolive_logger->error("Target IBF file                : " + ibf_target_file_name);
            //std::cerr<<"Error message : " + std::string(e.what())<<std::endl;
            nanolive_logger->error("Error message : " + std::string(e.what()));
            nanolive_logger->flush();
           // throw;
        }
    }


    //if (parser.verbose)
    //{
        //std::cout << "Successfully loaded Interleaved Bloom Filter(s)!" << ::std::endl;
        //std::cout << "Trying to connect to MinKNOW" << std::endl;
        //std::cout << "Host : " << host_name << std::endl;
        //std::cout << "Port : " << port << std::endl;
    //}

    nanolive_logger->info("Successfully loaded Interleaved Bloom Filter(s)!");
    nanolive_logger->info("Trying to connect to MinKNOW");
    nanolive_logger->info("Host : " + host_name);
    std::stringstream sstr;
    sstr << "Port : " << port;
    std::cout<<sstr.str()<<std::endl;
    nanolive_logger->info(sstr.str());
    nanolive_logger->flush();



    // create ReadUntilClient object and connect to specified device
    readuntil::ReadUntilClient& client = readuntil::ReadUntilClient::getClient();
    client.setHost(host_name);
    client.setPort(port);

   // TODO: throw exception if connection could not be established
    try
    {
        if (client.connect(device_name))
        {
            std::cout << "Connection successfully established!" << std::endl;
            std::cout << "You can start live-depletion using these settings." << std::endl;
        }
    }
    catch (readuntil::DeviceServiceException& e)
    {
        std::cerr << "Connection to MinKNOW successfully established." << std::endl;
        std::cerr << "But could not detect given device/flowcell" << std::endl;
        std::cerr << "Please check whether the Flowcell has already been inserted. " << std::endl;
        //throw;
    }
    catch (readuntil::ReadUntilClientException& e)
    {
        std::cerr << "Could not establish connection to MinKNOW." << std::endl;
        std::cerr << "Please check the given host IP address and TCP port. " << std::endl;
        //throw;
    }

    std::cout << "Waiting for device to start sequencing!" << ::std::endl;

    std::cout << "Please start the sequencing run now!" << ::std::endl;

    // this line!
   //readuntil::Acquisition* acq = (readuntil::Acquisition*)client.getMinKnowService(readuntil::MinKnowServiceType::ACQUISITION);
    /*if (acq->hasStarted())
    {
        //if (parser.verbose)
           // std::cout << "Sequencing has begun. Starting live signal processing!" << ::std::endl;

       nanolive_logger->info("Sequencing has begun. Starting live signal processing!");
       nanolive_logger->flush();

    }*/

    // create Data Service object
    // used for streaming live nanopore signals from MinKNOW and sending action messages back
    data = (readuntil::Data*)client.getMinKnowService(readuntil::MinKnowServiceType::DATA);

    // start live streaming of data
   /* try
    {
        data->startLiveStream();
    }
    catch (readuntil::DataServiceException& e)
    {
        //std::cerr<<"Could not start streaming signals from device("<<device_name<<")"<<std::endl;
        //std::cerr<<"Error message: "<<std::string(e.what())<<std::endl;
        nanolive_logger->error("Could not start streaming signals from device (" + device_name + ")");
        nanolive_logger->error("Error message : " + std::string(e.what()));
        nanolive_logger->flush();
        //throw;
    }

    // thread safe queue storing reads ready for basecalling
   SafeQueue<readuntil::SignalRead> basecall_queue{};
    // thread safe queue storing basecalled reads ready for classification
    SafeQueue<interleave::Read> classification_queue{};
    // thread safe queue storing classified reads ready for action creation
    SafeQueue<readuntil::ActionResponse> action_queue{};
    // thread safe queue storing for every read the duration for the different tasks to complete
    SafeQueue<Durations> duration_queue{};

    // start live signal streaming from ONT MinKNOW
    std::vector< std::future< void > > tasks;

    //if (parser.verbose)
    //{
        std::cout << "Start receiving live signals thread" << std::endl;
        std::cout << "Start basecalling thread" << std::endl;
        std::cout << "Start read classification thread" << std::endl;
        std::cout << "Start sending unblock messages thread" << std::endl;
    //}


    // create thread for receiving signals from MinKNOW
    tasks.emplace_back(std::async(std::launch::async, &readuntil::Data::getLiveSignals, data, std::ref(basecall_queue)));

    // create thread for live basecalling
    tasks.emplace_back(std::async(std::launch::async, &basecall_live_reads, std::ref(basecall_queue),
        std::ref(classification_queue), std::ref(weights_name), std::ref(weights_file), acq));

    // create thread/task for classification
    tasks.emplace_back(std::async(std::launch::async, &classify_live_reads, std::ref(classification_queue),
        std::ref(action_queue), std::ref(DepletionFilters), std::ref(TargetFilters),
        kmer_significance, error_rate, acq));

    // create thread/task for sending action messages back to MinKNOW
    tasks.emplace_back(std::async(std::launch::async, &readuntil::Data::sendActions, data, std::ref(action_queue), std::ref(duration_queue)));

    // create task for calculating average times needed to complete the different tasks
    tasks.emplace_back(std::async(std::launch::async, &compute_average_durations, std::ref(duration_queue), acq));

    // create task for writing classified and unclassified reads to output fasta files
    tasks.emplace_back(std::async(std::launch::async, &writeReads, std::ref(classifiedReads), acq, "classifiedReads.fasta"));
    tasks.emplace_back(std::async(std::launch::async, &writeReads, std::ref(unclassifiedReads), acq, "unclassifiedReads.fasta"));


    for (auto& task : tasks)
    {
        task.get();
    }

    data->stopLiveStream();*/

}
