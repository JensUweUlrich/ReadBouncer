
/*
 * depletionQt.hpp
 *
 *  Created on: 31.03.2021
 *      Author: jens-uwe.ulrich
 */

/**
 *	core method for live read depletion
 *	@parser: input from the command line
 *  @throws: IBFBuildException
 */
void live_read_depletion_qt(live_parser_qt& parser, bool target_sequencing)
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
        nanolive_logger->error("Could not find DeepNano weights file : " + weights_file.string());
        nanolive_logger->flush();
        throw;
    }

    std::vector<interleave::TIbf> DepletionFilters{};
    std::vector<interleave::TIbf> TargetFilters{};
    // first load IBFs of host reference sequence
    if (parser.verbose)
        std::cout << "Loading Depletion Interleaved Bloom Filter(s)!" << ::std::endl;

    if (parser.ibf_deplete_file.length() > 0)
    {
        try
        {
            interleave::IBFConfig config{};
            interleave::IBF filter{};
            config.input_filter_file = parser.ibf_deplete_file;
            interleave::FilterStats stats = filter.load_filter(config);
            if (parser.verbose)
                interleave::print_load_stats(stats);
            DepletionFilters.emplace_back(filter.getFilter());
        }
        catch (interleave::IBFBuildException& e)
        {
            nanolive_logger->error("Could not load IBF File : " + std::string(e.what()));
            nanolive_logger->flush();
            throw;
        }
    }

    if (parser.verbose)
        std::cout << "Loading Target Interleaved Bloom Filter(s)!" << ::std::endl;
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
            if (parser.verbose)
                interleave::print_load_stats(stats);
            withTarget = true;
        }
        catch (interleave::ParseIBFFileException& e)
        {
            nanolive_logger->error("Error parsing target IBF using the following parameters");
            nanolive_logger->error("Target IBF file                : " + parser.ibf_target_file);
            nanolive_logger->error("Error message : " + std::string(e.what()));
            nanolive_logger->flush();
            throw;
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
    client.setRootPath(NanoLiveRoot);

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

    // set chunk size by changing break_reads_after_seconds
    // seems to be overturned by TOML file configuration
    readuntil::AnalysisConfiguration* ana_conf = (readuntil::AnalysisConfiguration*)client.getMinKnowService(readuntil::MinKnowServiceType::ANALYSIS_CONFIGURATION);
    ana_conf->set_break_reads_after_seconds(0.4);
    if (parser.verbose)
    {
        nanolive_logger->info("Set break_reads_after_seconds = 0.4");
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
    // thread safe set of reads which were too short after first basecalling
    SafeMap<std::string, std::pair<interleave::Read, uint8_t>> once_seen{};
    // thread safe map storing number of send reads for every channel
    SafeMap<uint16_t, uint32_t> channelStats{};
    for (uint16_t ch = 1; ch < 513; ++ch)
    {
        channelStats.insert(std::pair(ch, 0));
    }

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

    // create DeepNano2 caller object
    //std::string f = weights_file.string();
    // TODO: check if thread safe
    //caller = create_caller(parser.weights.c_str(), f.c_str(), 5, 0.01);


    // create threads for live basecalling
    for (uint8_t t = 0; t < parser.basecall_threads; ++t)
    {
        tasks.emplace_back(std::async(std::launch::async, &basecall_live_reads, std::ref(basecall_queue),
            std::ref(classification_queue), std::ref(channelStats), std::ref(parser.weights),
            std::ref(weights_file), acq));
    }

    // create classification config
    interleave::ClassifyConfig conf{};
    conf.strata_filter = -1;
    conf.significance = parser.kmer_significance;
    conf.error_rate = parser.error_rate;

    // create thread/task for classification
    for (uint8_t t = 0; t < parser.classify_threads; ++t)
    {
        if (target_sequencing)
            tasks.emplace_back(std::async(std::launch::async, &classify_target_reads, std::ref(classification_queue),
                std::ref(action_queue), std::ref(once_seen), std::ref(DepletionFilters), std::ref(TargetFilters),
                std::ref(conf), acq));
        else
            tasks.emplace_back(std::async(std::launch::async, &classify_deplete_reads, std::ref(classification_queue),
                std::ref(action_queue), std::ref(once_seen), std::ref(DepletionFilters), std::ref(TargetFilters),
                std::ref(conf), acq));
    }

    // create thread/task for sending action messages back to MinKNOW
    tasks.emplace_back(std::async(std::launch::async, &readuntil::Data::sendActions, data, std::ref(action_queue), std::ref(duration_queue)));


    // create task for calculating average times needed to complete the different tasks
    tasks.emplace_back(std::async(std::launch::async, &compute_average_durations, std::ref(duration_queue),
                       std::ref(basecall_queue), std::ref(classification_queue), std::ref(channelStats), acq));

    // create task for writing classified and unclassified reads to output fasta files
    tasks.emplace_back(std::async(std::launch::async, &writeReads, std::ref(classifiedReads), acq, "classifiedReads.fasta"));
    tasks.emplace_back(std::async(std::launch::async, &writeReads, std::ref(unclassifiedReads), acq, "unclassifiedReads.fasta"));


    for (auto& task : tasks)
    {
        task.get();
    }

    data->stopLiveStream();


}
