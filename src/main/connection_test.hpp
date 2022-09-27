/**
*	shift incoming signals directly as unblock response to action queue
*	needed only for unblock all
*	@signal_queue	:	thread safe queue with signal-only reads coming in from the sequencer
*	@action_queue	:	thread safe queue with unblock actions
*	@acq			:	Acquisition service checking if sequencing run is already finished
*/
void fill_action_queue(SafeQueue<RTPair>& signal_queue,
	SafeQueue<RTPair>& action_queue,
	readuntil::Acquisition* acq)
{
	while (true)
	{
		if (!signal_queue.empty())
		{
			RTPair rp = std::move(signal_queue.pop());
			rp.first.unblock = true;
			action_queue.push(std::move(rp));
		}

		if (acq->isFinished())
			break;
	}
}

/**
*	core function for testing connection to MinKNOW software and testing unblock all reads
*	@parser : input from the command line
*/
void test_connection(ConfigReader config)
{

	std::cout << "Trying to connect to MinKNOW" << std::endl;
	std::cout << "Host : " << config.MinKNOW_Parsed.host << std::endl;
	std::cout << "Port : " << config.MinKNOW_Parsed.port << std::endl;

	std::stringstream sstr;
	sstr << "Port : " << config.MinKNOW_Parsed.port;

	// create ReadUntilClient object and connect to specified device
	readuntil::ReadUntilClient& client = readuntil::ReadUntilClient::getClient();
	client.setHost(config.MinKNOW_Parsed.host);
	client.setPort(config.MinKNOW_Parsed.port); 
    client.setRootPath(ReadBouncerRoot);

	// TODO: throw exception if connection could not be established
	try
	{
		if (client.connect(config.MinKNOW_Parsed.flowcell))
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
		throw;
	}
	catch (readuntil::ReadUntilClientException& e)
	{
		std::cerr << "Could not establish connection to MinKNOW." << std::endl;
		std::cerr << "Please check the given host IP address and TCP port. " << std::endl;
		throw;
	}

	bool unblock_all = false;// as default and no changes in toml file! 

	if (unblock_all)
	{
		
		readuntil::AnalysisConfiguration* an_conf = (readuntil::AnalysisConfiguration*)client.getMinKnowService(readuntil::MinKnowServiceType::ANALYSIS_CONFIGURATION);
		an_conf->set_break_reads_after_seconds(0.4);
		// wait until sequencing run has been started
		//if (parser.verbose)
		std::cout << "Waiting for device to start sequencing!" << ::std::endl;

		std::cout << "Please start the sequencing run now!" << ::std::endl;

		readuntil::Acquisition* acq = (readuntil::Acquisition*)client.getMinKnowService(readuntil::MinKnowServiceType::ACQUISITION);
		if (acq->hasStarted())
		{
			//if (parser.verbose)
			std::cout << "Sequencing has begun. Starting live signal processing!" << ::std::endl;

            readbouncer_logger->info("Sequencing has begun. Starting live signal processing!");
            readbouncer_logger->flush();

		}

		// create Data Service object
		// used for streaming live nanopore signals from MinKNOW and sending action messages back
		data = (readuntil::Data*)client.getMinKnowService(readuntil::MinKnowServiceType::DATA);

		// set unblock all reads

		//(*data).setUnblockAll(true);
        readbouncer_logger->info("Unblocking all reads without basecalling or classification!");
        readbouncer_logger->flush();

		// start live streaming of data
		try
		{
			data->startLiveStream();
		}
		catch (readuntil::DataServiceException& e)
		{
            readbouncer_logger->error("Could not start streaming signals from device (" + config.MinKNOW_Parsed.flowcell + ")");
            readbouncer_logger->error("Error message : " + std::string(e.what()));
            readbouncer_logger->flush();
			throw;
		}

		

		// thread safe queue storing reads ready for basecalling
		SafeQueue<RTPair> read_queue{};
		// thread safe queue storing classified reads ready for action creation
		SafeQueue<RTPair> action_queue{};
		// thread safe queue storing for every read the duration for the different tasks to complete
		SafeQueue<Durations> duration_queue{};

		// start live signal streaming from ONT MinKNOW
		std::vector< std::future< void > > tasks;

		std::cout << "Start receiving live signals thread" << std::endl;
		std::cout << "Start sending unblock messages thread" << std::endl;


		// create thread for receiving signals from MinKNOW
		tasks.emplace_back(std::async(std::launch::async, &readuntil::Data::getLiveSignals, data, std::ref(read_queue)));

		// create thread for live basecalling
		tasks.emplace_back(std::async(std::launch::async, &fill_action_queue, std::ref(read_queue),
			std::ref(action_queue), acq));

		// create thread/task for sending action messages back to MinKNOW
		tasks.emplace_back(std::async(std::launch::async, &readuntil::Data::sendActions, data, std::ref(action_queue), std::ref(duration_queue)));

		for (auto& task : tasks)
		{
			task.get();
		}

		data->stopLiveStream();
	}
	
}

