/*
 * live-depletion.hpp
 *
 *  Created on: 31.03.2021
 *      Author: jens-uwe.ulrich
 */

#include <ont_read.hpp>
#include <algorithm>

using namespace interfaces;

 // global variables
std::filesystem::path NanoLiveRoot;

readuntil::Data* data;
Runner runner{};
//Caller* caller;
std::mutex callerMutex;
SafeMap<std::string, std::pair<uint8_t, RTPair> > pending{};
// too short reads that need more data
SafeQueue<interleave::Read> DepletedReads{};
SafeQueue<interleave::Read> TargetReads{};

// write access only via mutex because several classification threads
// modify these variables
std::mutex avgReadLenMutex;
double avgReadLen = 0.0;
uint64_t rCounter = 0;

//-------------------------------------------------------------------
/*
*/
bool check_unblock( interleave::Read& read,
					interleave::ClassifyConfig& conf,
					std::vector<interleave::IBFMeta>& DepletionFilters,
					std::vector<interleave::IBFMeta>& TargetFilters)
{
	bool withTarget = !(TargetFilters.empty());
	bool withDepletion = !(DepletionFilters.empty());
	bool unblock = false;

	// depletion and target filters given
	if (withDepletion && withTarget)
	{
		std::pair<uint64_t, uint64_t> p = read.classify(DepletionFilters, TargetFilters, conf);
		// read matches with one of the depletion filters
		if (p.first > 0)
		{
			// read matches with one of the target filters as well
			if (p.second > 0)
			{
				// assume less sequencing errors to increase specificity
				conf.error_rate -= 0.02;
				p = read.classify(DepletionFilters, TargetFilters, conf);
				conf.error_rate += 0.02;


				if (p.first > 0 && p.second == 0)
				{
					// now read only matches to depletion filters
					unblock = true;
				}
				else
				{
					// read still matches to both => continue sequencing
					unblock = false;
				}
			}
			// read only matches to depletion filters
			else
				unblock = true;
		}
		// read does not match depletion filters
		else
			unblock = false;
	}
	// in depletion mode only => only unblock reads matching depletion filters
	else if (withDepletion)
		unblock = read.classify(DepletionFilters, conf) > -1;
	// if only target filters given => unblock reads not mapping the target
	else
	{
		int best_filter_index = read.classify(TargetFilters, conf);
		unblock = best_filter_index < 0;
		/*if (best_filter_index != -1)
		{
			TargetFilters[best_filter_index].classified += 1;
			if (parser.output_dir.length() > 0)
				seqan::writeRecord(TargetFilters[best_filter_index].outfile, id, seq);
		}
		*/
	}
	return unblock;
}

/*void add_target_reads_to_action_queue(bool unblock,
										interleave::Read& read,
										SafeQueue<RTPair>& classification_queue,
										SafeQueue<RTPair>& action_queue,
										SafeMap<std::string, std::pair< interleave::Read, uint8_t>>& once_seen,
										std::vector<interleave::IBFMeta>& DepletionFilters,
										std::vector<interleave::IBFMeta>& TargetFilters,
										interleave::ClassifyConfig& conf,
										RTPair& rp)
{
	if (!unblock)
	{
		std::unordered_map<std::string, std::pair<interleave::Read, uint8_t>>::iterator it = once_seen.find(rp.first.id);
		// concat with sequence from earlier seen read chunks that have not been classified
		//uint32_t readlen = read.getReadLength();
		if (it != once_seen.end())
		{
			//readlen += (*it).second.first.getReadLength();
			std::stringstream sstr;
			sstr << (*it).second.first.sequence << read.sequence;
			read.sequence = std::move((seqan::Dna5String)sstr.str());
			once_seen.erase(it);
		}
		// store all read data and time measures for classified read
		rp.first.unblock = false;
		action_queue.push(std::move(rp));
		TargetReads.push(read);
		avgReadLenMutex.lock();
		avgReadLen += ((double)read.getReadLength() - avgReadLen) / (double) ++rCounter;
		avgReadLenMutex.unlock();


	}
	else
	{
		// check if we already marked read as unclassified
		// if read unclassified for the third time => stop receiving further data from this read 
		std::unordered_map<std::string, std::pair<interleave::Read, uint8_t>>::iterator it = once_seen.find(rp.first.id);
		if (it != once_seen.end())
		{
			// prolong sequence with sequence data from earlier seen chunks
			uint8_t iterstep = (*it).second.second;
			std::stringstream sstr;
			sstr << (*it).second.first.sequence << read.sequence;
			read.sequence = std::move((seqan::Dna5String)sstr.str());
			unblock = check_unblock(read, conf, DepletionFilters, TargetFilters);

			if (unblock)
			{
				rp.first.unblock = true;
				action_queue.push(std::move(rp));
				DepletedReads.push(read);
				avgReadLenMutex.lock();
				avgReadLen += ((double)read.getReadLength() - avgReadLen) / (double) ++rCounter;
				avgReadLenMutex.unlock();
				once_seen.erase(it);
			}
			else
			{
				// we don't try to classify further if read is longer than 1.5 kb
				// we assume read to be on target
				if (read.getReadLength() > 1500) //iterstep >= 5)
				{
					once_seen.erase(it);
					rp.first.unblock = false;
					action_queue.push(rp);
					TargetReads.push(read);

					avgReadLenMutex.lock();
					avgReadLen += ((double)read.getReadLength() - avgReadLen) / (double) ++rCounter;
					avgReadLenMutex.unlock();
				}
				else
				{
					// if read unclassified for the second time => try another chunk
					once_seen.assign(std::pair(rp.first.id, std::pair(read, ++iterstep)));
				}
			}
		}
		else
		{
			// read unclassified for the first time => insert in SafeMap of already seen reads
			once_seen.assign(std::pair(rp.first.id, std::pair(read, 1)));
		}
	}
}
*/

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
void classify_live_reads(SafeQueue<RTPair>& classification_queue,
	SafeQueue<RTPair>& action_queue,
	SafeMap<std::string, std::pair< interleave::Read, uint8_t>>& once_seen,
	std::vector<interleave::IBFMeta>& DepletionFilters,
	std::vector<interleave::IBFMeta>& TargetFilters,
	interleave::ClassifyConfig& conf,
	readuntil::Acquisition* acq)
{
	std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
	
	bool target_only = !(TargetFilters.empty()) && DepletionFilters.empty();
	uint64_t read_counter = 0;
	
	while (true)
	{
		if (!classification_queue.empty())
		{
			RTPair rp = std::move(classification_queue.pop());
			seqan::Dna5String seq = (seqan::Dna5String) rp.first.sequence;
			interleave::Read read = interleave::Read(rp.first.id, seq);
			try
			{
				rp.second.timeClassifyRead.start();
				bool unblock = check_unblock(read, conf, DepletionFilters, TargetFilters);
				rp.second.timeClassifyRead.stop();

				// push classified reads to action queue
				if (unblock)
				{
					std::unordered_map<std::string, std::pair<interleave::Read, uint8_t>>::iterator it = once_seen.find(rp.first.id);
					// concat with sequence from earlier seen read chunks that have not been classified
					//uint32_t readlen = read.getReadLength();
					if (it != once_seen.end())
					{
						//readlen += (*it).second.first.getReadLength();
						std::stringstream sstr;
						sstr << (*it).second.first.sequence << read.sequence;
						read.sequence = std::move((seqan::Dna5String)sstr.str());
						once_seen.erase(it);
					}
					// store all read data and time measures for classified read
					//ReadHolder rp_new{ true, rp.first };
					rp.first.unblock = true;
					action_queue.push(std::move(rp));
					DepletedReads.push(read);
					avgReadLenMutex.lock();
					avgReadLen += ((double)read.getReadLength() - avgReadLen) / (double) ++rCounter;
					avgReadLenMutex.unlock();
					
					
				}
				else
				{
					// check if we already marked read as unclassified
					// if read unclassified for the third time => stop receiving further data from this read 
					std::unordered_map<std::string, std::pair<interleave::Read, uint8_t>>::iterator it = once_seen.find(rp.first.id);
					if (it != once_seen.end())
					{
						// prolong sequence with sequence data from earlier seen chunks
						uint8_t iterstep = (*it).second.second;
						std::stringstream sstr;
						sstr << (*it).second.first.sequence << read.sequence;
						read.sequence = std::move((seqan::Dna5String)sstr.str());
						unblock = check_unblock(read, conf, DepletionFilters, TargetFilters);

						if (unblock)
						{
							rp.first.unblock = true;
							action_queue.push(std::move(rp));
							DepletedReads.push(read);
							avgReadLenMutex.lock();
							avgReadLen += ((double)read.getReadLength() - avgReadLen) / (double) ++rCounter;
							avgReadLenMutex.unlock();
							once_seen.erase(it);
						}
						else
						{
							// we don't try to classify further if read is longer than 1.5 kb
							// we assume read to be on target
							if (read.getReadLength() > 1500) //iterstep >= 5)
							{
								once_seen.erase(it);
								rp.first.unblock = false;
								action_queue.push(rp);
								TargetReads.push(read);

								avgReadLenMutex.lock();
								avgReadLen += ((double)read.getReadLength() - avgReadLen) / (double) ++rCounter;
								avgReadLenMutex.unlock();
							}
							else
							{
								// if read unclassified for the second time => try another chunk
								once_seen.assign(std::pair(rp.first.id, std::pair(read, ++iterstep)));
							}
						}
					}
					else
					{
						// read unclassified for the first time => insert in SafeMap of already seen reads
						once_seen.assign(std::pair(rp.first.id, std::pair(read, 1)));
					}
				}
			}
			catch (std::exception& e)
			{
				std::stringstream estr;
				estr << "Error classifying Read : " << rp.first.id << "(Len=" << read.getReadLength() << ")";
				nanolive_logger->error(estr.str());
				estr.str("");
				estr << "Error message          : " << e.what();
				nanolive_logger->error(estr.str());
				nanolive_logger->flush();
			}
		}

		if (acq->isFinished())
			break;

	}
}


/**
*	compute average duration for complete processing time, basecalling time and classification time per read
*	using online average computation
*	@duration_queue	        :	safe queue elapsed time for reads that finished sequencing
* *	@basecall_queue			: safe queue with reads ready for basecalling
*	@classification_queue	: safe queue with basecalled reads
*   @channel_stats	        : safe map with number of received reads for every flow cell channel
*	@acq			        :	Acquisition service checking if sequencing run is already finished
*/
void compute_average_durations(SafeQueue<Durations>& duration_queue, 
							   SafeQueue<RTPair>& basecall_queue,
							   SafeQueue<RTPair>& classification_queue,
	                           SafeMap<uint16_t, uint32_t>& channel_stats,
	                           readuntil::Acquisition* acq)
{
	std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
	uint64_t totalClassifiedReadCounter = 0;
	uint64_t totalUnclassifiedReadCounter = 0;
	uint64_t currentClassifiedReadCounter = 0;
	uint64_t currentUnclassifiedReadCounter = 0;
	double avgDurationCompleteClassifiedRead = 0;
	double currentAvgDurCompleteClassifiedReads = 0;
	double avgDurationCompleteUnClassifiedRead = 0;
	double currentAvgDurCompleteUnClassifiedRead = 0;
	double avgDurationBasecallRead = 0;
	double currentAvgDurationBasecallRead = 0;
	double avgDurationClassifyRead = 0;
	double currentAvgDurationClassifyRead = 0;
	StopClock::TimePoint begin = StopClock::Clock::now();
	while (true)
	{
		if (!duration_queue.empty())
		{
			Durations dur = duration_queue.pop();
			if (dur.completeClassified > -1)
			{
				currentClassifiedReadCounter++;
				avgDurationCompleteClassifiedRead += (dur.completeClassified - avgDurationCompleteClassifiedRead) / 
														(totalClassifiedReadCounter + currentClassifiedReadCounter);
				currentAvgDurCompleteClassifiedReads += (dur.completeClassified - currentAvgDurCompleteClassifiedReads) /
															currentClassifiedReadCounter;
			}
			else
			{
				currentUnclassifiedReadCounter++;
				avgDurationCompleteUnClassifiedRead += (dur.completeUnclassified - avgDurationCompleteUnClassifiedRead) / 
														(totalUnclassifiedReadCounter + currentUnclassifiedReadCounter);
				currentAvgDurCompleteUnClassifiedRead += (dur.completeClassified - currentAvgDurCompleteUnClassifiedRead) /
														currentUnclassifiedReadCounter;
			}

			avgDurationBasecallRead += (dur.basecalling - avgDurationBasecallRead) / 
				(totalClassifiedReadCounter + totalUnclassifiedReadCounter + currentClassifiedReadCounter + currentUnclassifiedReadCounter);
			currentAvgDurationBasecallRead += (dur.basecalling - currentAvgDurationBasecallRead) /
												(currentClassifiedReadCounter + currentUnclassifiedReadCounter);
			avgDurationClassifyRead += (dur.classification - avgDurationClassifyRead) / 
				(totalClassifiedReadCounter + totalUnclassifiedReadCounter + currentClassifiedReadCounter + currentUnclassifiedReadCounter);
			currentAvgDurationClassifyRead += (dur.classification - currentAvgDurationClassifyRead) / 
												(currentClassifiedReadCounter + currentUnclassifiedReadCounter);

			StopClock::TimePoint end = StopClock::Clock::now();
			std::chrono::duration< StopClock::Seconds > elapsed = end - begin;
			if (elapsed.count() > 60.0)
			{
				// calculate number of channels that were actively sequencing
				// during the last 60 seconds
				uint16_t activeChannels = 0;
				for (uint16_t ch = 1; ch < 513; ++ch)
				{
					if (channel_stats[ch] > 0)
					{
						activeChannels++;
						channel_stats.assign(std::pair(ch, 0));
					}
				}
				totalClassifiedReadCounter += currentClassifiedReadCounter;
				totalUnclassifiedReadCounter += currentUnclassifiedReadCounter;
				nanolive_logger->info("----------------------------- Intermediate Results -------------------------------------------------------");
				std::stringstream sstr;
				sstr << "Total Number of classified reads                            :	" << totalClassifiedReadCounter;
				nanolive_logger->info(sstr.str());
				sstr.str("");
				sstr << "Total Number of unclassified reads                          :	" << totalUnclassifiedReadCounter;
				nanolive_logger->info(sstr.str());
				sstr.str("");
				sstr << "Number of active sequencing channels                        :	" << activeChannels;
				nanolive_logger->info(sstr.str());
				sstr.str("");
				sstr << "Number of classified reads during last interval             :	" << currentClassifiedReadCounter;
				nanolive_logger->info(sstr.str());
				sstr.str("");
				sstr << "Number of unclassified reads during last interval           :	" << currentUnclassifiedReadCounter;
				nanolive_logger->info(sstr.str());
				sstr.str("");
				avgReadLenMutex.lock();
				sstr << "Total Average Read Length                                   :	" << avgReadLen;
				avgReadLenMutex.unlock();
				nanolive_logger->info(sstr.str());
				sstr.str("");
				sstr << "Average Processing Time for classified Reads (interval)     :	" << currentAvgDurCompleteClassifiedReads;
				nanolive_logger->info(sstr.str());
				sstr.str("");
				sstr << "Average Processing Time for unclassified Reads (interval)   :	" << currentAvgDurCompleteUnClassifiedRead;
				nanolive_logger->info(sstr.str());
				sstr.str("");
				sstr << "Average Processing Time Read Basecalling (interval)         :	" << currentAvgDurationBasecallRead;
				nanolive_logger->info(sstr.str());
				sstr.str("");
				sstr << "Average Processing Time Read Classification (interval)      :	" << currentAvgDurationClassifyRead;
				nanolive_logger->info(sstr.str());
				sstr.str("");
				sstr << "Size of Basecall Queue                            :	" << basecall_queue.size();
				nanolive_logger->debug(sstr.str());
				sstr.str("");
				sstr << "Size of Classification Queue                      :	" << classification_queue.size();
				nanolive_logger->debug(sstr.str());
				nanolive_logger->info("----------------------------------------------------------------------------------------------------------");
				nanolive_logger->flush();
				begin = end;
				currentClassifiedReadCounter = 0;
				currentUnclassifiedReadCounter = 0;
			}

		}

		if (acq->isFinished() && duration_queue.empty())
			break;
	}
	totalClassifiedReadCounter += currentClassifiedReadCounter;
	totalUnclassifiedReadCounter += currentUnclassifiedReadCounter;
	// print average duration times 
	std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "Number of classified reads						:	" << totalClassifiedReadCounter << std::endl;
	std::cout << "Number of unclassified reads						:	" << totalUnclassifiedReadCounter << std::endl;
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
				seqan::writeRecord(SeqOut, read.id, read.sequence);

			}
			catch (seqan::Exception const& e)
			{
				std::cerr << "ERROR: " << e.what() << " [@" << read.id << "]" << std::endl;
				break;
			}
		}

		if (acq->isFinished() && read_queue.empty())
			break;
	}

	seqan::close(SeqOut);
}

void checkRunning(Runner& runner, readuntil::Acquisition* acq)
{
	while (true)
	{
		if (acq->isFinished())
			runner.isRunning = false;

		std::this_thread::sleep_for(std::chrono::seconds(5));

		if (!runner.isRunning)
			break;
	}
}

/**
 *	core method for live read depletion
 *	@parser: input from the command line
 *  @throws: IBFBuildException
 */
void adaptive_sampling(live_parser& parser)
{
	std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
	bool withTarget = false;
#if !defined(ARM_BUILD)
	// first check if basecalling file exists
	std::filesystem::path weights_file = NanoLiveRoot;
	weights_file.append("data");
	weights_file /= "rnn48.txt";
	if (!std::filesystem::exists(weights_file))
	{
		nanolive_logger->error("Could not find DeepNano weights file : " + weights_file.string());
		nanolive_logger->flush();
		throw;
	}
#endif
	std::vector<interleave::IBFMeta> DepletionFilters{};
	std::vector<interleave::IBFMeta> TargetFilters{};
	// first load IBFs of host reference sequence
	if (parser.verbose)
		std::cout << "Loading Depletion Interleaved Bloom Filter(s)!" << ::std::endl;

	if (parser.ibf_deplete_file.length() > 0)
	{
		std::vector<std::string> vector_files = split(parser.ibf_deplete_file, ',');
		for (std::string n : vector_files)
		{
			interleave::IBFMeta filter{};
			filter.name = std::filesystem::path(n).stem().string();
			interleave::IBF tf{};
			interleave::IBFConfig DepleteIBFconfig{};
			try
			{
				DepleteIBFconfig.input_filter_file = n;
				interleave::FilterStats stats = std::move(tf.load_filter(DepleteIBFconfig));
				filter.filter = std::move(tf.getFilter());
				if (parser.verbose)
					interleave::print_load_stats(stats);
			}
			catch (interleave::ParseIBFFileException& e)
			{
				nanolive_logger->error("Error parsing depletion IBF using the following parameters");
				nanolive_logger->error("Depletion IBF file                : " + n);
				nanolive_logger->error("Error message : " + std::string(e.what()));
				nanolive_logger->flush();
				throw;
			}

			DepletionFilters.emplace_back(std::move(filter));
		}
	}

	if (parser.verbose)
		std::cout << "Loading Target Interleaved Bloom Filter(s)!" << ::std::endl;
	// parse target IBF if given as parameter
	if (parser.ibf_target_file.length() > 0)
	{
		std::vector<std::string> vector_files = split(parser.ibf_target_file, ',');
		for (std::string n : vector_files)
		{
			interleave::IBFMeta filter{};
			filter.name = std::filesystem::path(n).stem().string();
			interleave::IBF tf{};
			interleave::IBFConfig TargetIBFconfig{};
			try
			{
				TargetIBFconfig.input_filter_file = n;
				interleave::FilterStats stats = std::move(tf.load_filter(TargetIBFconfig));
				filter.filter = std::move(tf.getFilter());
				if (parser.verbose)
					interleave::print_load_stats(stats);
			}
			catch (interleave::ParseIBFFileException& e)
			{
				nanolive_logger->error("Error parsing target IBF using the following parameters");
				nanolive_logger->error("Target IBF file                : " + n);
				nanolive_logger->error("Error message : " + std::string(e.what()));
				nanolive_logger->flush();
				throw;
			}

			TargetFilters.emplace_back(std::move(filter));
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
	
	if (runner.isRunning = acq->hasStarted())
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
	SafeQueue<RTPair> basecall_queue{};
	// thread safe queue storing basecalled reads ready for classification
	SafeQueue<RTPair> classification_queue{};
	// thread safe queue storing classified reads ready for action creation
	SafeQueue<RTPair> action_queue{};
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

	tasks.emplace_back(std::async(std::launch::async, &checkRunning, std::ref(runner), acq));

	// create thread for receiving signals from MinKNOW
	tasks.emplace_back(std::async(std::launch::async, &readuntil::Data::getLiveSignals, data, std::ref(basecall_queue)));

	// create threads for live basecalling
	
	
	basecall::Basecaller* caller;
	if (strcmpi(parser.caller.c_str(), "guppy") == 0)
	{
		std::string basecall_host = parser.guppy_host + ":" + parser.guppy_port;
		std::string config_name = "dna_r9.4.1_450bps_fast";
		caller = new basecall::GuppyBasecaller(basecall_host, config_name);
	}
	else
		caller = new basecall::DeepNanoBasecaller(weights_file, parser.basecall_threads);
	
	tasks.emplace_back(std::async(std::launch::async, &basecall::Basecaller::basecall_live_reads, std::move(caller), std::ref(basecall_queue),
			std::ref(classification_queue), std::ref(channelStats), std::ref(runner)));
	

	// create classification config
	interleave::ClassifyConfig conf{};
	conf.strata_filter = -1;
	conf.significance = parser.kmer_significance;
	conf.error_rate = parser.error_rate;

	// create thread/task for classification
	for (uint8_t t = 0; t < parser.classify_threads; ++t)
	{
		
		tasks.emplace_back(std::async(std::launch::async, &classify_live_reads, std::ref(classification_queue),
			std::ref(action_queue), std::ref(once_seen), std::ref(DepletionFilters), std::ref(TargetFilters),
			std::ref(conf), acq));
	}

	// create thread/task for sending action messages back to MinKNOW
	tasks.emplace_back(std::async(std::launch::async, &readuntil::Data::sendActions, data, std::ref(action_queue), std::ref(duration_queue)));
	

	// create task for calculating average times needed to complete the different tasks
	tasks.emplace_back(std::async(std::launch::async, &compute_average_durations, std::ref(duration_queue), 
		               std::ref(basecall_queue), std::ref(classification_queue), std::ref(channelStats), acq));

	// create task for writing classified and unclassified reads to output fasta files
	tasks.emplace_back(std::async(std::launch::async, &writeReads, std::ref(DepletedReads), acq, "DepletedReads.fasta"));
	tasks.emplace_back(std::async(std::launch::async, &writeReads, std::ref(TargetReads), acq, "TargetReads.fasta"));


	for (auto& task : tasks)
	{
		task.get();
	}

	data->stopLiveStream();


}
