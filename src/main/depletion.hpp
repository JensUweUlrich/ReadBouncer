/*
 * live-depletion.hpp
 *
 *  Created on: 31.03.2021
 *      Author: jens-uwe.ulrich
 */

struct Basecaller
{
	osprey::Decoder decoder;
	osprey::CallerDP caller;
};

 // global variables
std::filesystem::path NanoLiveRoot;

readuntil::Data* data;
//Caller* caller;
std::mutex callerMutex;
SafeMap<std::string, std::pair<uint8_t, interleave::Read> > pending{};
// too short reads that need more data
SafeQueue<interleave::Read> classifiedReads{};
SafeQueue<interleave::Read> unclassifiedReads{};

// write access only via mutex because several classification threads
// modify these variables
std::mutex avgReadLenMutex;
double avgReadLen = 0.0;
uint64_t rCounter = 0;

//-------------------------------------------------------------------



//--------------------------------------------------------------------
// functions referenced as asynchronous tasks
//----------------------------------------------------------------------------------------------------------

/**
*	take read from the basecalling queue, perform basecalling and push that read on the classification queue
*	@basecall_queue			: safe queue with reads ready for basecalling
*	@classification_queue	: safe queue with basecalled reads
*	@channel_stats			: safe map with number of received reads for every flow cell channel
*	@weights				: weights file path needed by DeepNano to perform basecalling
*	@acq					: Acquisition service checking if sequencing run is already finished
*/
void basecall_live_reads(SafeQueue<readuntil::SignalRead>& basecall_queue,
	SafeQueue<interleave::Read>& classification_queue,
	SafeMap<uint16_t, uint32_t>& channel_stats,
	Basecaller& bc,
	readuntil::Acquisition* acq)
{
	std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
	// create DeepNano2 caller object
	

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
				char* sequence = osprey::call_raw_signal(bc.caller, bc.decoder, read.raw_signals, read.processingTimes);
				//call_raw_signal(caller, read.raw_signals.data(), read.raw_signals.size());
				read.processingTimes.timeBasecallRead.stop();
				interleave::Read r = interleave::Read(read.id, sequence, read.channelNr, read.readNr, read.processingTimes);


				// only classify reads with length >= 300
				// reads with length < 300 will wait for the next read chunks and combine them
				if (r.getSeqLength() < 250)
				{
					//stop clock for overall processing time of the read
					// add readid and read entry to pending map
					//read.processingTimes.timeCompleteRead.stop();
					//r.setProcessingTimes(read.processingTimes);
					pending.insert({ read.id , std::make_pair(0, r) });
				}
				else
				{
					classification_queue.push(r);
				}
				uint32_t chNr = (uint32_t)read.channelNr;
				channel_stats.assign(std::pair(chNr, channel_stats[chNr] + 1));
				
			}
			else
			{
				// if read was not classified after first and second try
				// concatenate sequence data from actual signals with former basecalled sequence data
				// use stop clock from former tries again

				TimeMeasures m = pending[read.id].second.getProcessingTimes();
				read.processingTimes.timeBasecallRead.start();
				char* sequence = osprey::call_raw_signal(bc.caller, bc.decoder, read.raw_signals, read.processingTimes);
				read.processingTimes.timeBasecallRead.stop();

				// add elapsed basecall time of first chunks to second chunk
				read.processingTimes.timeBasecallRead.decrementStart(m.timeBasecallRead.runtime());
				// add elapsed overall time of first chunks to second chunk
				read.processingTimes.timeCompleteRead.setBegin(m.timeCompleteRead.begin());

				std::stringstream sstr;
				sstr << pending[read.id].second.getSeq() << sequence;
				// push prolonged sequence to classification queue
				// longer read may be classified better
				interleave::Read r = interleave::Read(read.id, sstr.str(), read.channelNr, read.readNr, read.processingTimes);
				if (r.getSeqLength() < 250)
				{
					// add readid and read entry to pending map
					//read.processingTimes.timeCompleteRead.stop();
					//r.setProcessingTimes(read.processingTimes);
					pending.insert({ read.id , std::make_pair(0, r) });
				}
				else
				{
					classification_queue.push(r);
					pending.erase(read.id);
				}
			}
			
		}

		if (acq->isFinished())
			break;

	}
}

/**
*	take basecalled reads from classification queue, try to find read in target IBF
*	push read on action queue if classified as target read with stop_further_data label
*   if read is not in target but in depletion filter, push read on action queue with unblock label
*	if 3 chunks were not classified target or depletion, push read on action queue with unblock label
*	@classifcation_queue	: safe queue with basecalled reads ready for classification
*	@action_queue			: safe queue with reads for which action messages shall be sent to MinKNOW
*	@filters				: vector of host IBFs
*	@acq					: Acquisition service checking if sequencing run is already finished
*
*/
void classify_target_reads(SafeQueue<interleave::Read>& classification_queue,
	SafeQueue<readuntil::ActionResponse>& action_queue,
	SafeMap<std::string, std::pair< interleave::Read, uint8_t>>& once_seen,
	std::vector<interleave::TIbf>& DepletionFilters,
	std::vector<interleave::TIbf>& TargetFilters,
	interleave::ClassifyConfig& conf,
	readuntil::Acquisition* acq)
{
	std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
	csvfile csv("read_until_decision_stats.csv");
	// header
	csv << " " << "read_id" << "channel_nr" << "read_nr" << "sequence_length" << "decision" << endrow;

	uint64_t read_counter = 0;
	bool withDepletion = !(DepletionFilters.empty());

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
				m.timeClassifyRead.start();
				bool stop_further = read.classify(TargetFilters, conf);
				bool unblock = false;
				
				if (withDepletion && !stop_further)
				{
					unblock = read.classify(DepletionFilters, conf);
				}
				m.timeClassifyRead.stop();

				// unblock all reads that not match the target filter
				// stop receiving further data for reads that belong to the target
				if (stop_further)
				{
					// store all read data and time measures for classified read
					std::unordered_map<std::string, std::pair<interleave::Read, uint8_t>>::iterator it = once_seen.find(readID);
					uint32_t readlen = read.getSeqLength();
					if (it != once_seen.end())
					{
						readlen += (*it).second.first.getSeqLength();
						once_seen.erase(it);
					}
					action_queue.push(readuntil::ActionResponse{ read.getChannelNr(), read.getReadNr(),
									seqan::toCString(read.getID()), readlen, m, false });
					classifiedReads.push(read);
					avgReadLenMutex.lock();
					avgReadLen += ((double)read.getSeqLength() - avgReadLen) / (double) ++rCounter;
					avgReadLenMutex.unlock();
				}
				// if read is not in target filter but in depletion filter => unblock read
				else if (unblock)
				{
					// store all read data and time measures for classified read
					std::unordered_map<std::string, std::pair<interleave::Read, uint8_t>>::iterator it = once_seen.find(readID);
					uint32_t readlen = read.getSeqLength();
					if (it != once_seen.end())
					{
						readlen += (*it).second.first.getSeqLength();
						once_seen.erase(it);
					}
					action_queue.push(readuntil::ActionResponse{ read.getChannelNr(), read.getReadNr(),
									seqan::toCString(read.getID()), readlen, m, true });
					unclassifiedReads.push(read);
					avgReadLenMutex.lock();
					avgReadLen += ((double)read.getSeqLength() - avgReadLen) / (double) ++rCounter;
					avgReadLenMutex.unlock();
				}
				else
				{
					// check if we already marked read as unclassified
					// if read is neither in target nor in depletion filter for second time => unblock this read 
					std::unordered_map<std::string, std::pair<interleave::Read, uint8_t>>::iterator it = once_seen.find(readID);
					if (it != once_seen.end())
					{
						//if ((*it).second.second == 1)
						//{
							uint32_t readlen = read.getSeqLength() + (*it).second.first.getSeqLength();
							once_seen.erase(it);
							action_queue.push(readuntil::ActionResponse{ read.getChannelNr(), read.getReadNr(),
											seqan::toCString(read.getID()), readlen, m , true });
							unclassifiedReads.push(read);
							
							avgReadLenMutex.lock();
							avgReadLen += ((double)read.getSeqLength() - avgReadLen) / (double) ++rCounter;
							avgReadLenMutex.unlock();
						/*}
						else
						{
							// if read unclassified for the second time => try another chunk
							once_seen.assign(std::pair(readID, std::pair(read, 2)));
						}*/
					}
					else
					{
						// read unclassified for the first time => insert in Ste of already seen reads
						// but erase from pendin if first chunks were too short for classification
						once_seen.assign(std::pair(readID, std::pair(read, 1)));
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
void classify_deplete_reads(SafeQueue<interleave::Read>& classification_queue,
	SafeQueue<readuntil::ActionResponse>& action_queue,
	SafeMap<std::string, std::pair< interleave::Read, uint8_t>>& once_seen,
	std::vector<interleave::TIbf>& DepletionFilters,
	std::vector<interleave::TIbf>& TargetFilters,
	interleave::ClassifyConfig& conf,
	readuntil::Acquisition* acq)
{
	std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
	
	
	uint64_t read_counter = 0;
	bool withTarget = !(TargetFilters.empty());

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
				m.timeClassifyRead.start();
				bool classified = false;

				
				// deplete reads if they match the depletion filter
				
				// if additional target filter is given
				// read is classified for depletion iff it was found in depletion filters but not in target filters
				if (withTarget)
					classified = read.classify(DepletionFilters, conf) && !read.classify(TargetFilters, conf);
				else
					classified = read.classify(DepletionFilters, conf);
				m.timeClassifyRead.stop();

				if (classified)
				{
					std::unordered_map<std::string, std::pair<interleave::Read, uint8_t>>::iterator it = once_seen.find(readID);
					uint32_t readlen = read.getSeqLength();
					if (it != once_seen.end())
					{
						readlen += (*it).second.first.getSeqLength();
						once_seen.erase(it);
					}
					// store all read data and time measures for classified read
					action_queue.push(readuntil::ActionResponse{ read.getChannelNr(), read.getReadNr(),
									seqan::toCString(read.getID()), readlen, m, true });
					classifiedReads.push(read);
					avgReadLenMutex.lock();
					avgReadLen += ((double)read.getSeqLength() - avgReadLen) / (double) ++rCounter;
					avgReadLenMutex.unlock();
					
					
				}
				else
				{
					// check if we already marked read as unclassified
					// if read unclassified for the third time => stop receiving further data from this read 
					std::unordered_map<std::string, std::pair<interleave::Read, uint8_t>>::iterator it = once_seen.find(readID);
					if (it != once_seen.end())
					{
						uint8_t iterstep = (*it).second.second;
						// after 10 kb of sequencing the read => stop receiving further data
						if (iterstep >= 5)
						{
							uint32_t readlen = read.getSeqLength() + (*it).second.first.getSeqLength();
							once_seen.erase(it);
							action_queue.push(readuntil::ActionResponse{ read.getChannelNr(), read.getReadNr(),
											seqan::toCString(read.getID()), readlen, m , false });
							unclassifiedReads.push(read);
							
							avgReadLenMutex.lock();
							avgReadLen += ((double)read.getSeqLength() - avgReadLen) / (double) ++rCounter;
							avgReadLenMutex.unlock();
						}
						else
						{
							// if read unclassified for the second time => try another chunk
							once_seen.assign(std::pair(readID, std::pair(read, ++iterstep)));
						}
					}
					else
					{
						// read unclassified for the first time => insert in Ste of already seen reads
						// but erase from pendin if first chunks were too short for classification
						once_seen.assign(std::pair(readID, std::pair(read, 1)));
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
							   SafeQueue<readuntil::SignalRead>& basecall_queue,
							   SafeQueue<interleave::Read>& classification_queue,
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
	double avgDurationRescaleRead = 0;
	double currentAvgDurationRescaleRead = 0;
	double avgDurationCallRead = 0;
	double currentAvgDurationCallRead = 0;
	double avgDurationBeamSearchRead = 0;
	double currentAvgDurationBeamSearchRead = 0;
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

			avgDurationRescaleRead += (dur.rescale - avgDurationRescaleRead) /
				(totalClassifiedReadCounter + totalUnclassifiedReadCounter + currentClassifiedReadCounter + currentUnclassifiedReadCounter);
			currentAvgDurationRescaleRead += (dur.rescale - currentAvgDurationRescaleRead) /
				(currentClassifiedReadCounter + currentUnclassifiedReadCounter);

			avgDurationCallRead += (dur.call - avgDurationCallRead) /
				(totalClassifiedReadCounter + totalUnclassifiedReadCounter + currentClassifiedReadCounter + currentUnclassifiedReadCounter);
			currentAvgDurationCallRead += (dur.call - currentAvgDurationCallRead) /
				(currentClassifiedReadCounter + currentUnclassifiedReadCounter);

			avgDurationBeamSearchRead += (dur.beam_search - avgDurationBeamSearchRead) /
				(totalClassifiedReadCounter + totalUnclassifiedReadCounter + currentClassifiedReadCounter + currentUnclassifiedReadCounter);
			currentAvgDurationBeamSearchRead += (dur.beam_search - currentAvgDurationBeamSearchRead) /
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
				sstr << "Average Processing Time Beam Search (interval)              :	" << currentAvgDurationBeamSearchRead;
				nanolive_logger->info(sstr.str());
				sstr.str("");
				sstr << "Average Processing Time Call (interval)                     :	" << currentAvgDurationCallRead;
				nanolive_logger->info(sstr.str());
				sstr.str("");
				sstr << "Average Processing Rescale (interval)                       :	" << currentAvgDurationRescaleRead;
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
void live_read_depletion(live_parser& parser, bool target_sequencing)
{
	std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
	bool withTarget = false;
	// first check if basecalling file exists
	//std::filesystem::path weights_file{ "C:\\ReadBouncer" };//NanoLiveRoot;
        std::filesystem::path weights_file{ "/home/jens/software/ReadBouncer/src/weights" };
//	weights_file.append("data");
	weights_file /= "net24dp.txt";
	if (!std::filesystem::exists(weights_file))
	{
		nanolive_logger->error("Could not find Basecaller weights file : " + weights_file.string());
		nanolive_logger->flush();
		throw;
	}

	if (parser.verbose)
		std::cout << "Initializing Basecaller!" << ::std::endl;

	std::vector<Basecaller> basecaller_vector{};

	std::string f = weights_file.string();
	std::string tables = f + ".tabs.txt";

	for (uint8_t t = 0; t < parser.basecall_threads; ++t)
	{
		osprey::Decoder decoder{ tables };
		//std::cout << "Decoder initialized" << std::endl;
		osprey::CallerDP caller{ weights_file.string() };
		//std::cout << "Caller initialized" << std::endl;
		basecaller_vector.push_back(Basecaller{decoder, caller});
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
	

	// create threads for live basecalling
	for (uint8_t t = 0; t < parser.basecall_threads; ++t)
	{
		tasks.emplace_back(std::async(std::launch::async, &basecall_live_reads, std::ref(basecall_queue),
			std::ref(classification_queue), std::ref(channelStats),
			std::ref(basecaller_vector[t]), acq));
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
