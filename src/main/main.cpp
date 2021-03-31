#include <string>
#include <vector>
#include <math.h>
#include <chrono>
#include <csignal>
#include <iostream>
#include <sstream>
#include <fstream>
#include <future>
#include <filesystem>

#include <SafeQueue.hpp>
#include <SafeMap.hpp>
#include <StopClock.hpp>
#include <NanoLiveExceptions.hpp>

// spdlog library
#include "spdlog/spdlog.h"
#include "spdlog/sinks/rotating_file_sink.h"

// ReadUntil library
#include "ReadUntilClient.hpp"
#include "Data.hpp"
// IBF library
#include "IBF.hpp"

// Basecalling library
#include "DeepNano2.h"

// command line parser
#include "parser.hpp"


// global variables
//-------------------------------------------------------------------

readuntil::Data *data;

std::shared_ptr<spdlog::logger> nanolive_logger;
std::filesystem::path NanoLiveRoot;

double avgDurationCompleteClassifiedRead = 0;
double avgDurationCompleteUnClassifiedRead = 0;
double avgDurationBasecallRead = 0;
double avgDurationClassifyRead = 0;

// reads not classified after first or second try
SafeMap<std::string, std::pair<uint8_t, interleave::Read> > pending{};
SafeQueue<interleave::Read> classifiedReads{};
SafeQueue<interleave::Read> unclassifiedReads{};

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
void classify_live_reads(	SafeQueue<interleave::Read>& classification_queue,
							SafeQueue<readuntil::ActionResponse>& action_queue,
							std::vector<interleave::TIbf>& DepletionFilters,
							std::vector<interleave::TIbf>& TargetFilters,
							const double significance,
							const double error_rate,
							readuntil::Acquisition* acq)
{
	
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
						action_queue.push(readuntil::ActionResponse{read.getChannelNr(), read.getReadNr(), 
											seqan::toCString(read.getID()), m, true});
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
void live_read_depletion(live_depletion_parser& parser)
{
	bool withTarget = false;
	// first check if basecalling file exists
	std::filesystem::path weights_file = NanoLiveRoot;
	weights_file.append("data");
	weights_file /= "rnn" + parser.weights + ".txt";
	if (!std::filesystem::exists(weights_file))
	{
		nanolive_logger->error("Could not find DeepNano weights file : " + weights_file.string());
		nanolive_logger->flush();
		throw ;
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
		throw;
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
    readuntil::ReadUntilClient &client = readuntil::ReadUntilClient::getClient();
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

	readuntil::Acquisition *acq = (readuntil::Acquisition*) client.getMinKnowService(readuntil::MinKnowServiceType::ACQUISITION);
	if (acq->hasStarted())
	{
		if (parser.verbose)
			std::cout << "Sequencing has begun. Starting live signal processing!" << ::std::endl;

		nanolive_logger->info("Sequencing has begun. Starting live signal processing!");
		nanolive_logger->flush();
		
	}

	// create Data Service object
	// used for streaming live nanopore signals from MinKNOW and sending action messages back
	data = (readuntil::Data*) client.getMinKnowService(readuntil::MinKnowServiceType::DATA);

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

	
}

void signalHandler(int signum)
{
	//TODO: shutdown the gRPC stream smoothly
    if (data != nullptr)
   	{
		data->getContext()->TryCancel();
    }
	exit(signum);
}


void parse_reads( std::string const& 	reads_file,
                  interleave::TReads& 	reads,
				  uint64_t prefixLength)
{
	seqan::SeqFileIn seqFileIn;
    if ( !seqan::open( seqFileIn, seqan::toCString( reads_file ) ) )
    {
        std::cerr << "ERROR: Unable to open the file: " << reads_file << std::endl;
        return;
    }

    seqan::CharString id;
    seqan::CharString seq;
    
    while ( !seqan::atEnd( seqFileIn ) )
    {
        try
        {
            seqan::readRecord( id, seq, seqFileIn );

			uint64_t fragend = prefixLength;
			// make sure that last fragment ends at last position of the reference sequence
			if (fragend > length(seq)) fragend = length(seq);
			seqan::Infix< seqan::CharString >::Type fragment = seqan::infix(seq, 0, fragend);
			reads.emplace_back(interleave::Read(id, fragment));
            
        }
        catch ( seqan::Exception const& e )
        {
            std::cerr << "ERROR: " << e.what() << " [@" << id << "]" << std::endl;
			break;
        }
    }
	seqan::close( seqFileIn );
}

/**
*	initialize config for and build IBF
*	@parser	: input from the command line for "build" command
*	@throws	: IBFBuildException
*/
void buildIBF(ibf_build_parser& parser)
{
	interleave::IBFConfig config{};

	config.reference_files.emplace_back(parser.reference_file);
	config.output_filter_file = parser.bloom_filter_output_path;
	config.kmer_size = parser.size_k;
	config.threads_build = parser.threads;
	config.fragment_length = parser.fragment_size;
	config.filter_size = parser.filter_size;
	config.verbose = parser.verbose;

	interleave::IBF filter {};
	try
	{
		interleave::FilterStats stats = filter.create_filter(config);
		interleave::print_stats(stats);
	}
	catch (const interleave::IBFBuildException& e)
	{
		nanolive_logger->error("Error building IBF using the following parameters");
		nanolive_logger->error("Input reference file                : " + parser.reference_file);
		nanolive_logger->error("Output IBF file                     : " + parser.bloom_filter_output_path);
		nanolive_logger->error("Kmer size                           : " + parser.size_k);
		nanolive_logger->error("Size of reference fragments per bin : " + parser.fragment_size);
		nanolive_logger->error("IBF file size in MegaBytes          : " + parser.filter_size);
		nanolive_logger->error("Building threads                    : " + parser.threads);
		nanolive_logger->error("Error message : " + std::string(e.what()));
		nanolive_logger->error("---------------------------------------------------------------------------------------------------");
		nanolive_logger->flush();
		throw;
	}
	
}

/**
*	classify reads from an input file based on given depletion and/or target filters
*	@parser	: command line input parameters
*/
void classify_reads(read_classify_parser& parser)
{
	// initialize depletion and target filters
	interleave::IBFConfig DepleteIBFconfig{};
	interleave::IBFConfig TargetIBFconfig{};
	interleave::IBF DepleteFilter {};
	interleave::IBF TargetFilter{};
	std::vector<interleave::TIbf> DepletionFilters{};
	std::vector<interleave::TIbf> TargetFilters{};

	bool deplete = false;
	bool target = false;

	// parse depletion IBF if given as parameter
	if (parser.ibf_deplete_file.length() > 0)
	{
		try
		{
			DepleteIBFconfig.input_filter_file = parser.ibf_deplete_file;
			interleave::FilterStats stats = DepleteFilter.load_filter(DepleteIBFconfig);
			DepletionFilters.emplace_back(DepleteFilter.getFilter());
			interleave::print_stats(stats);
			deplete = true;
		}
		catch (interleave::ParseIBFFileException& e)
		{
			nanolive_logger->error("Error parsing depletion IBF using the following parameters");
			nanolive_logger->error("Depletion IBF file                : " + parser.ibf_deplete_file);
			nanolive_logger->error("Error message : " + std::string(e.what()));
			nanolive_logger->flush();
			throw;
		}
	}

	// parse target IBF if given as parameter
	if (parser.ibf_target_file.length() > 0)
	{
		try
		{
			TargetIBFconfig.input_filter_file = parser.ibf_target_file;
			interleave::FilterStats stats = TargetFilter.load_filter(TargetIBFconfig);
			TargetFilters.emplace_back(TargetFilter.getFilter());
			interleave::print_stats(stats);
			target = true;
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

	// parse input reads
	interleave::TReads reads;
	parse_reads(parser.read_file, reads, parser.preLen);
	
	// create classification config
	interleave::ClassifyConfig Conf{};
	Conf.strata_filter = -1;
	Conf.significance = parser.kmer_significance; 
	Conf.error_rate = parser.error_rate;

	uint64_t found = 0;
	uint16_t failed = 0;
	uint64_t too_short = 0;
	uint64_t readCounter = 0;
	StopClock::Seconds avgClassifyduration = 0;
	// start stop clock
	StopClock::TimePoint begin = StopClock::Clock::now();


	// initialize classification output files
	seqan::SeqFileOut ClassifiedOut;
	seqan::SeqFileOut UnclassifiedOut;

	// only print classified reads to file if option given via command line
	if (parser.classified_file.length() > 0)
	{
		if (!seqan::open(ClassifiedOut, seqan::toCString(parser.classified_file)))
		{
			std::cerr << "ERROR: Unable to open the file: " << parser.classified_file << std::endl;
			return;
		}
	}
	// only print unclassified reads to file if option given via command line
	if (parser.unclassified_file.length() > 0)
	{
		if (!seqan::open(UnclassifiedOut, seqan::toCString(parser.unclassified_file)))
		{
			std::cerr << "ERROR: Unable to open the file: " << parser.unclassified_file << std::endl;
			return;
		}
	}
	

	for (interleave::Read r : reads)
	{
		readCounter++;
		StopClock classifyRead;
		classifyRead.start();
		try
		{
			// read length has to be at least the size of the prefix used for read classification
			if (r.getSeqLength() < parser.preLen)
			{
				too_short++;
				continue;
			}

			// only classify if read is in depletion filter but NOT in target filter
			if (deplete && target)
			{
				if (r.classify(DepletionFilters, Conf) && !r.classify(TargetFilters, Conf))
				{
					found++;
					if (parser.classified_file.length() > 0)
					{
						std::stringstream sstr;
						sstr << r.getID() << " channelNr=" << r.getChannelNr() << " readNr=" << r.getReadNr();
						seqan::writeRecord(ClassifiedOut, sstr.str(), r.getSeq());
					}
					//std::cout << ">" << r.getID() << " kmers=" << r.getMaxKmerCount() << " len=" << r.getSeqLength() - config.kmer_size + 1 << std::endl;
					//std::cout << r.getSeq() << std::endl;
				}
				else if (parser.unclassified_file.length() > 0)
				{
					std::stringstream sstr;
					sstr << r.getID() << " channelNr=" << r.getChannelNr() << " readNr=" << r.getReadNr();
					seqan::writeRecord(UnclassifiedOut, sstr.str(), r.getSeq());
				}
			}
			// only classify if read is in depletion filter
			else if (deplete)
			{
				if (r.classify(DepletionFilters, Conf))
				{
					found++;
					if (parser.classified_file.length() > 0)
					{
						std::stringstream sstr;
						sstr << r.getID() << " channelNr=" << r.getChannelNr() << " readNr=" << r.getReadNr();
						seqan::writeRecord(ClassifiedOut, sstr.str(), r.getSeq());
					}
					//std::cout << ">" << r.getID() << " kmers=" << r.getMaxKmerCount() << " len=" << r.getSeqLength() - config.kmer_size + 1 << std::endl;
					//std::cout << r.getSeq() << std::endl;
				}
				else if (parser.unclassified_file.length() > 0)
				{
					std::stringstream sstr;
					sstr << r.getID() << " channelNr=" << r.getChannelNr() << " readNr=" << r.getReadNr();
					seqan::writeRecord(UnclassifiedOut, sstr.str(), r.getSeq());
				}
			}
			// only classify if read is in target filter
			else
			{
				if (r.classify(TargetFilters, Conf))
				{
					found++;
					if (parser.classified_file.length() > 0)
					{
						std::stringstream sstr;
						sstr << r.getID() << " channelNr=" << r.getChannelNr() << " readNr=" << r.getReadNr();
						seqan::writeRecord(ClassifiedOut, sstr.str(), r.getSeq());
					}
					//std::cout << ">" << r.getID() << " kmers=" << r.getMaxKmerCount() << " len=" << r.getSeqLength() - TargetIBFconfig.kmer_size + 1 << std::endl;
					//std::cout << r.getSeq() << std::endl;
				}
				else if (parser.unclassified_file.length() > 0)
				{
					std::stringstream sstr;
					sstr << r.getID() << " channelNr=" << r.getChannelNr() << " readNr=" << r.getReadNr();
					seqan::writeRecord(UnclassifiedOut, sstr.str(), r.getSeq());
				}
			}
		}
		catch (std::exception& e)
		{
			failed++;
			std::stringstream estr;
			estr << "Error classifying Read : " << r.getID() << "(Len=" << r.getSeqLength() << ")";
			nanolive_logger->error(estr.str());
			estr.str("");
			estr << "Error message          : " << e.what();
			nanolive_logger->error(estr.str());
			nanolive_logger->flush();
		}
		classifyRead.stop();
		avgClassifyduration += (classifyRead.elapsed() - avgClassifyduration) / readCounter;
		std::chrono::duration< StopClock::Seconds > elapsed = classifyRead.end() - begin;
		if (elapsed.count() > 60.0)
		{
			std::stringstream sstr;
			nanolive_logger->info("------------------------------- Intermediate Results -------------------------------");
			sstr << "Number of classified reads                         :   " << found;
			nanolive_logger->info(sstr.str());
			sstr.str("");
			sstr << "Number of of too short reads (len < " << parser.preLen << ")   :   " << too_short;
			nanolive_logger->info(sstr.str());
			sstr.str("");
			sstr << "Number of all reads                                :   " << readCounter;
			nanolive_logger->info(sstr.str());
			sstr.str("");
			sstr << "Average Processing Time Read Classification        :   " << avgClassifyduration;
			nanolive_logger->info(sstr.str());
			nanolive_logger->info("-----------------------------------------------------------------------------------");
			nanolive_logger->flush();
			begin = classifyRead.end();
		}

	}
	if (parser.classified_file.length() > 0)
		seqan::close(ClassifiedOut);
	if (parser.unclassified_file.length() > 0)
		seqan::close(UnclassifiedOut);
	std::stringstream sstr;
	std::cout << "------------------------------- Final Results -------------------------------" << std::endl;
	std::cout << "Number of classified reads                         :   " << found << std::endl;
	std::cout << "Number of of too short reads (len < " << parser.preLen << ")          :   " << too_short << std::endl;
	std::cout << "Number of all reads                                :   " << readCounter << std::endl;
	std::cout << "Average Processing Time Read Classification        :   " << avgClassifyduration << std::endl;
	std::cout << "-----------------------------------------------------------------------------------" << std::endl;
}

/**
*	setup global Logger for NanoLIVE
*/
void initializeLogger()
{
	try
	{
		nanolive_logger = spdlog::rotating_logger_mt("NanoLiveLog", "logs/NanoLiveLog.txt", 1048576 * 5, 100);
		nanolive_logger->set_level(spdlog::level::debug);
	}
	catch (const spdlog::spdlog_ex& e)
	{
		std::cerr << "Log initialization failed: " << e.what() << std::endl;
	}
}

/**
*	shift incoming signals directly as unblock response to action queue
*	needed only for unblock all
*	@signal_queue	:	thread safe queue with signal-only reads coming in from the sequencer
*	@action_queue	:	thread safe queue with unblock actions
*	@acq			:	Acquisition service checking if sequencing run is already finished
*/
void fill_action_queue(SafeQueue<readuntil::SignalRead>& signal_queue, 
					   SafeQueue<readuntil::ActionResponse>& action_queue,
					   readuntil::Acquisition* acq)
{
	while (true)
	{
		if (!signal_queue.empty())
		{
			readuntil::SignalRead read = signal_queue.pop();
			action_queue.push(readuntil::ActionResponse{ read.channelNr, read.readNr,
										read.id, read.processingTimes, true });
		}

		if (acq->isFinished())
			break;
	}
}

/**
*	core function for testing connection to MinKNOW software and testing unblock all reads
*	@parser : input from the command line
*/
void test_connection(connection_test_parser& parser)
{
	std::cout << "Trying to connect to MinKNOW" << std::endl;
	std::cout << "Host : " << parser.host << std::endl;
	std::cout << "Port : " << parser.port << std::endl;

	std::stringstream sstr;
	sstr << "Port : " << parser.port;

	// create ReadUntilClient object and connect to specified device
	readuntil::ReadUntilClient& client = readuntil::ReadUntilClient::getClient();
	client.setHost(parser.host);
	client.setPort(parser.port);

	// TODO: throw exception if connection could not be established
	try
	{
		if (client.connect(parser.device))
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

	if (parser.unblock_all)
	{
		
		
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

		// set unblock all reads

		//(*data).setUnblockAll(true);
		nanolive_logger->info("Unblocking all reads without basecalling or classification!");
		nanolive_logger->flush();

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
		SafeQueue<readuntil::SignalRead> read_queue{};
		// thread safe queue storing classified reads ready for action creation
		SafeQueue<readuntil::ActionResponse> action_queue{};
		// thread safe queue storing for every read the duration for the different tasks to complete
		SafeQueue<Durations> duration_queue{};

		// start live signal streaming from ONT MinKNOW
		std::vector< std::future< void > > tasks;

		if (parser.verbose)
		{
			std::cout << "Start receiving live signals thread" << std::endl;
			std::cout << "Start sending unblock messages thread" << std::endl;
		}


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

int main(int argc, char const **argv)
{
	
	std::signal(SIGINT, signalHandler);	

	initializeLogger();
	std::string binPath = argv[0];
	NanoLiveRoot = binPath.substr(0, binPath.find("bin"));

	auto cli = lyra::cli();
	std::string command;
	bool show_help = false;
	cli.add_argument(lyra::help(show_help));
	ibf_build_parser ibfbuild_parser{cli};
	read_classify_parser classify_parser{cli};
	live_depletion_parser deplete_parser{cli};
	connection_test_parser connect_parser{ cli };
	auto result = cli.parse({ argc, argv });
	if (!result)
    {
        std::cerr << result.errorMessage() << std::endl;
        std::cerr << cli;
        exit(1);
    }
	
	

    if(show_help)
    {
        std::cout << cli << std::endl;
        exit(0);
    }


	try
	{
		if (ibfbuild_parser.command)
			buildIBF(ibfbuild_parser);
		else if (connect_parser.command)
			test_connection(connect_parser);
		else if (classify_parser.command)
			classify_reads(classify_parser);
		else if (deplete_parser.command)
			live_read_depletion(deplete_parser);
		else
			std::cout << cli << std::endl;
			
			
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}
	
	return 0;
}

