#include <string>
#include <vector>
#include <math.h>
#include <chrono>
#include <csignal>
#include <iostream>
#include <sstream>
#include <fstream>
#include <future>

#include "SafeQueue.hpp"
#include <StopClock.hpp>

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

#include <lyra/lyra.hpp>


// global variables
//-------------------------------------------------------------------

readuntil::Data *data;

std::shared_ptr<spdlog::logger> nanolive_logger;

double avgDurationCompleteClassifiedRead = 0;
double avgDurationCompleteUnClassifiedRead = 0;
double avgDurationBasecallRead = 0;
double avgDurationClassifyRead = 0;

//--------------------------------------------------------------------

//command line parser
//--------------------------------------------------------------------
/**
	class for generating the IBF build parser group
*/
struct ibf_build_parser
{
	std::string bloom_filter_output_path{ };
	std::string reference_file{};
    bool command = false;
    bool show_help = false;
	// never use uint*_t => not supported by lyra
	// will lead to unrecognized tokens
	int size_k = 13;
	int threads = 1;
	int fragment_size=10000;
	int filter_size = 0;
	bool verbose = false;

	/**
		parser constructor
		creates the ibfbuild group and adds it to the lyra cli 
		@cli: lyra command line interface object
	*/
    ibf_build_parser(lyra::cli& cli) 
    {
        cli.add_argument(
            lyra::command("ibfbuild",
                [this](const lyra::group & g) { this->do_command(g); }) 
                .help("Build Interleaved Bloom Filter with given references sequences")
                .add_argument(lyra::help(show_help))
                .add_argument(
                    lyra::opt(verbose)
                        .name("-v")
                        .name("--verbose")
                        .optional()
                        .help(
                            "Show additional output as to what we are doing."))
                .add_argument(
					lyra::opt(bloom_filter_output_path, "output-file")
						.name("-o")
						.name("--output-file")
						.required()
						.help("Output file of Interleaved Bloom Filter"))
				.add_argument(
					lyra::opt(reference_file, "input-reference")
						.name("-i")
						.name("--input-reference")
						.required()
						.help("Reference sequence file (fasta format) used to build the IBF; reads matching this reference will be filtered out"))
				.add_argument(
					lyra::opt(size_k, "kmer-size")
						.name("-k")
						.name("--kmer-size")
						.optional()
						.help("Kmer size used for building the Interleaved Bloom Filter"))
				.add_argument(
					lyra::opt(threads, "threads")
						.name("-t")
						.name("--threads")
						.optional()
						.help("Number of building threads"))
				.add_argument(
					lyra::opt(fragment_size, "fragment-size")
						.name("-f")
						.name("--fragment-size")
						.optional()
						.help("Length of fragments from the reference that are put in one bin of the IBF"))
				.add_argument(
					lyra::opt(filter_size, "filter-size")
						.name("-s")
						.name("--filter-size")
						.optional()
						.help("IBF size in MB"))
		);
				
    }

	/**
		function is called after parsing the group parameters from the command line
		prints the help page or the parameter values if option verbose is set
	*/
    void do_command(const lyra::group & g)
    {
        if (show_help)
            std::cout << g; 
        else
        {
			// trigger for calling the correct function after parsing the group parameters
			command = true;
			if (verbose)
			{
				std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
            	std::cout << "Build Interleaved Bloom Filter      : " << "verbose=" << (verbose ? "true" : "false") << std::endl;
            	std::cout << "Input reference file                : " << reference_file << std::endl;
				std::cout << "Output IBF file                     : " << bloom_filter_output_path << std::endl;
				std::cout << "Kmer size                           : " << size_k << std::endl;
				std::cout << "Size of reference fragments per bin : " << fragment_size << std::endl;
				std::cout << "IBF file size in MegaBytes          : " << filter_size << std::endl;
				std::cout << "Building threads                    : " << threads << std::endl;
				std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
			}
        }
    }
};

/**
	class for generating the IBF build parser group
*/
struct read_classify_parser
{
	std::string ibf_input_file{ };
	std::string read_file{};
    bool command = false;
    bool show_help = false;
	double kmer_significance = 0.95;
	double error_rate = 0.1;
	int threads = 1;
	bool verbose = false;

	/**
		parser constructor
		creates the classify group and adds it to the lyra cli 
		@cli: lyra command line interface object
	*/
    read_classify_parser(lyra::cli& cli) 
    {
        cli.add_argument(
            lyra::command("classify",
                [this](const lyra::group & g) { this->do_command(g); }) 
                .help("classify nanopore reads based on a given IBF file")
                .add_argument(lyra::help(show_help))
                .add_argument(
                    lyra::opt(verbose)
                        .name("-v")
                        .name("--verbose")
                        .optional()
                        .help(
                            "Show additional output as to what we are doing."))
                .add_argument(
					lyra::opt(read_file, "read-file")
						.name("-r")
						.name("--read-file")
						.required()
						.help("File with reads to classify (FASTA or FASTQ format)"))
				.add_argument(
					lyra::opt(ibf_input_file, "ibf-file")
						.name("-i")
						.name("--ibf-file")
						.required()
						.help("Interleaved Bloom Filter file"))
				.add_argument(
					lyra::opt(kmer_significance, "probability")
						.name("-s")
						.name("--significance")
						.optional()
						.help("significance level for confidence interval of number of errorneous kmers (default is 0.95"))
				.add_argument(
					lyra::opt(error_rate, "err")
					.name("-e")
					.name("--error-rate")
					.optional()
					.help("exepected per read sequencing error rate (default is 0.1"))
				.add_argument(
					lyra::opt(threads, "threads")
						.name("-t")
						.name("--threads")
						.optional()
						.help("Number of classification threads"))
		);
				
    }

	/**
		function is called after parsing the group parameters from the command line
		prints the help page or the parameter values if option verbose is set
	*/
    void do_command(const lyra::group & g)
    {
        if (show_help)
            std::cout << g; 
        else
        {
			// trigger for calling the correct function after parsing the group parameters
			command = true;
			if (verbose)
			{
				std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
            	std::cout << "Classify Reads                               : " << "verbose=" << (verbose ? "true" : "false") << std::endl;
            	std::cout << "Input read file                              : " << read_file << std::endl;
				std::cout << "Input IBF file                               : " << ibf_input_file << std::endl;
				std::cout << "Significance level for confidence interval   : " << kmer_significance << std::endl;
				std::cout << "Expected sequencing error rate               : " << error_rate << std::endl;
				std::cout << "Building threads                             : " << threads << std::endl;
				std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
			}
        }
    }
};

/**
	class for generating the IBF build parser group
*/
struct live_depletion_parser
{
	// default host & port to communicate with MinKNOW
	std::string host = "127.0.0.1";
	int port = 9501;
	std::string device{};
	std::string weights{};
	std::string ibf_input_file{ };
	double kmer_significance = 0.95;
	double error_rate = 0.1;
    bool command = false;
    bool show_help = false;
	bool verbose = false;
	bool unblock_all = false;

	/**
		parser constructor
		creates the live-deplete group and adds it to the lyra cli 
		@cli: lyra command line interface object
	*/
    live_depletion_parser(lyra::cli& cli) 
    {
        cli.add_argument(
            lyra::command("live-deplete",
                [this](const lyra::group & g) { this->do_command(g); }) 
                .help("Live classification and rejection of nanopore reads")
                .add_argument(lyra::help(show_help))
                .add_argument(
					lyra::opt(verbose)
                        .name("-v")
                        .name("--verbose")
                        .optional()
                        .help(
                            "Show additional output as to what we are doing."))
                .add_argument(
					lyra::opt(device, "device")
						.name("-d")
						.name("--device")
						.required()
						.help("Device or FlowCell name for live analysis"))
				.add_argument(
					lyra::opt(host, "host")
						.name("-c")
						.name("--host")
						.optional()
						.help("IP address on which MinKNOW software runs"))
				.add_argument(
					lyra::opt(port, "port")
						.name("-p")
						.name("--port")
						.optional()
						.help("MinKNOW communication port"))
				.add_argument(
					lyra::opt(ibf_input_file, "ibf-file")
						.name("-i")
						.name("--ibf-file")
						.required()
						.help("Interleaved Bloom Filter file"))
				.add_argument(
					lyra::opt(kmer_significance, "probability")
					.name("-s")
					.name("--significance")
					.optional()
					.help("significance level for confidence interval of number of errorneous kmers (default is 0.95"))
				.add_argument(
					lyra::opt(error_rate, "err")
					.name("-e")
					.name("--error-rate")
					.optional()
					.help("exepected per read sequencing error rate (default is 0.1"))
				.add_argument(
					lyra::opt(weights, "weights")
						.name("-w")
						.name("--weights")
						.optional()
						.help("Weights file"))
				.add_argument(
                    lyra::opt(unblock_all)
                        .name("-u")
                        .name("--unblock-all")
                        .optional()
                        .help(
                            "Unblock all reads"))
		);
				
    }

	/**
		function is called after parsing the group parameters from the command line
		prints the help page or the parameter values if option verbose is set
	*/
    void do_command(const lyra::group & g)
    {
        if (show_help)
            std::cout << g; 
        else
        {
			// trigger for calling the correct function after parsing the group parameters
			command = true;
			if (verbose)
			{
				std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
            	std::cout << "Live Nanopore Read Depletion                 : " << "verbose=" << (verbose ? "true" : "false") << std::endl;
            	std::cout << "Host IP address                              : " << host << std::endl;
				std::cout << "MinKNOW communication port                   : " << port << std::endl;
				std::cout << "Device or Flowcell name                      : " << device << std::endl;
				std::cout << "Input IBF file                               : " << ibf_input_file << std::endl;
				std::cout << "Significance level for confidence interval   : " << kmer_significance << std::endl;
				std::cout << "Expected sequencing error rate               : " << error_rate << std::endl;
				std::cout << "Unblock all live reads                       : " << (unblock_all ? "yes" : "no") << std::endl;
				std::cout << "Weights file for Live Basecalling            : " << weights << std::endl;
				std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
			}
        }
    }
};


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
					readuntil::Acquisition* acq)
{
	// create DeepNano2 caller object
	Caller* caller = create_caller("48", weights.c_str(), 5, 0.01);

	StopClock::TimePoint begin = StopClock::Clock::now();

	while (true)
	{
		if (!basecall_queue.empty())
		{
			readuntil::SignalRead read = basecall_queue.pop();
			read.processingTimes.timeBasecallRead.start();
			char* sequence = call_raw_signal(caller, read.raw_signals.data(), read.raw_signals.size());
			read.processingTimes.timeBasecallRead.stop();
			classification_queue.push(interleave::Read(read.id, sequence, read.channelNr, read.readNr, read.processingTimes));
			StopClock::TimePoint end = StopClock::Clock::now();
			std::chrono::duration< StopClock::Seconds > elapsed = end - begin;
			if (elapsed.count() > 60.0)
			{
				nanolive_logger->debug("Size of Basecall Queue       :	" + basecall_queue.size());
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
							std::vector<interleave::TIbf>& filters,
							const double significance,
							const double error_rate,
							readuntil::Acquisition* acq)
{
	
	interleave::DepleteConfig deplConf{};
	deplConf.strata_filter = -1;
	deplConf.significance = significance;
	deplConf.error_rate = error_rate;
	uint16_t found = 0;
	uint16_t failed = 0;
	double avgReadLen = 0.0;
	uint64_t rCounter = 0;

	// for unclassified read => store number of chunks that were unclassified
	//							and Time for the first data chunk
	std::map<seqan::CharString, std::pair<uint8_t, TimeMeasures> > unclassified{};
	StopClock::TimePoint begin = StopClock::Clock::now();
	while (true)
	{
		if (!classification_queue.empty())
		{
			interleave::Read read = classification_queue.pop();
			TimeMeasures m = read.getProcessingTimes();
			avgReadLen += ((double)read.getSeqLength() - avgReadLen) / (double) ++rCounter;
			m.timeClassifyRead.start();
			try
			{
				if (read.classify(filters, deplConf))
				{
					m.timeClassifyRead.stop();
					// store all read data and time measures for classified read
					action_queue.push(readuntil::ActionResponse{read.getChannelNr(), read.getReadNr(), 
										seqan::toCString(read.getID()), m, true});
				}
				else
				{
					// add readid entry to unclassified map
					// if read was unclassified for the third time -> add to action_queue with stop_further_data
					// store only processing times of the first chunk
					// helps to calculate time from getting first chunk to sending stopFurtherData message
					if (unclassified.find(read.getID()) == unclassified.end())
					{
						m.timeClassifyRead.stop();
						unclassified.insert({ read.getID() , std::pair(1, m) });
					}
					else
					{
						if (unclassified[read.getID()].first == 2)
						{
							// push read on action queue if unclassified for the third time
							// using time measures from the first chunk of data
							action_queue.push(readuntil::ActionResponse{ read.getChannelNr(), read.getReadNr(),
												seqan::toCString(read.getID()), unclassified[read.getID()].second , false });
						}
						else
							unclassified[read.getID()].first += 1;
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
				nanolive_logger->debug("Size of Classification Queue       :	" + classification_queue.size());
				std::stringstream sstr;
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
*	@acq			: Acquisition service checking if sequencing run is already finished
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
				nanolive_logger->info("Number of classified reads                        :	" + classifiedReadCounter);
				nanolive_logger->info("Number of unclassified reads                      :	" + unclassifiedReadCounter);
				std::stringstream sstr;
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

/**
 *	core method for live read depletion
 *	@parser: input from the command line
 *  @throws: IBFBuildException
 */ 
void live_read_depletion(live_depletion_parser& parser)
{
	// first load IBFs of host reference sequence

	if (parser.verbose)
		std::cout << "Loading Interleaved Bloom Filter(s)!" << ::std::endl;

	interleave::IBFConfig config{};
	interleave::IBF filter{};
	config.input_filter_file = parser.ibf_input_file;
	try
	{
		interleave::FilterStats stats = filter.load_filter(config);
		if (parser.verbose)
			interleave::print_stats(stats);
	}
	catch (interleave::IBFBuildException& e)
	{
		nanolive_logger->error("Could not load IBF File : " + std::string(e.what()));
		nanolive_logger->flush();
		throw;
	}
	
	std::vector<interleave::TIbf> filters{};
	filters.emplace_back(filter.getFilter());

	if (parser.verbose)
	{
		std::cout << "Successfully loaded Interleaved Bloom Filter(s)!" << ::std::endl;
		std::cout << "Trying to connect to MinKNOW" << std::endl;
		std::cout << "Host : " << parser.host << std::endl;
		std::cout << "Port : " << parser.port << std::endl;
	}
	else
	{
		nanolive_logger->info("Successfully loaded Interleaved Bloom Filter(s)!");
		nanolive_logger->info("Trying to connect to MinKNOW");
		nanolive_logger->info("Host : " + parser.host);
		nanolive_logger->info("Port : " + parser.port);
		nanolive_logger->flush();
	}


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

	// set unblock all reads
	if (parser.unblock_all)
		(*data).setUnblockAll(true);

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

	tasks.emplace_back(std::async(std::launch::async, &readuntil::Data::getLiveSignals, data, std::ref(basecall_queue)));

	// start t basecalling tasks/threads
	//for (uint8_t t = 0; t < 2; ++t)
	//{
		tasks.emplace_back(std::async(std::launch::async, &basecall_live_reads, std::ref(basecall_queue),
			std::ref(classification_queue), std::ref(parser.weights), acq));
	//}

	// create thread/task for classification
	tasks.emplace_back(std::async(std::launch::async, &classify_live_reads, std::ref(classification_queue), 
									std::ref(action_queue), std::ref(filters), parser.kmer_significance, parser.error_rate, acq));

	// create thread/task for sending action messages back to MinKNOW
	tasks.emplace_back(std::async(std::launch::async, &readuntil::Data::sendActions, data, std::ref(action_queue), std::ref(duration_queue)));

	// create task for calculating average times needed to complete the different tasks
	tasks.emplace_back(std::async(std::launch::async, &compute_average_durations, std::ref(duration_queue), acq));

	for (auto& task : tasks)
	{
		task.get();
	}

	data->stopLiveStream();

	/*
	try
	{
		
		
	}
	catch (readuntil::DataServiceException ex)
	{
		std::cerr << "Could not get live reads : " << ex.what() << std::endl;
		if (data != nullptr)
		{
			data->getContext()->TryCancel();
   		}
	}
	*/
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
                  interleave::TReads& 	reads )
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
			int64_t fragIdx = 0;
			int64_t seqlen = (int64_t) (length(seq));
			int64_t fragstart = fragIdx * 500 + 100;
			while (fragstart < (seqlen - 1))
			{
				std::string newid = std::string(seqan::toCString(id));
				std::size_t pos = newid.find(" ");
				newid = newid.substr(0,pos) + "_" + std::to_string(fragIdx); 
				uint64_t fragend = (fragIdx+1) * 500 + 100;
                // make sure that last fragment ends at last position of the reference sequence
                if (fragend > length(seq)) fragend = length(seq);
				seqan::Infix< seqan::CharString >::Type fragment = seqan::infix( seq, fragstart, fragend );
				reads.emplace_back(interleave::Read(newid, fragment));
                fragstart = ++fragIdx * 500 + 100;
				if (fragIdx > 0)
					break;
			}
            
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
* 
*/
void classify_reads(read_classify_parser& parser)
{
	interleave::IBFConfig config{};
	interleave::IBF filter {};
	config.input_filter_file = parser.ibf_input_file;
	interleave::FilterStats stats = filter.load_filter(config);
	interleave::print_stats(stats);
	interleave::TReads reads;
	parse_reads(parser.read_file, reads);
	std::vector<interleave::TIbf> filters{};
	filters.emplace_back(filter.getFilter());
	interleave::DepleteConfig deplConf{};
	deplConf.strata_filter = -1;

	deplConf.significance = parser.kmer_significance; 
	deplConf.error_rate = parser.error_rate;
	uint64_t found = 0;
	uint16_t failed = 0;
	uint64_t readCounter = 0;
	StopClock::Seconds avgClassifyduration = 0;
	// start stop clock
	StopClock::TimePoint begin = StopClock::Clock::now();
	for (interleave::Read r : reads)
	{
		readCounter++;
		StopClock classifyRead;
		classifyRead.start();
		try
		{
			if (r.classify(filters, deplConf))
				found++;
			std::cout << r.getID() << "\t" << r.getMaxKmerCount() << "\t" << r.getSeqLength() - config.kmer_size + 1 << std::endl;
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
	std::cout<<found << "/" << reads.size()<<std::endl;
	std::cout<<failed << "/" << reads.size()<<std::endl;
}

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

int main(int argc, char const **argv)
{
	std::signal(SIGINT, signalHandler);	

	initializeLogger();

	auto cli = lyra::cli();
	std::string command;
	bool show_help = false;
	cli.add_argument(lyra::help(show_help));
	ibf_build_parser ibfbuild_parser{cli};
	read_classify_parser classify_parser{cli};
	live_depletion_parser deplete_parser{cli};
	auto result = cli.parse({ argc, argv });
	if (!result)
    {
        std::cerr << result.errorMessage() << std::endl;
        std::cerr << cli;
        exit(1);
    }

    if(show_help)
    {
        std::cout << cli << '\n';
        exit(0);
    }

	try
	{
		if (ibfbuild_parser.command)
			buildIBF(ibfbuild_parser);
		else if (classify_parser.command)
			classify_reads(classify_parser);
		else if (deplete_parser.command)
			live_read_depletion(deplete_parser);
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}
	
	return 0;
}

