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
#include <SafeSet.hpp>
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

// tomel parser
//#include "parsertoml.hpp"
#include "configReader.hpp"
// toml library 
#include "../toml11/toml.hpp"

// Basecalling library
#if !defined(ARM_BUILD)
	#include <DeepNanoBasecaller.hpp>
#endif

#include <GuppyBasecaller.hpp>

// command line parser
#include "parser.hpp"

// subcommand related functions
#include "ibfbuild.hpp"
#include "classify.hpp"
#include "adaptive_sampling.hpp"

#if defined(_WIN32)
	#include <windows.h>
	#include <psapi.h>
	#pragma comment( lib, "psapi.lib" )
	
	static const unsigned __int64 epoch = ((unsigned __int64)116444736000000000ULL);
#else
	#include <sys/resource.h>
	#include <sys/time.h>
#endif

std::shared_ptr<spdlog::logger> nanolive_logger;
//std::filesystem::path NanoLiveRoot;
//readuntil::Data* data;

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
	client.setRootPath(NanoLiveRoot);

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

	readuntil::AnalysisConfiguration* an_conf = (readuntil::AnalysisConfiguration*)client.getMinKnowService(readuntil::MinKnowServiceType::ANALYSIS_CONFIGURATION);
	an_conf->set_break_reads_after_seconds(0.4);

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
		SafeQueue<RTPair> read_queue{};
		// thread safe queue storing classified reads ready for action creation
		SafeQueue<RTPair> action_queue{};
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

void signalHandler(int signum)
{
	//TODO: shutdown the gRPC stream smoothly
	if (data != nullptr)
	{
		data->getContext()->TryCancel();
	}
	runner.isRunning = false;
	exit(signum);
}


/**
*	setup global Logger for NanoLIVE
*/
void initializeLogger(std::string toml_file)
{
	std::ifstream tomlFileReadBouncer(toml_file, std::ios_base::binary);
	const auto configuration_ = toml::parse(tomlFileReadBouncer, /*optional -> */ toml_file);
	auto log_file = toml::find<std::string>(configuration_, "log_directory");


	try
	{
		nanolive_logger = spdlog::rotating_logger_mt("NanoLiveLog", log_file + "logs/NanoLiveLog.txt", 1048576 * 5, 100);
		nanolive_logger->set_level(spdlog::level::debug);
	}
	catch (const spdlog::spdlog_ex& e)
	{
		std::cerr << "Log initialization failed: " << e.what() << std::endl;
	}
}


#if defined(_WIN32)
/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t getPeakRSS()
{
	/* Windows -------------------------------------------------- */
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
	return (size_t)info.PeakWorkingSetSize;
}

// taken from https://stackoverflow.com/questions/5272470/c-get-cpu-usage-on-linux-and-windows
double cputime()
{
	HANDLE hProcess = GetCurrentProcess();
	FILETIME ftCreation, ftExit, ftKernel, ftUser;
	SYSTEMTIME stKernel;
	SYSTEMTIME stUser;

	GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser);
	FileTimeToSystemTime(&ftKernel, &stKernel);
	FileTimeToSystemTime(&ftUser, &stUser);

	double kernelModeTime = ((stKernel.wHour * 60.) + stKernel.wMinute * 60.) + stKernel.wSecond * 1. + stKernel.wMilliseconds / 1000.;
	double userModeTime = ((stUser.wHour * 60.) + stUser.wMinute * 60.) + stUser.wSecond * 1. + stUser.wMilliseconds / 1000.;

	return kernelModeTime + userModeTime;
}
#else
	double cputime(void)
	{
		struct rusage r;
		getrusage(RUSAGE_SELF, &r);
		return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
	}

	long getPeakRSS(void)
	{
		struct rusage r;
		getrusage(RUSAGE_SELF, &r);
		return r.ru_maxrss * 1024;
	}
#endif

/**
* Set the configuration structs for classify/target
* @param : config object, Toml output file, usage, a list of target and deplete files
*/

void inline configurationReader(configReader config, std::string const tomlFile, std::string subcommand, std::fstream& tomlOutput,
	std::string target_files_,
	std::string deplete_files_) 
{

	if (subcommand == "build") {

		configReader::ibf_build_parser_ struct_ = config.ibfReader(tomlOutput, subcommand, target_files_, deplete_files_);

		// Create a log of toml usage 
		auto tbl = toml::value{ {

		{subcommand, toml::table{{
				{ "target_files", target_files_ },
				{ "deplete_files", deplete_files_ },
				{ "kmer-size", struct_.size_k },
				{ "threads", struct_.threads },
				{ "fragment-size", struct_.fragment_size }
				 }}
				},
		} };

		// chrono: https://en.cppreference.com/w/cpp/chrono
		auto start = std::chrono::system_clock::now();
		auto end = std::chrono::system_clock::now();

		std::chrono::duration<double> elapsed_seconds = end - start;
		std::time_t end_time = std::chrono::system_clock::to_time_t(end);

		tomlOutput << "# Computation time: " << std::ctime(&end_time) << '\n';

		tomlOutput << toml::format(tbl) << '\n';

		tomlOutput.close();
	}
	
	else if (subcommand == "classify") {

		configReader::read_classify_parser_ struct_ = config.classifyReader(tomlOutput, subcommand, target_files_, deplete_files_);
		read_classify_parser classify_command;

		std::string reads_;
		reads_ = struct_.read_file;
		classify_command.ibf_deplete_file = struct_.ibf_deplete_file;
		classify_command.ibf_target_file = struct_.ibf_target_file;
		classify_command.read_file = reads_;
		classify_command.out_dir = struct_.out_dir;
		classify_command.command = struct_.command;
		classify_command.show_help = struct_.show_help;
		classify_command.kmer_significance = struct_.kmer_significance;
		classify_command.error_rate = struct_.error_rate;
		classify_command.threads = struct_.threads;
		classify_command.preLen = struct_.preLen;
		classify_command.max_chunks = struct_.max_chunks;
		classify_command.verbose = struct_.verbose;

		classify_reads(classify_command);

		auto tbl = toml::value{ {

		{subcommand, toml::table{{
		{ "deplete_files", struct_.ibf_deplete_file  },
		{ "target_files", struct_.ibf_target_file  },
		{ "read_files", reads_},
		{ "exp_seq_error_rate", struct_.kmer_significance },
		{ "threads", struct_.threads },
		{ "chunk_length", struct_.preLen },
		{ "max_chunks", struct_.max_chunks }

		    }}
		   },
		} };

		// chrono: https://en.cppreference.com/w/cpp/chrono
		auto start = std::chrono::system_clock::now();
		auto end = std::chrono::system_clock::now();

		std::chrono::duration<double> elapsed_seconds = end - start;
		std::time_t end_time = std::chrono::system_clock::to_time_t(end);

		tomlOutput << "# Computation time: " << std::ctime(&end_time) << '\n';

		tomlOutput << toml::format(tbl) << '\n';

		tomlOutput.close();

	}

	else if (subcommand == "target") {

		configReader::live_target_parser_ struct_ = config.targetReader(tomlOutput, subcommand, target_files_, deplete_files_);

		target_parser  target_command;

		target_command.host = struct_.host;
		target_command.device = struct_.device;
		target_command.ibf_deplete_file = struct_.ibf_deplete_file;
		target_command.ibf_target_file = struct_.ibf_target_file;
		target_command.output_dir = struct_.output_dir;
		target_command.guppy_host = struct_.guppy_host;
		target_command.guppy_port = struct_.guppy_port;
		target_command.caller = struct_.caller;
		target_command.port = struct_.port;
		target_command.basecall_threads = struct_.basecall_threads;
		target_command.classify_threads = struct_.classify_threads;
		target_command.kmer_significance = struct_.kmer_significance;
		target_command.error_rate = struct_.error_rate;
		target_command.command = struct_.command;
		target_command.show_help = struct_.show_help;
		target_command.verbose = struct_.verbose;

		connection_test_parser cT = { target_command.host, target_command.device, target_command.port, false, false, true, false };
		test_connection(cT);
		

		auto tbl1 = toml::value{ {
		{subcommand, toml::table{{
			{ "flowcell", target_command.device },
			{ "host-ip ", target_command.host },
			{ "port", target_command.port },
			{ "depletion-file", target_command.ibf_deplete_file },
			{ "target-file", target_command.ibf_target_file },
			{ "significance", target_command.kmer_significance },
			{ "error-rate", target_command.error_rate },
			{ "basecall-threads", target_command.basecall_threads },
			{ "classification-th", target_command.classify_threads },
			{ "caller", target_command.caller },
			{ "guppy_host", target_command.guppy_host },
			{ "guppy_port", target_command.guppy_port }
		   }}
		   },
		} };

		// chrono: https://en.cppreference.com/w/cpp/chrono
		auto start = std::chrono::system_clock::now();
		auto end = std::chrono::system_clock::now();

		std::chrono::duration<double> elapsed_seconds = end - start;
		std::time_t end_time = std::chrono::system_clock::to_time_t(end);

		tomlOutput << "# Computation date: " << std::ctime(&end_time) << '\n';
		tomlOutput << toml::format(tbl1) << '\n';
		tomlOutput.close();
		adaptive_sampling(target_command);

	}


}


int main(int argc, char const **argv)
{

	StopClock NanoLiveTime;
	NanoLiveTime.start();

	std::signal(SIGINT, signalHandler);

	std::string binPath = argv[0];
	NanoLiveRoot = binPath.substr(0, binPath.find("bin"));

	std::string const tomlFile = argv[1];
	initializeLogger(tomlFile);

	std::ifstream tomlFileReadBouncer(tomlFile, std::ios_base::binary);
	const auto configuration_ = toml::parse(tomlFileReadBouncer, /*optional -> */ tomlFile);

	auto log_file = toml::find<std::string>(configuration_, "log_directory");
	auto output_fileTOML = toml::find<std::string>(configuration_, "output_directory");

	configReader config(tomlFile);

	//log files
	readuntil::CSVFile += output_fileTOML;
	interleave::InterleavedBloomFilterLog += log_file;
	interleave::IbfClassificationLog += log_file;
	readuntil::ReadUntilClientLog += log_file;
	
	std::string subcommand = config.usage();
	std::fstream tomlOutput = config.writeTOML();

	if (subcommand.length() > 1) {

		std::cout << "The usage is: " << subcommand << '\n';
		std::cout << "\n";
	}

	else {

		std::cerr << "No usage found in config.TOML file\nPlease define one of the usages:  [build, target, classify]" << '\n';
		exit(0);
	}

	const auto& IBF = toml::find(configuration_, "IBF");

	const auto  target_files = toml::find<std::string>(IBF, "target_files");
	std::string target_files_ = target_files;

	const auto  deplete_files = toml::find<std::string>(IBF, "deplete_files");
	std::string deplete_files_ = deplete_files;

	configurationReader(config, tomlFile, subcommand, tomlOutput, target_files_, deplete_files_);

	NanoLiveTime.stop();

	size_t peakSize = getPeakRSS();
	int peakSizeMByte = (int)(peakSize / (1024 * 1024));

	std::cout << "Real time : " << NanoLiveTime.elapsed() << " sec" << std::endl;
	std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
	std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;

	return 0;
}

