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
#include <string_view>

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

	

	if (parser.unblock_all)
	{
		
		readuntil::AnalysisConfiguration* an_conf = (readuntil::AnalysisConfiguration*)client.getMinKnowService(readuntil::MinKnowServiceType::ANALYSIS_CONFIGURATION);
		an_conf->set_break_reads_after_seconds(0.4);
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
*	setup global Logger for ReadBouncer
*/
void initializeLogger(const std::string& toml_file)
{
	std::ifstream tomlFileReadBouncer(toml_file, std::ios_base::binary);
	
	try
	{
		const toml::value configuration_ = toml::parse(tomlFileReadBouncer, /*optional -> */ toml_file);
		std::filesystem::path log_file = toml::find<std::string>(configuration_, "log_directory");
		
		log_file = log_file.make_preferred();
		
		if (!std::filesystem::is_directory(log_file) || !std::filesystem::exists(log_file))
		{
			std::filesystem::create_directories(log_file);
		}

		log_file /= "ReadBouncerLog.txt";
		nanolive_logger = spdlog::rotating_logger_mt("ReadBouncerLog",  log_file.string() , 1048576 * 5, 100);
		nanolive_logger->set_level(spdlog::level::debug);
	}
	catch (const toml::exception& e)
	{
		std::cerr << "Could not parse " << toml_file << std::endl;
		std::cerr << e.what() << std::endl;
		std::cerr << "Please check the correct syntax of the TOML file in the ReadBouncer User Guide!" << std::endl;
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
	std::vector<std::filesystem::path>& target_files_,
	std::vector<std::filesystem::path>& deplete_files_)
{

	toml::value target_files(toml::array{});
	for (std::filesystem::path file : target_files_)
		target_files.push_back(file.string());

	toml::value deplete_files(toml::array{});
	for (std::filesystem::path file : deplete_files_)
		deplete_files.push_back(file.string());

	if (subcommand == "build") {

		configReader::IBF_Build_Params struct_ = config.ibfReader(tomlOutput, subcommand, target_files_, deplete_files_);

		// Create a log of toml usage 

		auto tbl = toml::value{ {
		{subcommand, toml::table{{
				{ "target_files", target_files },
				{ "deplete_files", deplete_files },
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

		configReader::Classify_Params struct_{};
		try
		{
			struct_ = config.classifyReader(tomlOutput, subcommand, target_files_, deplete_files_);
		}
		catch (ConfigReaderException& e)
		{
			std::cerr << "Error in reading TOML configuration file!" << std::endl;
			std::cerr << e.what() << std::endl;
			throw;
		}
		classify_reads(struct_);

		toml::value read_files(toml::array{});
		for (std::filesystem::path file : struct_.read_files)
			read_files.push_back(file.string());

		toml::value tbl = toml::value{ {

		{subcommand, toml::table{{
		{ "deplete_files", deplete_files  },
		{ "target_files", target_files  },
		{ "read_files", read_files},
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

		configReader::Target_Params struct_{};
		try
		{
			 struct_ = config.targetReader(tomlOutput, subcommand, target_files_, deplete_files_);
		}
		catch (ConfigReaderException& e)
		{
			std::cerr << "Error in reading TOML configuration file!" << std::endl;
			std::cerr << e.what() << std::endl;
			throw;
		}
		
		// if caller==guppy => test guppy client connection before starting computation

#if defined(_WIN32)
		if (stricmp(struct_.caller.c_str(), "guppy") == 0)
#else
		if (strcasecmp(struct_.caller.c_str(), "guppy") == 0)
#endif
		{
			std::string basecall_host = struct_.guppy_host + ":" + struct_.guppy_port;
			try
			{
				basecall::GuppyBasecaller* caller = new basecall::GuppyBasecaller(basecall_host, struct_.guppy_config);
				caller->disconnect();
			}
			catch (basecall::BasecallerException& e)
			{
				nanolive_logger->error("Failed establishing connection to Guppy basecall server!");
				nanolive_logger->error("Error message : " + std::string(e.what()));
				nanolive_logger->flush();
				std::cerr << "[Error] Could not connect to Guppy Basecall Server!" << std::endl;
				std::cerr << e.what() << std::endl;
				throw;
			}
		}
		//connection_test_parser cT = { struct_.host, struct_.device, struct_.port, false, false, true, false };
		//test_connection(cT);

		auto tbl1 = toml::value{ {
		{subcommand, toml::table{{
			{ "flowcell", struct_.device },
			{ "host-ip ", struct_.host },
			{ "port", struct_.port },
			{ "minChannel", struct_.minChannel},
			{ "maxChannel", struct_.maxChannel},
			{ "depletion-files", deplete_files },
			{ "target-files", target_files },
			{ "significance", struct_.kmer_significance },
			{ "error-rate", struct_.error_rate },
			{ "basecall-threads", struct_.basecall_threads },
			{ "classification-th", struct_.classify_threads },
			{ "caller", struct_.caller },
			{ "guppy_host", struct_.guppy_host },
			{ "guppy_port", struct_.guppy_port },
			{ "guppy_config", struct_.guppy_config }
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
		
		try
		{
		    adaptive_sampling(struct_);
		}
		catch(std::exception& e)
		{
		    std::cerr << e.what() << std::endl;
		    return;
		}

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
	toml::value configuration_{};
	std::filesystem::path log_file{};
	std::filesystem::path output_fileTOML{};
	try
	{
		 configuration_ = toml::parse(tomlFileReadBouncer, /*optional -> */ tomlFile);
		 log_file = toml::find<std::string>(configuration_, "log_directory");
		 log_file = log_file.make_preferred();
		 output_fileTOML = toml::find<std::string>(configuration_, "output_directory");
		 output_fileTOML = output_fileTOML.make_preferred();
	}
	catch (toml::exception& e)
	{
		std::cerr << "Could not parse " << tomlFile << std::endl;
		std::cerr << e.what() << std::endl;
		return 1;
	}
	catch (std::out_of_range& e)
	{
		std::cerr << "Error in " << tomlFile << std::endl;
		std::cerr << e.what() << std::endl;

		return 1;
	}

	if (!std::filesystem::is_directory(output_fileTOML) || !std::filesystem::exists(output_fileTOML))
	{
		std::filesystem::create_directories(output_fileTOML); 
	}
	
	configReader config(tomlFile);

	//log files
	readuntil::CSVFile = std::filesystem::path(output_fileTOML);
	interleave::InterleavedBloomFilterLog = std::filesystem::path(log_file);
	interleave::IbfClassificationLog = std::filesystem::path(log_file);
	readuntil::ReadUntilClientLog = std::filesystem::path(log_file);
	
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

	std::vector<std::filesystem::path> target_files{};
	std::vector<std::filesystem::path> deplete_files{};
	try
	{
		const toml::value& IBF = toml::find(configuration_, "IBF");
		std::vector<std::string> tmp = toml::find<std::vector<std::string>>(IBF, "target_files");
		for (std::string s : tmp)
			target_files.emplace_back((std::filesystem::path(s)).make_preferred());
		tmp.clear();
		tmp = toml::find<std::vector<std::string>>(IBF, "deplete_files");
		for (std::string s : tmp)
			deplete_files.emplace_back((std::filesystem::path(s)).make_preferred()); 
	}
	catch (std::out_of_range& e)
	{
		std::cerr << "Error in " << tomlFile << std::endl;
		std::cerr << e.what() << std::endl;

		return 1;
	}
	
	
	configurationReader(config, tomlFile, subcommand, tomlOutput, target_files, deplete_files);
	
	NanoLiveTime.stop();

	size_t peakSize = getPeakRSS();
	int peakSizeMByte = (int)(peakSize / (1024 * 1024));

	std::cout << "Real time : " << NanoLiveTime.elapsed() << " sec" << std::endl;
	std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
	std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;

	return 0;
}

