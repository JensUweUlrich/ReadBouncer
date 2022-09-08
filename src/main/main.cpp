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
#include <ReadBouncerExceptions.hpp>

// spdlog library
#include "spdlog/spdlog.h"
#include "spdlog/sinks/rotating_file_sink.h"

// ReadUntil library
#include "ReadUntilClient.hpp"
#include "Data.hpp"
// IBF library
#include "IBF.hpp"

// tomel parser
#include "configReader.hpp"


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
#include "connection_test.hpp"

#if defined(_WIN32)
	#include <windows.h>
	#include <psapi.h>
	#pragma comment( lib, "psapi.lib" )
	
	static const unsigned __int64 epoch = ((unsigned __int64)116444736000000000ULL);
#else
	#include <sys/resource.h>
	#include <sys/time.h>
#endif


//std::filesystem::path ReadBouncerRoot;
//readuntil::Data* data;


/**
 * Close open connection if [ctrl + c] is performed
 */

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

void initializeLogger(ConfigReader config)
{	
	try
	{

		readuntil::CSVFile = std::filesystem::path(config.log_dir);
		interleave::InterleavedBloomFilterLog = std::filesystem::path(config.log_dir);
		interleave::IbfClassificationLog = std::filesystem::path(config.log_dir);
		readuntil::ReadUntilClientLog = std::filesystem::path(config.log_dir);
		std::filesystem::path ReadBouncerLog (config.log_dir);

		ReadBouncerLog /= "ReadBouncerLog.txt";
		readbouncer_logger = spdlog::rotating_logger_mt("ReadBouncerLog",  ReadBouncerLog.string() , 1048576 * 5, 100);
		readbouncer_logger->set_level(spdlog::level::debug);
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

/*
std::vector<interleave::IBFMeta> getIBF_old (ConfigReader config, bool targetFilter, bool depleteFilter){

	std::vector<interleave::IBFMeta> DepletionFilters{};
	std::vector<interleave::IBFMeta> TargetFilters{};

	if(depleteFilter){
		// parse depletion IBF if given as parameter
		for (std::filesystem::path deplete_file : config.IBF_Parsed.deplete_files)
		{
			interleave::IBFMeta filter{};
			filter.name = deplete_file.stem().string();
			interleave::IBF tf{};
			interleave::IBFConfig DepleteIBFconfig{};

			if (config.filterException(deplete_file)){
				try
				{
					DepleteIBFconfig.input_filter_file = deplete_file.string();
					interleave::FilterStats stats = tf.load_filter(DepleteIBFconfig);
					filter.filter = std::move(tf.getFilter());
					interleave::print_load_stats(stats);
				}
				catch (interleave::ParseIBFFileException& e)
				{
					readbouncer_logger->error("Error parsing depletion IBF using the following parameters");
					readbouncer_logger->error("Depletion IBF file                : " + deplete_file.string());
					readbouncer_logger->error("Error message : " + std::string(e.what()));
					readbouncer_logger->flush();
					throw;
				}

				DepletionFilters.emplace_back(std::move(filter));
			}
		
		    else
			{
				try
				{
					//ibf_build_parser params;
					std::filesystem::path out = std::filesystem::path(config.output_dir);
					out /= deplete_file.filename();
					out.replace_extension("ibf");
					ibf_build_parser params = { out.string(), deplete_file.string(), false, false, config.IBF_Parsed.size_k, config.IBF_Parsed.threads, config.IBF_Parsed.fragment_size, 0, true };
					filter.filter = buildIBF(params);
					}

				catch (std::out_of_range& e)
				{
					throw ConfigReaderException(e.what());
				}
			DepletionFilters.emplace_back(std::move(filter));
			}
		}
		return DepletionFilters;
	}

	if(targetFilter)
	{
		for (std::filesystem::path target_file : config.IBF_Parsed.target_files)
		{
			interleave::IBFMeta filter{};
			filter.name = target_file.stem().string();
			interleave::IBF tf{};
			interleave::IBFConfig TargetIBFconfig{};
			if (config.filterException(target_file))
			{
				try
				{
					TargetIBFconfig.input_filter_file = target_file.string();
					interleave::FilterStats stats = tf.load_filter(TargetIBFconfig);
					filter.filter = std::move(tf.getFilter());
					interleave::print_load_stats(stats);

				}
				catch (interleave::ParseIBFFileException& e)
				{
					nanolive_logger->error("Error building IBF for target file using the following parameters");
					nanolive_logger->error("Depletion IBF file                : " + target_file.string());
					nanolive_logger->error("Error message : " + std::string(e.what()));
					nanolive_logger->flush();
					throw;
				}

				TargetFilters.emplace_back(std::move(filter));
			}
		
			else
			{
				try
				{
					//ibf_build_parser params;
					std::filesystem::path out = std::filesystem::path(config.output_dir);
					out /= target_file.filename();
					out.replace_extension("ibf");
					ibf_build_parser params = { out.string(), target_file.string(), false, false, config.IBF_Parsed.size_k, config.IBF_Parsed.threads, config.IBF_Parsed.fragment_size, 0, true };
					filter.filter = buildIBF(params);
				}

				catch (std::out_of_range& e)
				{
					throw ConfigReaderException(e.what());
				}

				TargetFilters.emplace_back(std::move(filter));
			}
		}

		return TargetFilters;
	}

}
*/

/**
 * Run ReadBouncer using the provided parameters in config.toml file
 * @param  config ConfigReader object 
 */

void run_program(ConfigReader config){

	try
	{
		config.parse(); // parse all params from the different Moduls (one time parse and stores in struct)
	}
	catch (ConfigReaderException& e)
	{
		std::cerr << e.what() << std::endl;
	}
	std::string subcommand = config.usage;

	if (subcommand == "build") {

		//ibf_build_parser params;
		config.createLog(config.usage);
		
		for (std::filesystem::path file : config.IBF_Parsed.target_files)
		{
			if (!std::filesystem::exists(file))
			{
				// TODO: write message in log file
				std::cerr << "[Error] The following target file does not exist: " << file.string() << std::endl;
				return;
			}
			
			if (!config.filterException(file))
			{
				// TODO: write in log file
				std::cout << "The target file: " << file.filename() << " is a fasta file, start building ibf ......." << '\n';
				std::filesystem::path out = std::filesystem::path(config.output_dir);
				out /= file.filename();
				out.replace_extension("ibf");

				buildIBF(config, file.string(), out.string());
				std::cout <<'\n';
			}

			else 
			{
				std::cout<< "[INFO] The following target file is an IBF file: " << file.string() << '\n';
			}
		}

		for (std::filesystem::path file : config.IBF_Parsed.deplete_files)
		{
			if (!std::filesystem::exists(file))
			{
				// TODO: write message in log file
				std::cerr << "[Error] The following deplete file does not exist: " << file.string() << std::endl;
				return;
			}
			
			if (!config.filterException(file))
			{
				// TODO: write in log file
				std::cout << "The deplete file: " << file.filename() << " is a fasta file, start building ibf ......." << '\n';
				std::filesystem::path out = std::filesystem::path(config.output_dir);
				out /= file.filename();
				out.replace_extension("ibf");

				buildIBF(config, file.string(), out.string());
				std::cout <<'\n';
			}
			else 
			{
				std::cout<< "[INFO] The following deplete file is an IBF file: " << file.string() << '\n';
			}
		}

	}
	
	else if (subcommand == "classify") {

		try
		{
			
			config.createLog(config.usage);
			//std::vector<interleave::IBFMeta> DepletionFilters = getIBF(config, false, true);// avoid copying the IBF's 
			//std::vector<interleave::IBFMeta> TargetFilters = getIBF(config, true, false);// avoid copying the IBF's
			classify_reads(config, getIBF(config, true, false), getIBF(config, false, true));
			
		}
		catch(std::exception& e)
		{
		    std::cerr << e.what() << std::endl;
		    return;
		}

	}

	else if (subcommand == "target") {

		try
		{
		    config.createLog(config.usage);
			// adaptive_sampling(config, DepletionFilters, TargetFilters)
			adaptive_sampling(config, getIBF(config, true, false), getIBF(config, false, true));
		}
		catch(std::exception& e)
		{
		    std::cerr << e.what() << std::endl;
		    return;
		}
	}	


	else if( subcommand == "test") 
  {


		try
		{
			config.createLog(config.usage);
		    test_connection(config);
		}

		catch(std::exception& e)
		{
		    std::cerr << e.what() << std::endl;
		    return;
		}
		
		
	}	

	else{

		std::cerr << "Please define one of the usages:  [build, target, classify, test]" << '\n';
		exit(1);
	}

}
int main(int argc, char const **argv)
{

	StopClock ReadBouncerTime;
	ReadBouncerTime.start();

	std::signal(SIGINT, signalHandler);

	std::string binPath = argv[0];
	ReadBouncerRoot = binPath.substr(0, binPath.find("bin"));

	// parse configuration file 
	std::string const tomlFile = parse_config(argc, argv);

	ConfigReader config{};
	try
	{
		config = ConfigReader(tomlFile);
		config.parse_general();
	}
	catch (ConfigReaderException& e)
	{
		std::cerr << "Error in " << tomlFile << std::endl;
		std::cerr << e.what() << std::endl;
		return 1;
	}

	initializeLogger(config);
	run_program(config);

	ReadBouncerTime.stop();

	size_t peakSize = getPeakRSS();
	int peakSizeMByte = (int)(peakSize / (1024 * 1024));

	std::cout << "Real time : " << ReadBouncerTime.elapsed() << " sec" << std::endl;
	std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
	std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;

	return 0;
}

