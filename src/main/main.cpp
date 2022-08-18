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

// Qt
#include <QApplication>
#include <QDebug>
#include <QGuiApplication>
#include <QQmlApplicationEngine>
// Qt local
#include "mainwindow.h"
#include "ibf_mainwindow.h"

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
	client.setRootPath(NanoLiveRoot);

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
			nanolive_logger->error("Could not start streaming signals from device (" + config.MinKNOW_Parsed.flowcell + ")");
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
		nanolive_logger = spdlog::rotating_logger_mt("ReadBouncerLog",  ReadBouncerLog.string() , 1048576 * 5, 100);
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
 * Build or load target/deplete IBF's for classify or target usage
 * @param  config ConfigReader constructor 
 * @param  targetFilter bool to parse target filters/fasta files
 * @param  depleteFilter bool to parse deplete filters/fasta files
 * @return vector of loaded/constructed IBF's 
 */

std::vector<interleave::IBFMeta> getIBF (ConfigReader config, bool targetFilter, bool depleteFilter){

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
					nanolive_logger->error("Error parsing depletion IBF using the following parameters");
					nanolive_logger->error("Depletion IBF file                : " + deplete_file.string());
					nanolive_logger->error("Error message : " + std::string(e.what()));
					nanolive_logger->flush();
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

/**
 * Run ReadBouncer using the provided parameters in config.toml file
 * @param  config ConfigReader constructor 
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
				
				ibf_build_parser params = { out.string(), file.string(), false, false, config.IBF_Parsed.size_k, config.IBF_Parsed.threads, config.IBF_Parsed.fragment_size, 0, true };
                buildIBF(params);
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

				ibf_build_parser params = { out.string(), file.string(), false, false, config.IBF_Parsed.size_k, config.IBF_Parsed.threads, config.IBF_Parsed.fragment_size, 0, true };
                buildIBF(params);
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
			classify_reads(config, getIBF(config, false, true), getIBF(config, true, false));
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
			adaptive_sampling(config, getIBF(config, false, true), getIBF(config, true, false));
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

// QT pushes
// Build IBF
void IBF_mainwindow::on_pushButton_clicked()
{
    slot_control_std();
    StopClock NanoLiveTime;

    if (IBF_mainwindow::k < 10){

        std::string warning_kmer = "The selcted k-mer size is smaller than 10, we will use the default value 13";
        QMessageBox::warning(this , tr("Warning"), QString::fromUtf8(warning_kmer.c_str()));
        IBF_mainwindow::k = 13;
    }

    ibf_build_parser params = { IBF_mainwindow::output_path,
                                IBF_mainwindow::reference_file,
                                false, false,
                                IBF_mainwindow::k,
                                IBF_mainwindow::threads,
                                IBF_mainwindow::fragment_size,
                                IBF_mainwindow::filter_size, true };
    NanoLiveTime.start();
    buildIBF(params);
    NanoLiveTime.stop();
    size_t peakSize = getPeakRSS();
    int peakSizeMByte = (int)(peakSize / (1024 * 1024));

    std::cout<<"--------------------------------------------------------------"<<std::endl;
    std::cout << "Real time : " << NanoLiveTime.elapsed() << " sec" << std::endl;
    std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
    std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;
    std::cout<<"--------------------------------------------------------------"<<std::endl;

}


int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    MainWindow mw;

    /*
    int x = mw.width()*0.7;
    int y = mw.height()*0.7;

    mw.setFixedSize(x, y);*/
    /*QGuiApplication app(argc, argv);

    QCoreApplication::addLibraryPath("./gui/");

    QQmlApplicationEngine engine;
    engine.load(QUrl(QStringLiteral("qrc:/main.qml")));*/


    mw.show();

    return app.exec();
}

/*
int main(int argc, char const **argv)
{
	StopClock NanoLiveTime;
	NanoLiveTime.start();

	std::signal(SIGINT, signalHandler);

	std::string binPath = argv[0];
	NanoLiveRoot = binPath.substr(0, binPath.find("bin"));

	std::string const tomlFile = argv[1];
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

	NanoLiveTime.stop();

	size_t peakSize = getPeakRSS();
	int peakSizeMByte = (int)(peakSize / (1024 * 1024));

	std::cout << "Real time : " << NanoLiveTime.elapsed() << " sec" << std::endl;
	std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
    std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;

	return 0;
}*/

