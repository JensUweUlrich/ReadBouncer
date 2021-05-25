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

// subcommand related functions
#include "ibfbuild.hpp"
#include "classify.hpp"
#include "depletion.hpp"

#include <windows.h>
#include <psapi.h>

#include "mainwindow.h"
#include<QApplication>
#include "ibf_mainwindow.h"
#include "classify_mainwindow.h"
#include "connection_test_mainwindow.h"
#include "live_deplete_mainwindow.h"
#include "QDebugStream.h"

#pragma comment( lib, "psapi.lib" )

std::shared_ptr<spdlog::logger> nanolive_logger;

static const unsigned __int64 epoch = ((unsigned __int64)116444736000000000ULL);

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

void signalHandler(int signum)
{
	//TODO: shutdown the gRPC stream smoothly
	if (data != nullptr)
	{
		data->getContext()->TryCancel();
	}
	exit(signum);
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


/*
 *
 *Test connection using Qt
 *
 *
*/
void test_connection_qt(std::string host_name,int port,std::string device_name,
                        bool unblock_all = false)
{
    std::cout << "Trying to connect to MinKNOW" << std::endl;
    std::cout << "Host : " << host_name << std::endl;
    std::cout << "Port : " << port << std::endl;

    std::stringstream sstr;
    sstr << "Port : " << port;

    // create ReadUntilClient object and connect to specified device
    readuntil::ReadUntilClient& client = readuntil::ReadUntilClient::getClient();
    client.setHost(host_name);
    client.setPort(port);

    // TODO: throw exception if connection could not be established
    try
    {
        if (client.connect(device_name))
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
        //throw;
    }
    catch (readuntil::ReadUntilClientException& e)
    {
        std::cerr << "Could not establish connection to MinKNOW." << std::endl;
        std::cerr << "Please check the given host IP address and TCP port. " << std::endl;
        //throw;
    }

    if (unblock_all)
    {


        // wait until sequencing run has been started
        //if (parser.verbose)
            std::cout << "Waiting for device to start sequencing!" << ::std::endl;

        std::cout << "Please start the sequencing run now!" << ::std::endl;

        readuntil::Acquisition* acq = (readuntil::Acquisition*)client.getMinKnowService(readuntil::MinKnowServiceType::ACQUISITION);
        if (acq->hasStarted())
        {
            //if (parser.verbose)
                //std::cout << "Sequencing has begun. Starting live signal processing!" << ::std::endl;
                //std::cout <<"Sequencing has begun. Starting live signal processing!"<<std::endl;
            nanolive_logger->info("Sequencing has begun. Starting live signal processing!");
            nanolive_logger->flush();

        }

        // create Data Service object
        // used for streaming live nanopore signals from MinKNOW and sending action messages back
        data = (readuntil::Data*)client.getMinKnowService(readuntil::MinKnowServiceType::DATA);

        // set unblock all reads

        //(*data).setUnblockAll(true);
        //std::cout<<"Unblocking all reads without basecalling or classification!"<<std::endl;
        nanolive_logger->info("Unblocking all reads without basecalling or classification!");
        nanolive_logger->flush();

        // start live streaming of data
        try
        {
            data->startLiveStream();
        }
        catch (readuntil::DataServiceException& e)
        {
            //std::cerr<<"Could not start streaming signals from device ("<<device_name<<")"<<std::endl;
            //std::cerr<<"Error: " << std::string(e.what())<<std::endl;
            nanolive_logger->error("Could not start streaming signals from device (" + device_name + ")");
            nanolive_logger->error("Error message : " + std::string(e.what()));
            nanolive_logger->flush();
           // throw;
        }



        // thread safe queue storing reads ready for basecalling
        SafeQueue<readuntil::SignalRead> read_queue{};
        // thread safe queue storing classified reads ready for action creation
        SafeQueue<readuntil::ActionResponse> action_queue{};
        // thread safe queue storing for every read the duration for the different tasks to complete
        SafeQueue<Durations> duration_queue{};

        // start live signal streaming from ONT MinKNOW
        std::vector< std::future< void > > tasks;

       // if (parser.verbose)
        //{
          //  std::cout << "Start receiving live signals thread" << std::endl;
           // std::cout << "Start sending unblock messages thread" << std::endl;
        //}


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
// Qt6.0.3

//Argument: IBF Build Step!!!





void IBF_mainwindow::on_pushButton_3_clicked()
{
    StopClock NanoLiveTime;

    QMessageBox::StandardButton ask;
     ask = QMessageBox::question(this, "Building IBF", "The IBF will be built, are you sure?",
                                   QMessageBox::Yes|QMessageBox::No);

     if (ask == QMessageBox::Yes) {
         third_party();
         //proBar();
         NanoLiveTime.start();
         buildIBF_qt(k, 0, f, t, ref_file_Name, output_file_Name);
         NanoLiveTime.stop();
         size_t peakSize = getPeakRSS();
         int peakSizeMByte = (int)(peakSize / (1024 * 1024));
         std::cout<<"--------------------------------------------------------------"<<std::endl;
         std::cout << "Real time : " << NanoLiveTime.elapsed() << " sec" << std::endl;
         std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
         std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;
         std::cout<<"--------------------------------------------------------------"<<std::endl;
         QMessageBox::StandardButton ask_1;
          ask_1 = QMessageBox::question(this, "Clear Results", "Do you want to keep the last results in the output window?",
                                        QMessageBox::Yes|QMessageBox::No);

          if (ask_1 == QMessageBox::Yes) {}

          if (ask_1 == QMessageBox::No) {clearResults();}

          }
}


//Argument: Classify
/*void Classify_mainwindow::on_pushButton_4_clicked()
{
    if (ibf_deplete_file.length() < 1 && ibf_target_file.length() < 1)
    {
        QMessageBox::information(this, "Warning", "Please provide an IBF file for depletion and/or  an IBF file for targeted sequencing. ");
        return;
    }
    QMessageBox::StandardButton ask;
     ask = QMessageBox::question(this, "Classify", "Nanopore Reads will be classified, are you sure?",
                                   QMessageBox::Yes|QMessageBox::No);

     if (ask == QMessageBox::Yes) {
        third_party();
        classify_reads_qt(ibf_deplete_file_name, ibf_target_file_name, kmer_significance,
                          error_rate, classified_file_name, unclassified_file_name,
                          read_file_name, preLen, threads);
         //close();
     } else {}
}*/

void Classify_mainwindow::on_pushButton_4_clicked()
{
    StopClock NanoLiveTime;
    if (ibf_deplete_file.length() < 1 && ibf_target_file.length() < 1)
    {
        QMessageBox::information(this, "Warning", "Please provide an IBF file for depletion and/or  an IBF file for targeted sequencing. ");
        return;
    }
    QMessageBox::StandardButton ask;
     ask = QMessageBox::question(this, "Classify", "Nanopore Reads will be classified, are you sure?",
                                   QMessageBox::Yes|QMessageBox::No);

     if (ask == QMessageBox::Yes) {
         //std::string classified {"class_new.fasta"};
         //std::string unclassified {"unclass_new.fasta"};
        third_party();
        NanoLiveTime.start();
        read_classify_parser classify_parser{ibf_deplete_file_name, ibf_target_file_name, read_file_name, classified,
                                            unclassified, false, false, kmer_significance, error_rate, threads,
                                            preLen, true};
        classify_reads(classify_parser);
        NanoLiveTime.stop();
        size_t peakSize = getPeakRSS();
        int peakSizeMByte = (int)(peakSize / (1024 * 1024));
        std::cout<<"--------------------------------------------------------------"<<std::endl;
        std::cout << "Real time : " << NanoLiveTime.elapsed() << " sec" << std::endl;
        std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
        std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;
        std::cout<<"--------------------------------------------------------------"<<std::endl;
        QMessageBox::StandardButton ask_1;
         ask_1 = QMessageBox::question(this, "Clear Results", "Do you want to keep the last results in the output window?",
                                       QMessageBox::Yes|QMessageBox::No);

         if (ask_1 == QMessageBox::Yes) {}

         if (ask_1 == QMessageBox::No) {clearResults();}

         //close();
     } else {}
}


//Argument:live-deplete

/*void live_deplete_mainwindow::on_pushButton_4_clicked()
{
    if (ibf_deplete_file.length() < 1 && ibf_target_file.length() < 1)
    {
        QMessageBox::information(this, "Warning", "Please provide an IBF file for depletion and/or  an IBF file for targeted sequencing. ");
        return;
    }
    QMessageBox::StandardButton ask;
     ask = QMessageBox::question(this, "Live-Deplete", "Live classification and rejection of nanopore reads will start, Are you sure?",
                                   QMessageBox::Yes|QMessageBox::No);

     if (ask == QMessageBox::Yes) {
        third_party();
        live_read_depletion_qt(ibf_deplete_file_name, ibf_target_file_name, host_name,
                               port, device_name, weights_name, kmer_significance, error_rate);
         //close();
     } else {}
}*/



void live_deplete_mainwindow::on_pushButton_4_clicked()
{
    StopClock NanoLiveTime;
    if (ibf_deplete_file.length() < 1 && ibf_target_file.length() < 1)
    {
        QMessageBox::information(this, "Warning", "Please provide an IBF file for depletion and/or  an IBF file for targeted sequencing. ");
        return;
    }
    QMessageBox::StandardButton ask;
     ask = QMessageBox::question(this, "Live-Deplete", "Live classification and rejection of nanopore reads will start, Are you sure?",
                                   QMessageBox::Yes|QMessageBox::No);

     if (ask == QMessageBox::Yes) {
        live_depletion_parser deplete_parser{host_name, device_name, ibf_deplete_file_name,
                                            ibf_target_file_name, weights_name, port,
                                            kmer_significance, error_rate, false, false, true};

        third_party();
        NanoLiveTime.start();
        live_read_depletion(deplete_parser);
        NanoLiveTime.stop();
        size_t peakSize = getPeakRSS();
        int peakSizeMByte = (int)(peakSize / (1024 * 1024));
        std::cout<<"--------------------------------------------------------------"<<std::endl;
        std::cout << "Real time : " << NanoLiveTime.elapsed() << " sec" << std::endl;
        std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
        std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;
        std::cout<<"--------------------------------------------------------------"<<std::endl;
        /*live_read_depletion_qt(ibf_deplete_file_name, ibf_target_file_name, host_name,
                               port, device_name, weights_name, kmer_significance, error_rate);*/
         //close();
        QMessageBox::StandardButton ask_1;
         ask_1 = QMessageBox::question(this, "Clear Results", "Do you want to keep the last results in the output window?",
                                       QMessageBox::Yes|QMessageBox::No);

         if (ask_1 == QMessageBox::Yes) {}

         if (ask_1 == QMessageBox::No) {clearResults();}
     } else {}
}


//Argument:connection_test

void connection_test_mainwindow::on_pushButton_4_clicked()
{
    StopClock NanoLiveTime;
    QMessageBox::StandardButton ask;
     ask = QMessageBox::question(this, "Test Connection", "Test connection to a working MinKNOW instance, Are you sure?",
                                   QMessageBox::Yes|QMessageBox::No);

     if (ask == QMessageBox::Yes) {
        third_party();
        NanoLiveTime.start();
        test_connection_qt(host_name, port, device_name, unblock_all);
        NanoLiveTime.stop();
        size_t peakSize = getPeakRSS();
        int peakSizeMByte = (int)(peakSize / (1024 * 1024));
        std::cout<<"--------------------------------------------------------------"<<std::endl;
        std::cout << "Real time : " << NanoLiveTime.elapsed() << " sec" << std::endl;
        std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
        std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;
        std::cout<<"--------------------------------------------------------------"<<std::endl;
        QMessageBox::StandardButton ask_1;
         ask_1 = QMessageBox::question(this, "Clear Results", "Do you want to keep the last results in the output window?",
                                       QMessageBox::Yes|QMessageBox::No);

         if (ask_1 == QMessageBox::Yes) {}

         if (ask_1 == QMessageBox::No) {clearResults();}
     } else {}
}

/*void connection_test_mainwindow::on_pushButton_4_clicked()
{
    QMessageBox::StandardButton ask;
     ask = QMessageBox::question(this, "Test Connection", "Test connection to a working MinKNOW instance, Are you sure?",
                                   QMessageBox::Yes|QMessageBox::No);

     if (ask == QMessageBox::Yes) {
         third_party();
         //auto cli = lyra::cli();
         connection_test_parser connect_parser{host_name, device_name,
                                              port, false, false, false, unblock_all};

         test_connection(connect_parser);
        //test_connection_qt(host_name, port, device_name, unblock_all);
     } else {}
}
*/
int main(int argc, char *argv[])
{
    StopClock NanoLiveTime;
    NanoLiveTime.start();
	std::signal(SIGINT, signalHandler);	

	initializeLogger();
	std::string binPath = argv[0];
	NanoLiveRoot = binPath.substr(0, binPath.find("bin"));

    /*auto cli = lyra::cli();
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
    NanoLiveTime.stop();
	size_t peakSize = getPeakRSS();
	int peakSizeMByte = (int)(peakSize / (1024 * 1024));

	std::cout << "Real time : " << NanoLiveTime.elapsed() << " sec" << std::endl;
	std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
    std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;*/



    QApplication a(argc, argv);
    MainWindow w;
    //w.QWidget::showMaximized();
    //w.resize(800, 120);
    w.show();
    return a.exec();
}

