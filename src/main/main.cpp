#include <string>
#include <vector>
#include <math.h>
#include <chrono>
#include <csignal>
#include <iostream>
#include <fstream>

// ReadUntil library
#include "ReadUntilClient.hpp"
#include "Data.hpp"
// IBF library
#include "Config.hpp"
#include "IBFBuild.hpp"

#include <lyra/lyra.hpp>

/*
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/validators.hpp>

#include <seqan3/io/record.hpp>
#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include <seqan3/core/debug_stream.hpp>

#include <seqan3/std/ranges>

#include <seqan3/range/views/convert.hpp>

using namespace seqan3;
*/
readuntil::Data *data;

struct cmd_arguments
{
		//std::vector<std::filesystem::path> sequence_files{ };
		//std::filesystem::path bloom_filter_output_path{ };
		//std::filesystem::path query_read_file{ };
		std::string mode
		{ };
		uint8_t size_k
		{ 31 };
		float error_rate
		{ 0.05 };
		std::string host{"127.0.0.1"};
		uint16_t port{9501};
		std::string device{};
		uint8_t unblock_channels{2};
		uint8_t unblock_reads{2};
		uint16_t batch_size{512};
		bool unblock_all{false};
		std::string weights{};
};
/*
void initialize_main_argument_parser(argument_parser &parser, cmd_arguments &args)
{
	// TODO refine parser information
	parser.info.author = "Jens-Uwe Ulrich";
	parser.info.short_description = "Reject Nanopore Reads mapping to given reference sequences";
	parser.info.version = "0.0.2";
	parser.info.date = "03-March-2020";
	parser.info.email = "jens-uwe.ulrich@hpi.de";

	std::vector<std::string> description
	{ "This program receives read data from ONT's MinKNOW API and, decides wether the ",
	  "read is part of a given set of reference sequences and, if so, sends a reject request ", "for the according pore of the read." };
	parser.info.description = description;

	std::vector<std::string> synopsis
	{ "[bloom, minhash, unblock, connectiontest] [OPTIONS]" };
	parser.info.synopsis = synopsis;
	//TODO refine examples
	//parser.info.examples = "mhc bloom ";

	parser.add_positional_option(args.mode, "Modus to run NanoLiveTk : ", value_list_validator
	{ "bloom", "minhash", "unblock", "connectiontest" });

	//TODO add all working modes as options and provide all additional information
}



void initialize_unblock_argument_parser(argument_parser &parser, cmd_arguments &args)
{
	// TODO refine parser information
	parser.info.author = "Jens-Uwe Ulrich";
	parser.info.short_description = "unblock reads streamed from a ONT device";
	parser.info.version = "0.0.3";
	parser.info.date = "27-August-2020";
	parser.info.email = "jens-uwe.ulrich@hpi.de";

	// only for debugging
	parser.add_option(args.device, 'd', "device", "device used for unblocking", option_spec::REQUIRED);
	parser.add_option(args.host, 'c', "host", "ip address of machine hosting MinKNOW application", option_spec::DEFAULT);
	parser.add_option(args.port, 'p', "port", "port on which to communicate with host", option_spec::DEFAULT);
	//parser.add_option(args.unblock_channels, 'u',"unblock-channels", "channels to unblock",option_spec::DEFAULT);
	//parser.add_option(args.unblock_reads, 'r',"unblock-reads","unblock every r-th read of an unblock channel",option_spec::DEFAULT);
	//parser.add_option(args.batch_size, 'b', "batch-size", "number of actions send in one response to MinKNOW", option_spec::DEFAULT);
	parser.add_flag(args.unblock_all, 'a', "unblock-all", "unblock all reads in all channels", option_spec::DEFAULT);

}

void initialize_connectiontest_argument_parser(argument_parser &parser, cmd_arguments &args)
{
	// TODO refine parser information
	parser.info.author = "Jens-Uwe Ulrich";
	parser.info.short_description = "start a client to communicate with ONT's MinKNOW software";
	parser.info.version = "0.0.2";
	parser.info.date = "11-March-2020";
	parser.info.email = "jens-uwe.ulrich@hpi.de";

	// only for debugging
	parser.add_option(args.host, 'c', "host", "host IP address", option_spec::DEFAULT);
	parser.add_option(args.port, 'p', "port", "port on which to communicate with host", option_spec::DEFAULT);
	parser.add_option(args.device, 'd', "device", "device used for unblocking", option_spec::REQUIRED);
}








/**
 * core method to run the program depending on the switched mode given
 * @param : struct of command line arguments provided
 */
 
void run_program(cmd_arguments &args)
{
	
	
    readuntil::ReadUntilClient &client = readuntil::ReadUntilClient::getClient();
	client.setHost(args.host);
	client.setPort(args.port);
	if (client.connect(args.device))
	{
		std::cout << "Connection successfully established!" << ::std::endl;
	}
	else
	{
		std::cerr << "Could not establish connection to MinKNOW or MinION device" << std::endl;
	}

	std::cout << "Waiting for device to start sequencing!" << ::std::endl;

	readuntil::Acquisition *acq = (readuntil::Acquisition*) client.getMinKnowService(readuntil::MinKnowServiceType::ACQUISITION);

	if (acq->hasStarted())
	{
		std::cout << "Sequencing has begun. Starting live signal processing!" << ::std::endl;
	}


   	data = (readuntil::Data*) client.getMinKnowService(readuntil::MinKnowServiceType::DATA);

	try
	{
		if(args.unblock_all)
		{
			(*data).setUnblockAll(true);
		}
		else
		{
			(*data).setUnblockChannels(args.unblock_channels);		
			(*data).setUnblockReads(args.unblock_reads);
		}
		//(*data).setActionBatchSize(args.batch_size);
		(*data).getLiveReads(args.weights);
	}
	catch (readuntil::DataServiceException ex)
	{
		std::cerr << "Could not get live reads : " << ex.what() << std::endl;
		if (data != nullptr)
		{
			data->getContext()->TryCancel();
   		}
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

int main(int argc, char const **argv)
{
	std::signal(SIGINT, signalHandler);	

	//argument_parser parser("mhc", argc, argv);
	cmd_arguments args
	{ };

/*
	initialize_main_argument_parser(parser, args);
	if (std::string(argv[1]).compare("bloom") == 0)
	{
		args.mode = "bloom";
		argument_parser bloom_parser("bloom", --argc, argv + 1);
		initialize_bloom_argument_parser(bloom_parser, args);

		try
		{
			bloom_parser.parse();
		}
		catch (parser_invalid_argument const &ext)
		{
			std::cerr << "[PARSER ERROR] " << ext.what() << '\n';
			return -1;
		}
	}
	else if (std::string(argv[1]).compare("minhash") == 0)
	{
		args.mode = "minhash";
		argument_parser read_until_parser("minhash", --argc, argv + 1);
		initialize_read_until_argument_parser(read_until_parser, args);

		try
		{
			read_until_parser.parse();
		}
		catch (parser_invalid_argument const &ext)
		{
			std::cerr << "[PARSER ERROR] " << ext.what() << '\n';
			return -1;
		}
	}
	else if (std::string(argv[1]).compare("unblock") == 0)
	{
		args.mode = "unblock";
		argument_parser client_parser("unblock", --argc, argv + 1);
		initialize_unblock_argument_parser(client_parser, args);

		try
		{
			client_parser.parse();
		}
		catch (parser_invalid_argument const &ext)
		{
			std::cerr << "[PARSER ERROR] " << ext.what() << '\n';
			return -1;
		}
	}
	else if (std::string(argv[1]).compare("connectiontest") == 0)
	{
		args.mode = "connectiontest";
		argument_parser client_parser("connectiontest", --argc, argv + 1);
		initialize_connectiontest_argument_parser(client_parser, args);

		try
		{
			client_parser.parse();
		}
		catch (parser_invalid_argument const &ext)
		{
			std::cerr << "[PARSER ERROR] " << ext.what() << '\n';
			return -1;
		}
	}
	else
	{
		try
		{
			parser.parse();
		}
		catch (parser_invalid_argument const &ext)
		{
			std::cerr << "[PARSER ERROR] " << ext.what() << '\n';
			return -1;
		}
	}
	
*/

	interleave::Config config{};
	config.reference_files.emplace_back("C:\\NanoLiveTk\\testData\\phiX_reference.fasta");
	config.output_filter_file = "C:\\NanoLiveTk\\testData\\phiX.ibf";
	config.fragment_length = 500;
	config.filter_size = 128;
	config.verbose = true;
	interleave::IBF filter {};
	filter.build(config);

	std::cout<<"hier 295"<<std::endl;
	
	auto cli = lyra::cli();
	std::string command;
	bool show_help = false;
	cli.add_argument(lyra::help(show_help));
	cli.add_argument(lyra::opt(args.device, "device")
						.name("-d")
						.name("--device")
						.required()
						.help("Device or flowCell name for live analysis")
						);
	cli.add_argument(lyra::opt(args.host, "host")
						.name("-c")
						.name("--host")
						.optional()
						.help("IP address on which MinKNOW software runs")
						);
	cli.add_argument(lyra::opt(args.port, "port")
						.name("-p")
						.name("--port")
						.optional()
						.help("MinKNOW communication port")
						);
	cli.add_argument(lyra::opt(args.weights, "weights")
						.name("-w")
						.name("--weights")
						.optional()
						.help("Weights file")
						);					

    // Parse the program arguments:
    const auto result = cli.parse({ argc, argv });

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

	run_program(args);

	return 0;
}

