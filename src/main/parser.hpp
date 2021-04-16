#pragma once

#include <lyra/lyra.hpp>

#ifndef PARSER_HPP_
#define PARSER_HPP_

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
	int fragment_size = 100000;
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
				[this](const lyra::group& g) { this->do_command(g); })
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
				.help("Output file of Interleaved Bloom Filter (required)"))
			.add_argument(
				lyra::opt(reference_file, "input-reference")
				.name("-i")
				.name("--input-reference")
				.required()
				.help("Reference sequence file (fasta format) used to build the IBF (required)"))
			.add_argument(
				lyra::opt(size_k, "kmer-size")
				.name("-k")
				.name("--kmer-size")
				.optional()
				.help("Kmer size used for building the Interleaved Bloom Filter (default: 13)"))
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
				.help("Length of fragments from the reference that are put in one bin of the IBF (default: 100000)"))
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
	void do_command(const lyra::group& g)
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
	std::string ibf_deplete_file{ };
	std::string ibf_target_file{ };
	std::string read_file{};
	std::string classified_file{};
	std::string unclassified_file{};
	bool command = false;
	bool show_help = false;
	double kmer_significance = 0.95;
	double error_rate = 0.1;
	int threads = 1;
	int preLen = 360;
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
				[this](const lyra::group& g) { this->do_command(g); })
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
				.help("File with reads to classify in FASTA or FASTQ format (required)"))
			.add_argument(
				lyra::opt(ibf_deplete_file, "ibf-file")
				.name("-d")
				.name("--depletion-file")
				.optional()
				.help("Interleaved Bloom Filter file with depletion references"))
			.add_argument(
				lyra::opt(ibf_target_file, "ibf-file")
				.name("-t")
				.name("--target-file")
				.optional()
				.help("Interleaved Bloom Filter file with target references"))
			.add_argument(
				lyra::opt(classified_file, "file")
				.name("-c")
				.name("--classified-file")
				.optional()
				.help("File with classified reads in FASTA format"))
			.add_argument(
				lyra::opt(unclassified_file, "file")
				.name("-u")
				.name("--unclassified-file")
				.optional()
				.help("File with unclassified reads in FASTA format"))
			.add_argument(
				lyra::opt(kmer_significance, "probability")
				.name("-s")
				.name("--significance")
				.optional()
				.help("significance level for confidence interval of number of errorneous kmers (default: 0.95)"))
			.add_argument(
				lyra::opt(error_rate, "err")
				.name("-e")
				.name("--error-rate")
				.optional()
				.help("expected per read sequencing error rate (default: 0.1)"))
			.add_argument(
				lyra::opt(preLen, "length")
				.name("-p")
				.name("--prefix-length")
				.optional()
				.help("Length of read prefix used for classification (default: 360)"))
			.add_argument(
				lyra::opt(threads, "threads")
				.name("-n")
				.name("--num-threads")
				.optional()
				.help("Number of classification threads"))
		);

	}

	/**
		function is called after parsing the group parameters from the command line
		prints the help page or the parameter values if option verbose is set
	*/
	void do_command(const lyra::group& g)
	{
		if (show_help)
			std::cout << g;
		else
		{
			if (ibf_deplete_file.length() < 1 && ibf_target_file.length() < 1)
			{
				std::cerr << "Please provide an IBF file for depletion and/or  an IBF file for targeted sequencing." << std::endl;
				command = false;
				return;
			}
			// trigger for calling the correct function after parsing the group parameters
			command = true;
			if (verbose)
			{
				std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
				std::cout << "Classify Reads                                : " << "verbose=" << (verbose ? "true" : "false") << std::endl;
				std::cout << "Input read file                               : " << read_file << std::endl;
				if (ibf_deplete_file.length() > 0)
					std::cout << "Depletion IBF file                            : " << ibf_deplete_file << std::endl;
				if (ibf_target_file.length() > 0)
					std::cout << "Target IBF file                               : " << ibf_target_file << std::endl;
				if (classified_file.length() > 0)
					std::cout << "Classified reads file                         : " << classified_file << std::endl;
				if (unclassified_file.length() > 0)
					std::cout << "Unclassified reads file                       : " << unclassified_file << std::endl;
				std::cout << "Significance level for confidence interval    : " << kmer_significance << std::endl;
				std::cout << "Expected sequencing error rate                : " << error_rate << std::endl;
				std::cout << "Length of read prefix used for classification : " << preLen << std::endl;
				std::cout << "Classification threads                            : " << threads << std::endl;
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
	std::string device{};
	std::string ibf_deplete_file{ };
	std::string ibf_target_file{ };
	std::string weights = "48";
	int port = 9501;
	int basecall_threads = 1;
	int classify_threads = 1;
	double kmer_significance = 0.95;
	double error_rate = 0.1;
	bool command = false;
	bool show_help = false;
	bool verbose = false;

	/**
		parser constructor
		creates the live-deplete group and adds it to the lyra cli
		@cli: lyra command line interface object
	*/
	live_depletion_parser(lyra::cli& cli)
	{
		cli.add_argument(
			lyra::command("live-deplete",
				[this](const lyra::group& g) { this->do_command(g); })
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
				.name("-f")
				.name("--flowcell")
				.required()
				.help("Device or FlowCell name for live analysis (required)"))
			.add_argument(
				lyra::opt(host, "ip")
				.name("-i")
				.name("--host-ip")
				.optional()
				.help("IP address on which MinKNOW software runs (default: localhost)"))
			.add_argument(
				lyra::opt(port, "port")
				.name("-p")
				.name("--port")
				.optional()
				.help("MinKNOW communication port (default: 9501)"))
			.add_argument(
				lyra::opt(ibf_deplete_file, "ibf-file")
				.name("-d")
				.name("--depletion-file")
				.optional()
				.help("Interleaved Bloom Filter file with depletion references"))
			.add_argument(
				lyra::opt(ibf_target_file, "ibf-file")
				.name("-t")
				.name("--target-file")
				.optional()
				.help("Interleaved Bloom Filter file with target references"))
			.add_argument(
				lyra::opt(kmer_significance, "probability")
				.name("-s")
				.name("--significance")
				.optional()
				.help("significance level for confidence interval of number of errorneous kmers (default: 0.95)"))
			.add_argument(
				lyra::opt(error_rate, "err")
				.name("-e")
				.name("--error-rate")
				.optional()
				.help("expected per read sequencing error rate (default: 0.1)"))
			.add_argument(
				lyra::opt(weights, "weights")
				.name("-w")
				.name("--weights")
				.optional()
				.choices("48", "56", "64", "80", "96", "256")
				.help("Deep Nano Weights (default: 48; other choices: 56, 64, 80, 96, 256)"))
			.add_argument(
				lyra::opt(basecall_threads, "t")
				.name("-b")
				.name("--basecall-threads")
				.optional()
				.help("Number of threads used for base calling (default: 1)"))
			.add_argument(
				lyra::opt(classify_threads, "t")
				.name("-c")
				.name("--classification-threads")
				.optional()
				.help("Number of threads used for read classification (default: 1)"))
		);

	}

	/**
		function is called after parsing the group parameters from the command line
		prints the help page or the parameter values if option verbose is set
	*/
	void do_command(const lyra::group& g)
	{
		if (show_help)
			std::cout << g;
		else
		{
			if (ibf_deplete_file.length() < 1 && ibf_target_file.length() < 1)
			{
				std::cerr << "Please provide an IBF file for depletion and/or  an IBF file for targeted sequencing." << std::endl;
				command = false;
				return;
			}
			// trigger for calling the correct function after parsing the group parameters
			command = true;
			if (verbose)
			{
				std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
				std::cout << "Live Nanopore Read Depletion                 : " << "verbose=" << (verbose ? "true" : "false") << std::endl;
				std::cout << "Host IP address                              : " << host << std::endl;
				std::cout << "MinKNOW communication port                   : " << port << std::endl;
				std::cout << "Device or Flowcell name                      : " << device << std::endl;
				if (ibf_deplete_file.length() > 0)
					std::cout << "Depletion IBF file                           : " << ibf_deplete_file << std::endl;
				if (ibf_target_file.length() > 0)
					std::cout << "Target IBF file                              : " << ibf_target_file << std::endl;
				std::cout << "Significance level for confidence interval   : " << kmer_significance << std::endl;
				std::cout << "Expected sequencing error rate               : " << error_rate << std::endl;
				std::cout << "Deep Nano Weights for Live Basecalling       : " << weights << std::endl;
				std::cout << "Base calling threads                         : " << basecall_threads << std::endl;
				std::cout << "Classification threads                       : " << classify_threads << std::endl;
				std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
			}
		}
	}
};

/**
	class for generating the IBF build parser group
*/
struct connection_test_parser
{
	// default host & port to communicate with MinKNOW
	std::string host = "127.0.0.1";
	std::string device{};
	int port = 9501;
	bool command = false;
	bool show_help = false;
	bool verbose = false;
	bool unblock_all = false;

	/**
		parser constructor
		creates the connection-test group and adds it to the lyra cli
		@cli: lyra command line interface object
	*/
	connection_test_parser(lyra::cli& cli)
	{
		cli.add_argument(
			lyra::command("connection-test",
				[this](const lyra::group& g) { this->do_command(g); })
			.help("Test connection to a working MinKNOW instance")
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
				.help("Device or FlowCell name for live analysis (required)"))
			.add_argument(
				lyra::opt(host, "host")
				.name("-c")
				.name("--host")
				.optional()
				.help("IP address on which MinKNOW software runs (default: localhost)"))
			.add_argument(
				lyra::opt(port, "port")
				.name("-p")
				.name("--port")
				.optional()
				.help("MinKNOW communication port (default: 9501)"))
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
	void do_command(const lyra::group& g)
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
				std::cout << "NanoLIVE to MinKNOW connection test          : " << "verbose=" << (verbose ? "true" : "false") << std::endl;
				std::cout << "Host IP address                              : " << host << std::endl;
				std::cout << "MinKNOW communication port                   : " << port << std::endl;
				std::cout << "Device or Flowcell name                      : " << device << std::endl;
				std::cout << "Unblock all live reads                       : " << (unblock_all ? "yes" : "no") << std::endl;
				std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
			}
		}
	}
};


#endif //PARSER_HPP_