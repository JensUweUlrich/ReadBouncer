#include <string>
#include <vector>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/validators.hpp>

#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

using namespace seqan3;

struct cmd_arguments
{
		std::vector<std::filesystem::path> sequence_files
		{ };
		std::filesystem::path bloom_filter_output_path
		{ };
		std::string mode
		{ };
};

void initialize_mhc_argument_parser(argument_parser &parser, cmd_arguments &args)
{
	// TODO refine parser information
	parser.info.author = "Jens-Uwe Ulrich";
	parser.info.short_description = "Reject Nanopore Reads mapping to given reference sequences";
	parser.info.version = "0.0.1";
	parser.info.date = "26-MAR-2019";
	parser.info.email = "ulrichj@rki.de";

	std::vector<std::string> description
	{ "This program receives read data from ONT's ReadUntil API and, decides wether the ",
	  "read is part of a given set of reference sequences and, if so, sends a reject request ", "for the according pore of the read." };
	parser.info.description = description;

	std::vector<std::string> synopsis
	{ "[bloom, read-until] [OPTIONS]" };
	parser.info.synopsis = synopsis;
	//TODO refine examples
	//parser.info.examples = "mhc bloom ";

	parser.add_positional_option(args.mode, "Modus to run mhc : ", value_list_validator(
	{ "bloom", "read-until" }));

	//TODO add all working modes as options and provide all additional information
}

void initialize_bloom_argument_parser(argument_parser &parser, cmd_arguments &args)
{
	// TODO refine parser information
	parser.info.author = "Jens-Uwe Ulrich";
	parser.info.short_description = "Create Bloom Filter of given reference sequences";
	parser.info.version = "0.0.1";
	parser.info.date = "26-MAR-2019";
	parser.info.email = "ulrichj@rki.de";

	parser.add_option(args.bloom_filter_output_path, 'b', "bloom-output", "output file path to bloom filter", option_spec::REQUIRED);
	parser.add_positional_option(args.sequence_files, "reference file(s) to create bloom filter for");
}

void initialize_read_until_argument_parser(argument_parser &parser, cmd_arguments &args)
{
	// TODO refine parser information
	parser.info.author = "Jens-Uwe Ulrich";
	parser.info.short_description = "Create Bloom Filter of given reference sequences";
	parser.info.version = "0.0.1";
	parser.info.date = "26-MAR-2019";
	parser.info.email = "ulrichj@rki.de";
}

int main(int argc, char const ** argv)
{
	argument_parser parser("mhc", argc, argv);
	cmd_arguments args
	{ };
	initialize_mhc_argument_parser(parser, args);
	if (std::string(argv[1]).compare("bloom") == 0)
	{
		args.mode = "bloom";
		argument_parser bloom_parser("mhc bloom", --argc, argv + 1);
		initialize_bloom_argument_parser(bloom_parser, args);

		try
		{
			bloom_parser.parse();
		}
		catch (parser_invalid_argument const & ext)
		{
			std::cerr << "[PARSER ERROR] " << ext.what() << '\n';
			return -1;
		}
	}
	else if (std::string(argv[1]).compare("read-until") == 0)
	{
		args.mode = "read-until";
		argument_parser read_until_parser("mhc read-until", --argc, argv + 1);
		initialize_read_until_argument_parser(read_until_parser, args);

		try
		{
			read_until_parser.parse();
		}
		catch (parser_invalid_argument const & ext)
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
		catch (parser_invalid_argument const & ext)
		{
			std::cerr << "[PARSER ERROR] " << ext.what() << '\n';
			return -1;
		}
	}

	//run_program(args);
	return 0;
}
