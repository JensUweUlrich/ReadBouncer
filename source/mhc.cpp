#include <string>
#include <vector>

#include "minimizer3.h"
#include "bloomFilter.hpp"

#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/validators.hpp>

#include <seqan3/range/view/kmer_hash.hpp>

#include <seqan3/io/record.hpp>
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
		std::filesystem::path query_read_file
		{ };
		std::string mode
		{ };
		uint8_t size_k{19};
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

	// only for debugging
	// TODO delete after implementing client architecture
	parser.add_option(args.query_read_file, 'q', "query", "query read file");
	std::string k{};
	parser.add_option(k, 'k', "kmer-size", "k-mer size used for bottom up sketching reads", option_spec::ADVANCED, arithmetic_range_validator{1,31});
	if (k.size() > 0)
	{
		args.size_k = stoi(k);
	}
}

void build_ref_bloom_filter(std::vector<uint64_t> &sketch, BloomFilter &bf)
{
	for (uint64_t value : sketch)
	{
		bf.addHashValue(value);
	}

}
/**
 * create a bloom filter from a set of reference sequences
 * @param refFilePaths : vector of file paths
 * @param output : bloom filter output file
 */
void create_bloom_filter(std::vector<std::filesystem::path> &refFilePaths, std::filesystem::path &output)
{
	// parse all provided reference files
	std::vector<std::vector<dna5>> ref_store
	{ };
	for (std::filesystem::path file : refFilePaths)
	{
		// parse all sequences from file
		sequence_file_input fin
		{ file };
		for (auto & record : fin)
		{
			ref_store.push_back(get<field::SEQ>(record));
		}
	}

	// compute minimizer for all reference sequences
	Minimizer minimizer{};
	BloomFilter bf(100000, 200);
	for (std::vector<dna5> ref : ref_store)
	{
		std::vector<uint64_t> sketch = minimizer.getMinimizer(ref);
		build_ref_bloom_filter(sketch, bf);

	}

	// TODO store bloom filter
}

void load_bloom_filter()
{
	// TODO load bloom filter from file
	// TODO store bloom filter in a given data structure
}

void load_query_reads(std::filesystem::path &input, std::vector<dna5_vector> &queries)
{
	sequence_file_input fin{input};
	for (auto & record : fin)
	{
		queries.push_back(get<field::SEQ>(record));
	}
}

void bottom_up_sketching(dna5_vector &read, uint8_t &k)
{
	Minimizer minimizer{};
	std::vector<uint64_t> sketch = minimizer.getMinimizer(read);

}

/**
 * core method to run the program depending on the switched mode given
 * @param : struct of command line arguments provided
 */
void run_program(cmd_arguments &args)
{
	if (std::string("bloom").compare(args.mode) == 0)
	{
		create_bloom_filter(args.sequence_files, args.bloom_filter_output_path);
	}
	else if (std::string("read-until").compare(args.mode) == 0)
	{
		load_bloom_filter();

		// TODO Exchange this method when implementing client architecture for pulling reads from the event sampler
		// method only used for debugging with provided sequence file
		std::vector<dna5_vector> queries{};
		load_query_reads(args.query_read_file, queries);
		// TODO compute bottom up minhash sketch for every read provided
		for (dna5_vector query : queries)
		{
			bottom_up_sketching(query, args.size_k);
		}
		// TODO calculate containment of sketches in reference bloom filter
	}
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

	run_program(args);
	return 0;
}
