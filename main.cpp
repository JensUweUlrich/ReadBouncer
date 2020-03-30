//Assertion failed: (size() > 0), function front
#include <string>//
#include <vector>
#include <math.h>
#include <chrono>
#include <csignal>
#include <algorithm> //
#include <iostream> // 
#include<iterator> //
#include "customBloomFilter.hpp"
#include "bloom_filter.hpp"
#include "bloomFilterException.hpp"
#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>//
#include "minimizer3.hpp"

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
//using namespace shape;



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
		uint8_t size_k
		{ 31 };
		float error_rate
		{ 0.05 };
};

 

// shape generator 
 
//shape generator as vector -> kein Unterschied ? 
//ideen : 
// 1) unsigned long long 
// 2) 

void shape_generator(unsigned long long n, unsigned long long str[], unsigned long long i)
{

   while (i == n)
 {
    std::vector<unsigned long long> t2 ( str , str + n  ) ;// generiere shapes als Array dann defneire ein Vektor von der Shapes-Array

	
    unsigned long long x = count(t2.begin(), t2.end(), 0);// Zaehle wie viele Fehlerraten 5 0's
	for (std::vector<unsigned long long>::iterator it = t2.begin() ; it != t2.end(); ++it) 
	//seqan3::shape t2;
	//seqan3::shape t2{seqan3::bin_literal{}};
	seqan3::shape t3 {seqan3::bin_literal{*it}};
	//std::cout<<*it<<""
	
	
	for(unsigned long long i = 0 ; i < n; i++ )
         if (x > 5 ){
             return;// tue nix
         }
         else {
			 t2[i];// speicher in Vektor von Shapes
		  seqan3::debug_stream <<t2[i]<<"";
		  
		  	
		   std::vector<unsigned long long>::iterator shapes_vektor;// Declaring iterator to a vector
	       //for (auto shapes_vektor = t2.begin(); shapes_vektor != t2.end(); shapes_vektor++);
		   //auto& shape_g = shapes_vektor;
		 //  seqan3::shape t3 {seqan3::bin_literal{shape_g}};
		   //seqan3::debug_stream <<*shapes_vektor<<"";
		   

		   
		
			
		 
		   //seqan3::shape seqan3::bin_literal shapes_vektor ;  
 
         }

    std::cout<<std::endl;
    return; }
    unsigned long long zeros = 0;
    str[i] = zeros;
    shape_generator(n, str, i + 1);
    unsigned long long ones = 1; 
    str[i] = ones;
    shape_generator(n, str, i + 1);


}

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
	{ "[bloom, minhash] [OPTIONS]" };
	parser.info.synopsis = synopsis;
	//TODO refine examples
	//parser.info.examples = "mhc bloom ";

	parser.add_positional_option(args.mode, "Modus to run NanoLiveTk : ", value_list_validator
	{ "bloom", "minhash"});

	//TODO add all working modes as options and provide all additional information
}

void initialize_bloom_argument_parser(argument_parser &parser, cmd_arguments &args)
{
	// TODO refine parser information
	parser.info.author = "Jens-Uwe Ulrich";
	parser.info.short_description = "Create Bloom Filter of given reference sequences";
	parser.info.version = "0.0.2";
	parser.info.date = "03-March-2020";
	parser.info.email = "jens-uwe.ulrich@hpi.de";

	parser.add_option(args.bloom_filter_output_path, 'b', "bloom-output", "output file path to bloom filter", option_spec::REQUIRED);
	parser.add_option(args.error_rate, 'e', "false-positive-rate", "target false positive rate for bloom filter construction [default: 0.05]");
	parser.add_option(args.size_k, 'k', "kmer-size", "k-mer size used for bottom up sketching reads", option_spec::DEFAULT, arithmetic_range_validator
	{ 1, 31 });
	parser.add_positional_option(args.sequence_files, "reference file(s) to create bloom filter for");
}

void initialize_read_until_argument_parser(argument_parser &parser, cmd_arguments &args)
{
	// TODO refine parser information
	parser.info.author = "Jens-Uwe Ulrich";
	parser.info.short_description = "Try finding reads of a given read set in a bloom filter";
	parser.info.version = "0.0.2";
	parser.info.date = "03-March-2020";
	parser.info.email = "jens-uwe.ulrich@hpi.de";

	// only for debugging
	// TODO delete after implementing client architecture
	parser.add_option(args.query_read_file, 'q', "query", "query read file");
	parser.add_option(args.bloom_filter_output_path, 'b', "bloom-filter", "path to bloom filter file", option_spec::REQUIRED);

}

/**
 * check write access for given filepath
 * @param file : file path to write/create
 * return true if file could be created, false otherwise
 */
bool checkWriteAccessRights(std::filesystem::path &file)
{
	::std::ofstream fout;
	fout.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	try
	{
		fout.open(file, std::ofstream::out);
		fout.write("t", sizeof(unsigned char));
		fout.close();
		std::filesystem::remove(file);
	}
	catch (std::ofstream::failure &e)
	{
		return false;
	}

	return true;
}

/**
 **
 **/
uint64_t computeMinimizer(const std::vector<std::filesystem::path> &refFilePaths, const uint16_t &kMerSize, std::vector<std::vector<uint64_t>> &sketch_vector)
{
	Minimizer minimizer
	{ };
	minimizer.setKmerSize(kMerSize);
	minimizer.setWindowSize(50);

///////////////////////////////////////////////////// Changed //////////////////////////////////////////////
	/*int shape_length = 10;
    int str[shape_length];
	std::vector<int> t2 ( str , str + shape_length ) ;
    shape_generator(shape_length,str,0);*/
	// Korrektur Iterator zum Vektor 
	seqan3::shape t2{seqan3::bin_literal{0b1110101000101001010010011001111}};
	/*seqan3::shape t2
    { 0b1110101000101001010010011001111_shape };*/
	
	

	minimizer.setGappedShape(t2); // setGappedShape WAS: setGappedShape(seqan3::shape& newShape) NOW : setGappedShape(std::vector<int>  )
	 
///////////////////////////////////////////////////// Changed //////////////////////////////////////////////

	uint64_t minimizer_number = 0;
	debug_stream << "start loading references ....\n";
	for (std::filesystem::path file : refFilePaths)
	{
		// load ref sequences and compute minimizer one after another
		sequence_file_input fin
		{ file };
		for (auto &record : fin)
		{
			debug_stream << "compute minimizer for " << get<field::ID>(record) << "\n";

			dna5_vector seq = get<field::SEQ>(record);
			// split references by stretches of N into many sequences
			for (auto v : std::views::split(seq, 'N'_dna5))
			{
				dna4_vector sub = v | seqan3::views::convert<seqan3::dna4> | ranges::to<std::vector<seqan3::dna4>>();
				if (sub.size() >= kMerSize)
				{
					// dna4_vector as parameters for getMinimizer
					std::vector<uint64_t> sketch = minimizer.getMinimizerHashValues(sub);
					minimizer_number += sketch.size();
					sketch_vector.push_back(sketch);
				}
			}

		}
	}
	return minimizer_number;
}

void create_bloom_filter(std::vector<std::filesystem::path> &refFilePaths, std::filesystem::path &output, const float error_rate, uint16_t kMerSize)
{
	if (!checkWriteAccessRights(output))
	{
		std::cerr << "ERROR: No access right to create or write to " << output.u8string() << std::endl;
		return;
	}

	std::vector<std::vector<uint64_t>> sketch_vector
	{ };

	// compute optimal parameters for the bloom filter creation

	bloom_parameters parameters;

	// How many elements roughly do we expect to insert?
	parameters.projected_element_count = computeMinimizer(refFilePaths, kMerSize, sketch_vector);

	// Maximum tolerable false positive probability? (0,1)
	parameters.false_positive_probability = error_rate; // 1 in 10000

	// Simple randomizer (optional)
	parameters.random_seed = 0xA5A5A5A5;

	if (!parameters)
	{
		throw BloomFilterException("Error - Invalid set of bloom filter parameters!");
	}

	parameters.compute_optimal_parameters();

	//Instantiate Bloom Filter
	CustomBloomFilter filter(parameters, kMerSize);
	//std::cout << "Open Bloom Filter size in bits: " << filter.size() << std::endl;
	//std::cout << "Open Bloom Filter hash number: " << filter.hash_count() << std::endl;

	for (std::vector<uint64_t> sketch : sketch_vector)
	{
		filter.insert(sketch.begin(), sketch.end());
	}

	// store kmerSize, hash seeds and bloom filter in a file
	filter.writeToFile(output);
}

bool bottom_up_sketching(dna4_vector &read, CustomBloomFilter &bf)
{
	std::chrono::high_resolution_clock::time_point begin, end;
	begin = std::chrono::high_resolution_clock::now();
	Minimizer minimizer
	{ };
	minimizer.setKmerSize(bf.kMerSize);
	minimizer.setWindowSize(50);

///////////////////////////////////////////////////// Changed //////////////////////////////////////////////
	/*int shape_length = 10;
    int str[shape_length];
    shape_generator(shape_length,str,0);
	std::vector<int> t2 ( str , str + shape_length ) ;*/
	 
	seqan3::shape t2{seqan3::bin_literal{0b1110101000101001010010011001111}};
	minimizer.setGappedShape(t2);
///////////////////////////////////////////////////// Changed //////////////////////////////////////////////



	std::vector<uint64_t> sketch = minimizer.getMinimizerHashValues(read);
	int num_containments
	{ 0 };
	//debug_stream << sketch.size() << "\n";
	for (uint64_t minimizer : sketch)
	{
		//debug_stream << minimizer << " ";
		if (bf.contains(minimizer))
		{
			++num_containments;
		}
	}
	//debug_stream << std::endl;
	end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
	debug_stream << "used time: " << duration << "\n";
	debug_stream << "Number of minimizer Containments: " << num_containments << "/" << sketch.size() << std::endl;
	return (double(num_containments) / double(sketch.size())) > 0.15;
}

/**
 * core method to run the program depending on the switched mode given
 * @param : struct of command line arguments provided
 */


void run_program(cmd_arguments &args)
{
	if (std::string("bloom").compare(args.mode) == 0)
	{
		create_bloom_filter(args.sequence_files, args.bloom_filter_output_path, args.error_rate, args.size_k);
	}
	else if (std::string("minhash").compare(args.mode) == 0)
	{
		debug_stream << "run program..." << std::endl;
		CustomBloomFilter bf
		{ };
		try
		{
			bf.readFromFile(args.bloom_filter_output_path);
		}
		catch (BloomFilterException &ex)
		{
			std::cerr << "ERROR: Failed to read Bloom Filter from " + args.bloom_filter_output_path.u8string() << std::endl;
			std::cerr << ex.what() << std::endl;
			return;
		}
		debug_stream << "bloom filter read" << std::endl;
		// TODO Exchange this method when implementing client architecture for pulling reads from the event sampler
		// method only used for debugging with provided sequence file

		// TODO compute bottom up minhash sketch for every read provided
		int num_contained_reads
		{ 0 };
		int num_query_reads
		{ 0 };
		sequence_file_input fin
		{ args.query_read_file };
		int k = 0;
		for (auto &record : fin)
		{
			num_query_reads++;
			dna4_vector query = get<field::SEQ>(record) | seqan3::views::convert<seqan3::dna4> | ranges::to<std::vector<seqan3::dna4>>();
			for (int i = 1; i <= 3; ++i)
			{
				std::vector<dna4> read(query.begin() + 100 + (i * 500), query.begin() + 100 + ((i + 1) * 500));
				if (bottom_up_sketching(read, bf))
				{
					num_contained_reads++;
					break;
					//	debug_stream << read << std::endl;
				}
				/*if (num_contained_reads > 10)
				 {
				 break;
				 }*/
			}

		}
		debug_stream << "Number of contained reads: " << num_contained_reads << "/" << num_query_reads;
		// TODO calculate containment of sketches in reference bloom filter
	}
}

int main(int argc, char const **argv)
{
	
	
        unsigned long long shape_length = 10;
        unsigned long long str[shape_length];
        shape_generator(shape_length,str,0);
		
		
		
	//seqan3::shape t2 
    //{ 0b1110101000101001010010011001111_shape };



	argument_parser parser("mhc", argc, argv);
	cmd_arguments args
	{ };
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

	run_program(args);
	return 0;
}

