/*
 * ibfbuild.hpp
 *
 *  Created on: 31.03.2021
 *      Author: jens-uwe.ulrich
 */


/**
* initialize config forand build IBF
* @parser	: input from the command line for "build" command
* @throws : IBFBuildException
*/

//ibf_build_parser ibfbuild_parser{cli};
//buildIBF(ibfbuild_parser);
void buildIBF(ibf_build_parser & parser)
{
	std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
	interleave::IBFConfig config{};

	config.reference_files.emplace_back(parser.reference_file);
	config.output_filter_file = parser.bloom_filter_output_path;
	config.kmer_size = parser.size_k;
	config.threads_build = parser.threads;
	config.fragment_length = parser.fragment_size;
	config.filter_size = parser.filter_size;
	config.verbose = parser.verbose;

	interleave::IBF filter{};
	try
	{
		interleave::FilterStats stats = filter.create_filter(config);
		interleave::print_stats(stats);
	}
	catch (const interleave::IBFBuildException& e)
	{
		nanolive_logger->error("Error building IBF using the following parameters");
		nanolive_logger->error("Input reference file                : " + parser.reference_file);
		nanolive_logger->error("Output IBF file                     : " + parser.bloom_filter_output_path);
		nanolive_logger->error("Kmer size                           : " + parser.size_k);
		nanolive_logger->error("Size of reference fragments per bin : " + parser.fragment_size);
		nanolive_logger->error("IBF file size in MegaBytes          : " + parser.filter_size);
		nanolive_logger->error("Building threads                    : " + parser.threads);
		nanolive_logger->error("Error message : " + std::string(e.what()));
		nanolive_logger->error("---------------------------------------------------------------------------------------------------");
		nanolive_logger->flush();
		throw;
	}

}



void buildIBF_qt (int k, int s, int f, int t, std::string ref_file_Name, std::string output_file_Name){
    std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
    interleave::IBFConfig config{};

    config.reference_files.emplace_back(ref_file_Name);
    config.output_filter_file =  output_file_Name;
    config.kmer_size = k;
    config.threads_build = t;
    config.fragment_length = f;
    config.filter_size = s;
    //config.verbose = parser.verbose;

    interleave::IBF filter{};
    try
    {
        interleave::FilterStats stats = filter.create_filter(config);
        interleave::print_stats(stats);
    }
    catch (const interleave::IBFBuildException& e)
    {
        nanolive_logger->error("Error building IBF using the following parameters");
        nanolive_logger->error("Input reference file                : " +ref_file_Name);
        nanolive_logger->error("Output IBF file                     : " + output_file_Name);
        nanolive_logger->error("Kmer size                           : " + k);
        nanolive_logger->error("Size of reference fragments per bin : " + f);
        nanolive_logger->error("IBF file size in MegaBytes          : " + s);
        nanolive_logger->error("Building threads                    : " + t);
        nanolive_logger->error("Error message : " + std::string(e.what()));
        nanolive_logger->error("---------------------------------------------------------------------------------------------------");
        nanolive_logger->flush();
        throw;
    }
    //std::cout<< "the new f" << f <<std::endl;
}
