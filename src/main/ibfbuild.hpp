/*
 * ibfbuild.hpp
 *
 *  Created on: 31.03.2021
 *      Author: jens-uwe.ulrich
 */


std::shared_ptr<spdlog::logger> readbouncer_logger;

/**
* initialize config forand build IBF
* @parser	: input from the command line for "build" command
* @throws : IBFBuildException
*/

/**
 * @todo change the input parameters to build ibf to ConfigReader config, input reference, output ibf
 * 
 */
interleave::TIbf buildIBF(ConfigReader config_reader, const std::string reference_file, const std::string bloom_filter_output_path)
{
	std::shared_ptr<spdlog::logger> readbouncer_logger = spdlog::get("ReadBouncerLog");
	
	interleave::IBFConfig config{};

	config.reference_files.emplace_back(reference_file);
	config.output_filter_file = bloom_filter_output_path;
	config.kmer_size = config_reader.IBF_Parsed.size_k;
	config.threads_build = config_reader.IBF_Parsed.threads;
	config.fragment_length = config_reader.IBF_Parsed.fragment_size;
	//config.filter_size = parser.filter_size;
	//config.verbose = parser.verbose;

	interleave::IBF filter{};
	//interleave::TIbf filter_out;
	try
	{
		interleave::FilterStats stats = filter.create_filter(config);
		interleave::print_build_stats(stats);
	}
	catch (const interleave::IBFBuildException& e)
	{
		readbouncer_logger->error("Error building IBF using the following parameters");
		readbouncer_logger->error("Input reference file                : " + reference_file);
		readbouncer_logger->error("Output IBF file                     : " + bloom_filter_output_path);
		readbouncer_logger->error("Kmer size                           : " + config_reader.IBF_Parsed.size_k);
		readbouncer_logger->error("Size of reference fragments per bin : " + config_reader.IBF_Parsed.fragment_size);
		readbouncer_logger->error("IBF file size in MegaBytes          : " + config.filter_size);
		readbouncer_logger->error("Building threads                    : " + config_reader.IBF_Parsed.threads);
		readbouncer_logger->error("Error message : " + std::string(e.what()));
		readbouncer_logger->error("---------------------------------------------------------------------------------------------------");
		readbouncer_logger->flush();
		throw;
	}

	return filter.getFilter();

}

/**
 * Build or load target/deplete IBF's for classify or target usage
 * @param  config ConfigReader object 
 * @param  targetFilter bool to parse target filters/fasta files
 * @param  depleteFilter bool to parse deplete filters/fasta files
 * @return vector of loaded/constructed IBF's 
 */

std::vector<interleave::IBFMeta> getIBF (ConfigReader config, bool depleteFilter, bool targetFilter){

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
					
					filter.filter = buildIBF(config, deplete_file.string(), out.string());
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
					readbouncer_logger->error("Error building IBF for target file using the following parameters");
					readbouncer_logger->error("Depletion IBF file                : " + target_file.string());
					readbouncer_logger->error("Error message : " + std::string(e.what()));
					readbouncer_logger->flush();
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

					filter.filter = buildIBF(config, target_file.string(), out.string());
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
