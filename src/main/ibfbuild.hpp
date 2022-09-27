/*
 * ibfbuild.hpp
 *
 *  Created on: 31.03.2021
 *      Author: jens-uwe.ulrich
 */

#include <QMessageBox>
#include <QTimer>
#include "ibf_mainwindow.h"

std::shared_ptr<spdlog::logger> readbouncer_logger;
/**
* initialize config forand build IBF
* @parser	: input from the command line for "build" command
* @throws : IBFBuildException
*/
interleave::TIbf buildIBF(ibf_build_parser & parser)
{
    std::shared_ptr<spdlog::logger> readbouncer_logger = spdlog::get("ReadBouncerLog");
	interleave::IBFConfig config{};

	config.reference_files.emplace_back(parser.reference_file);
	config.output_filter_file = parser.bloom_filter_output_path;
	config.kmer_size = parser.size_k;
	config.threads_build = parser.threads;
	config.fragment_length = parser.fragment_size;
	config.filter_size = parser.filter_size;
	config.verbose = parser.verbose;

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
        readbouncer_logger->error("Input reference file                : " + parser.reference_file);
        readbouncer_logger->error("Output IBF file                     : " + parser.bloom_filter_output_path);
        readbouncer_logger->error("Kmer size                           : " + parser.size_k);
        readbouncer_logger->error("Size of reference fragments per bin : " + parser.fragment_size);
        readbouncer_logger->error("IBF file size in MegaBytes          : " + parser.filter_size);
        readbouncer_logger->error("Building threads                    : " + parser.threads);
        readbouncer_logger->error("Error message : " + std::string(e.what()));
        readbouncer_logger->error("---------------------------------------------------------------------------------------------------");
        readbouncer_logger->flush();

        throw;
    }

	return filter.getFilter();

}


/**
 * Build or load target/deplete IBF's for classify or target usage
 * @param  config ConfigReader constructor
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
                    QMessageBox::critical(NULL , "Error", QString::fromUtf8(e.what()));
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
                    QMessageBox::critical(NULL , "Error", QString::fromUtf8(e.what()));
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
                    QMessageBox::critical(NULL , "Error", QString::fromUtf8(e.what()));
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
                    QMessageBox::critical(NULL , "Error", QString::fromUtf8(e.what()));
                }

                TargetFilters.emplace_back(std::move(filter));
            }
        }
        return TargetFilters;
    }

}
