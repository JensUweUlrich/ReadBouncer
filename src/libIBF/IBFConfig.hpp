#pragma once

#include <cinttypes>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

//#include <dirent.h>
#include <stdio.h>

// spdlog
#include "spdlog/spdlog.h"
#include "spdlog/sinks/rotating_file_sink.h"

namespace interleave
{

    class DepleteConfig
    {
        public:
            DepleteConfig()
            {
                try
                {
                    classification_logger = spdlog::rotating_logger_mt("ClassifyLog", "logs/IbfClassificationLog.txt", 1048576 * 5, 100);
                }
                catch (const spdlog::spdlog_ex& e)
                {
                    std::cerr << "Classification Log initialization failed: " << e.what() << std::endl;
                }
                classification_logger->set_level(spdlog::level::debug);
                classification_logger->flush_on(spdlog::level::debug);
            }
            ~DepleteConfig() {};
            double      significance;
            double      error_rate;
            uint16_t    max_error;
            uint16_t    strata_filter;
            std::shared_ptr<spdlog::logger> classification_logger;
    };

    class IBFConfig
    {

        public:
            static constexpr uint32_t MBinBits = 8388608;

            std::vector< std::string > reference_files;
            std::string                directory_reference_files = "";
            std::string                extension                 = "";

            
            std::string output_filter_file = "";
            std::string input_filter_file  = "";
            std::string update_filter_file = "";
            bool        update_complete    = false;

            uint64_t filter_size      = 0;
            uint64_t filter_size_bits = 0;

            uint64_t fragment_length  = 0;

            uint16_t kmer_size      = 19;
            uint16_t hash_functions = 4;

            uint16_t threads   = 2;
            uint32_t n_refs    = 400;
            uint32_t n_batches = 500000;

            double   max_fp = 0.01;

            bool     verbose   = false;
            bool     quiet     = false;

            uint16_t threads_build = 1;

            bool hasEnding( std::string const& fullString, std::string const& ending )
            {
                if ( fullString.length() >= ending.length() )
                {
                    return ( 0 == fullString.compare( fullString.length() - ending.length(), ending.length(), ending ) );
                }
                else
                {
                    return false;
                }
            }

            bool validate()
            {
            // both are not mandatory anymore
            /*if ( seqid_bin_file.empty() || output_filter_file.empty() )
                {
                    std::cerr << "--seqid-bin-file and --output-filter-file are mandatory" << std::endl;
                    return false;
                }
            */
                // add references from folder
                // TODO: move this part to global config
            /*    if ( !directory_reference_files.empty() && !extension.empty() )
                {
                    struct dirent* entry = nullptr;
                    DIR*           dp    = nullptr;
                    dp                   = opendir( directory_reference_files.c_str() );
                    if ( dp != nullptr )
                    {
                        while ( ( entry = readdir( dp ) ) )
                        {
                            if ( hasEnding( entry->d_name, extension ) )
                                reference_files.push_back( directory_reference_files + "/" + entry->d_name );
                        }
                    }
                    closedir( dp );
                }
            
                if ( reference_files.empty() )
                {
                    std::cerr << "Please provide reference sequence files with the parameters --reference-files or/and with "
                                "--directory-reference-files and --extension"
                            << std::endl;
                    return false;
                }
            */
                if ( threads <= 2 )
                {
                    threads_build = 1;
                }
                else
                {
                    threads_build = threads - 1; // extra reading thread
                }

                if ( n_batches < 1 )
                    n_batches = 1;

                if ( n_refs < 1 )
                    n_refs = 1;

                // Skip variables if updating, loads from existing filter file
                if ( !update_filter_file.empty() )
                {
                    if ( verbose )
                    {
                        std::cerr << "WARNING: --filter-size[-bits], --kmer-size --hash-funtions ignored, using metadata from "
                                    "--update-filter-file"
                                << std::endl;
                    }

                    kmer_size        = 0;
                    hash_functions   = 0;
                    filter_size      = 0;
                    filter_size_bits = 0;
                }
                
                else
                {   
                    if (filter_size_bits != 0)
                    {
                        filter_size = filter_size_bits / MBinBits;
                    }
                    else
                    {
                        if ( filter_size != 0 )
                        {
                            filter_size_bits = filter_size * MBinBits;
                        }
                    }
                }
                return true;
            }
    };

    inline std::ostream& operator<<( std::ostream& stream, const IBFConfig& config )
    {
        constexpr auto newl{ "\n" };
        constexpr auto separator{ "----------------------------------------------------------------------" };

        stream << separator << newl;
        stream << "--reference-files     " << newl;
        for ( const auto& file : config.reference_files )
        {
            stream << "                      " << file << newl;
        }
        //stream << "--seqid-bin-file      " << config.seqid_bin_file << newl;
        stream << "--output-filter-file  " << config.output_filter_file << newl;
        stream << "--update-filter-file  " << config.update_filter_file << newl;
        stream << "--update-complete     " << config.update_complete << newl;
        stream << "--filter-size         " << config.filter_size << newl;
        stream << "--filter-size-bits    " << config.filter_size_bits << newl;
        stream << "--hash-functions      " << config.hash_functions << newl;
        stream << "--kmer-size           " << config.kmer_size << newl;
        stream << "--threads             " << config.threads << newl;
        stream << "--n-refs              " << config.n_refs << newl;
        stream << "--n-batches           " << config.n_batches << newl;
        stream << "--verbose             " << config.verbose << newl;
        stream << "--quiet               " << config.quiet << newl;
        stream << separator << newl;

        return stream;
    }
    
    inline void logIBFConfig(const IBFConfig& config)
    {
        constexpr auto newl{ "\n" };
        constexpr auto separator{ "----------------------------------------------------------------------" };

        std::shared_ptr<spdlog::logger> logger = spdlog::get("IbfLog");
       
        logger->debug("--reference-files     ");
        for (const auto& file : config.reference_files)
        {
            logger->debug("                      " + file);
        }
        logger->debug("--output-filter-file  " + config.output_filter_file);
        logger->debug("--update-filter-file  " + config.update_filter_file);
        logger->debug("--update-complete     " + config.update_complete);
        logger->debug("--filter-size         " + config.filter_size);
        logger->debug("--filter-size-bits    " + config.filter_size_bits);
        logger->debug("--hash-functions      " + config.hash_functions);
        logger->debug("--kmer-size           " + config.kmer_size);
        logger->debug("--threads             " + config.threads);
        logger->debug("--n-refs              " + config.n_refs);
        logger->debug("--n-batches           " + config.n_batches);
        logger->debug("--verbose             " + config.verbose);
        logger->debug("--quiet               " + config.quiet);
        logger->debug(separator);
        logger->flush();

    }

} // namespace ibf