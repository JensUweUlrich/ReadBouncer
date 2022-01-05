
/*
 * configReader.hpp
 *
 *  Created on: 22.10.1
 *  
 */

#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include "../toml11/toml.hpp"
#include "IBF.hpp"



#ifdef __cplusplus
extern "C"
{
#endif


    
class configReader {

public:

    // Define many structs as Lyre to keep the reproducibility of Lyra command line
    struct ibf_build_parser_
    {
        std::string bloom_filter_output_path{ };
        std::string reference_file{};
        bool command = false;
        bool show_help = false;
        int size_k = 13;
        int threads = 1;
        int fragment_size = 100000;
        int filter_size = 0;
        bool verbose = false;
    };

    struct read_classify_parser_
    {
        std::string ibf_deplete_file{ };
        std::string ibf_target_file{ };
        std::string read_file{};
        std::string out_dir{};
        bool command = false;
        bool show_help = false;
        double kmer_significance = 0.95;
        double error_rate = 0.1;
        int threads = 1;
        int preLen = 360;
        int max_chunks = 1;
        bool verbose = false;
    };

    struct live_parser_
    {
        std::string host = "127.0.0.1";
        std::string device{};
        std::string ibf_deplete_file{ };
        std::string ibf_target_file{ };
        std::string output_dir{ };
        std::string guppy_host = "127.0.0.1";
        std::string guppy_port = "5555";
        std::string caller = "DeepNano";
        int port = 9501;
        int basecall_threads = 1;
        int classify_threads = 1;
        double kmer_significance = 0.95;
        double error_rate = 0.1;
        bool command = false;
        bool show_help = false;
        bool verbose = false;
    };

    /*struct live_depletion_parser_ : live_parser_
    {};*/
    struct live_target_parser_ : live_parser_
    {};
    

    configReader(std::string const);

    std::string  usage( );

    std::fstream writeTOML();

    bool filterException(std::string file);

    ibf_build_parser_  ibfReader(std::fstream& tomlOutput, std::string usage, std::string target_files_, std::string deplete_files_);

    void buildIBF_(ibf_build_parser_& parser);

    read_classify_parser_  classifyReader(std::fstream& tomlOutput, std::string usage, std::string target_files_, std::string deplete_files_);

    //live_depletion_parser_  depleteReader(std::fstream& tomlOutput, std::string usage, std::string target_files_, std::string deplete_files_);

    live_target_parser_  targetReader(std::fstream& tomlOutput, std::string usage, std::string target_files_, std::string deplete_files_);

   


private:

    std::string tomlFile;
    toml::basic_value<struct toml::discard_comments, std::unordered_map, std::vector> toml;

};


#ifdef __cplusplus
}


#endif
