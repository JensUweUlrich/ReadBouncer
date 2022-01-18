

#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include <exception>
#include "../toml11/toml.hpp"
#include "IBF.hpp"



#ifdef __cplusplus
extern "C"
{
#endif

    class ConfigReaderException : public std::exception
    {
    private:
        std::string error_message
        { };

    public:

        explicit ConfigReaderException() {};

        explicit ConfigReaderException(const std::string& msg) :
            error_message(msg)
        {
        }
        virtual ~ConfigReaderException() throw ()
        {
        }

        virtual const char* what() const throw ()
        {
            return error_message.c_str();
        }
    };

    
class configReader {

public:

    // Define many structs as Lyre to keep the reproducibility of Lyra command line
    struct IBF_Build_Params
    {
        std::filesystem::path bloom_filter_output_path{ };
        std::filesystem::path reference_file{};
        bool command = false;
        bool show_help = false;
        int size_k = 13;
        int threads = 1;
        int fragment_size = 100000;
        int filter_size = 0;
        bool verbose = false;
    };

    struct Classify_Params
    {
        std::vector<std::filesystem::path> ibf_deplete_files{ };
        std::vector<std::filesystem::path> ibf_target_files{ };
        std::vector<std::filesystem::path> read_files{};
        std::filesystem::path out_dir{};
        bool command = false;
        bool show_help = false;
        double kmer_significance = 0.95;
        double error_rate = 0.1;
        int threads = 1;
        int preLen = 360;
        int max_chunks = 1;
        bool verbose = false;
    };

    struct Target_Params
    {
        std::string host = "127.0.0.1";
        std::string device{};
        std::vector<std::filesystem::path> ibf_deplete_files{ };
        std::vector<std::filesystem::path> ibf_target_files{ };
        std::filesystem::path out_dir{};
        std::string guppy_host = "127.0.0.1";
        std::string guppy_port = "5555";
        std::string guppy_config = "dna_r9.4.1_450bps_fast";
        std::string caller = "DeepNano";
        std::string port = "9501";
        int basecall_threads = 3;
        int classify_threads = 3;
        double kmer_significance = 0.95;
        double error_rate = 0.1;
        bool command = false;
        bool show_help = false;
        bool verbose = false;
        uint8_t minChannel = 1;
        uint8_t maxChannel = 512;
    };

    configReader(std::string const);

    std::string  usage( );

    std::fstream writeTOML();

    bool filterException(std::filesystem::path& file);

    IBF_Build_Params ibfReader(std::fstream& tomlOutput, std::string& usage, 
                                 std::vector<std::filesystem::path>& target_files_, 
                                 std::vector<std::filesystem::path>& deplete_files_);

    void buildIBF_(IBF_Build_Params& parser);

    Classify_Params  classifyReader(std::fstream& tomlOutput, std::string& usage, 
                                          std::vector<std::filesystem::path>& target_files_, 
                                          std::vector<std::filesystem::path>& deplete_files_);

    //live_depletion_parser_  depleteReader(std::fstream& tomlOutput, std::string usage, std::string target_files_, std::string deplete_files_);

    Target_Params targetReader(std::fstream& tomlOutput, std::string& usage, 
                                      std::vector<std::filesystem::path>& target_files_, 
                                      std::vector<std::filesystem::path>& deplete_files_);

   


private:

    std::string tomlFile;
    toml::basic_value<struct toml::discard_comments, std::unordered_map, std::vector> toml;

};


#ifdef __cplusplus
}


#endif
