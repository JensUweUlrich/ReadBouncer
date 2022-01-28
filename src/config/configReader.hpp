

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

    
class ConfigReader {

public:

   toml::basic_value<struct toml::discard_comments, std::unordered_map, std::vector> configuration_ {};
   std::filesystem::path output_dir{};
   std::filesystem::path log_dir{};
   std::string usage; 
   
   struct Target_Params// TODO
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
        uint16_t minChannel = 1;
        uint16_t maxChannel = 512;
    };

   struct IBF_Params 
    {
        int size_k = 13;
        int fragment_size = 100000;
        int threads = 1;
        std::vector<std::filesystem::path> target_files{};
	    std::vector<std::filesystem::path> deplete_files{};
        std::vector<std::filesystem::path> read_files{};
        double error_rate = 0.1;
        int chunk_length = 360;
        int max_chunks = 1;
    }IBF_Parsed;

    struct MinKNOW_Params
    {
        std::string host = "127.0.0.1";
        std::string port = "9501";
        std::string flowcell{};
        uint8_t minChannel = 1;
        uint8_t maxChannel = 512;
    }MinKNOW_Parsed;

    struct Basecaller_Params 
    {
        std::string caller = "DeepNano";
        std::string guppy_host = "127.0.0.1";
        std::string guppy_port = "5555";
        int basecall_threads = 3;
        std::string guppy_config = "dna_r9.4.1_450bps_fast";
    }Basecaller_Parsed;

    ConfigReader(std::string const);

    void parse_general();
    bool filterException(std::filesystem::path& file);
    void parse();
    void createLog(std::string& usage);

    /*Target_Params targetReader( std::string& usage, 
                                      std::vector<std::filesystem::path>& target_files_, 
                                      std::vector<std::filesystem::path>& deplete_files_);*/

private:

    std::string tomlInputFile{};
    void readIBF(), readMinKNOW(), readBasecaller();

};


#ifdef __cplusplus
}


#endif
