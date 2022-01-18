

#include "configReader.hpp"
#include <string.h>
#include <filesystem>
#include <sstream>
#include <string>
// seqan libraries
#include <seqan/binning_directory.h>
#include <IBF.hpp>
//#include <ranges>


namespace toml
{
    inline namespace literals
    {
        inline namespace toml_literals
        {
            toml::value operator"" _toml(const char* str, std::size_t len);
        } // toml_literals
    } // literals
} // toml

typedef seqan::BinningDirectory< seqan::InterleavedBloomFilter,
    seqan::BDConfig< seqan::Dna5, seqan::Normal, seqan::Uncompressed > >
    TIbf_;


configReader::configReader(std::string const fileName) {

    this->tomlFile = fileName;
};

/**
 * Find the usage from toml file
 * @return: One usage of [build, classify, deplete, target]
 */

std::string configReader::usage() {

    std::ifstream tomlFileReadBouncer(this->tomlFile, std::ios_base::binary);
    assert(tomlFileReadBouncer.good());

    if (tomlFileReadBouncer.is_open()) {

        std::cout << "We could open and read the toml file: " << this->tomlFile << '\n';
    }
    else {

        std::cerr << "We couldn't parse or read the toml file: " << this->tomlFile << '\n';
    }

    const auto ReadBouncer = toml::parse(tomlFileReadBouncer, /*optional -> */ this->tomlFile);

    this->toml = ReadBouncer;

    const auto usage = toml::find<std::string>(ReadBouncer, "usage");

    return usage;

}


/**
 * Write a log from toml file
 * @return: toml file
 */


std::fstream configReader::writeTOML() {

    auto output_fileTOML = toml::find<std::string>(this->toml, "output_directory");
    output_fileTOML += "/configLog.toml";
    std::fstream tomlOutput(output_fileTOML, std::ios::app | std::ios::out | std::ios::in);

    return tomlOutput;
    
}


/**
 * Check if the target/deplete input files are IBF or not
 * @param : File name
 * @return: Bool due to decision
 */

bool configReader::filterException(std::filesystem::path& file) {

    TIbf_ filter;

    try
    {
        seqan::retrieve(filter, seqan::toCString(file.string()));
    }
    catch (seqan::Exception const& e)
    {
        return false;
    }

    return true;
}

// Copied from buildIBF (src/main/ibfbuild.hpp)

void configReader::buildIBF_(IBF_Build_Params& parser)
{
    std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("ReadBouncerLog");
    interleave::IBFConfig config{};

    config.reference_files.emplace_back(parser.reference_file.string());
    config.output_filter_file = parser.bloom_filter_output_path.string();
    config.kmer_size = parser.size_k;
    config.threads_build = parser.threads;
    config.fragment_length = parser.fragment_size;
    config.filter_size = parser.filter_size;
    config.verbose = parser.verbose;

    interleave::IBF filter{};
    try
    {
        interleave::FilterStats stats = filter.create_filter(config);
        interleave::print_build_stats(stats);
    }
    catch (const interleave::IBFBuildException& e)
    {
        nanolive_logger->error("Error building IBF using the following parameters");
        nanolive_logger->error("Input reference file                : " + parser.reference_file.string());
        nanolive_logger->error("Output IBF file                     : " + parser.bloom_filter_output_path.string());
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


/**
 * Parse parameters from toml file to build IBF from the reference sequence(s)
 * @param : Toml output file, usage, a list of target and deplete files 
 * @return: Struct with the needed parameters to construct IBF (all IBF will be constructed within ibfReader)
 */

configReader::IBF_Build_Params configReader::ibfReader(std::fstream& tomlOutput, std::string& usage,
                                                        std::vector<std::filesystem::path>& target_files_, 
                                                        std::vector<std::filesystem::path>& deplete_files_) {
    IBF_Build_Params build_IBF;
    int k, t, f;
    std::filesystem::path output_fileTOML{};
    try
    {
        output_fileTOML = toml::find<std::string>(this->toml, "output_directory");
	output_fileTOML = output_fileTOML.make_preferred();
        toml::value IBF = toml::find(this->toml, "IBF");

        k = toml::find_or<int>(this->toml, "IBF", "kmer_size", 13);
        t = toml::find_or<int>(this->toml, "IBF", "threads", 1);
        f = toml::find_or<int>(this->toml, "IBF", "fragment_size", 100000);
    }
    catch (std::out_of_range& e)
    {
        // TODO: write message in log file
        throw ConfigReaderException(e.what());
    }
    
    for (std::filesystem::path file : target_files_)
    {
        if (!std::filesystem::exists(file))
        {
            // TODO: write message in log file
            throw ConfigReaderException("[Error] The following target file does not exist: " + file.string());
        }

        if (!configReader::filterException(file))
        {

            // TODO: write in log file
            std::cout << "The target file is a fasta file, start building ibf ......." << '\n';

            std::filesystem::path target = std::filesystem::path(output_fileTOML);
            target /= file;
            target.replace_extension("ibf");

            build_IBF = { target, file, false, false, k, t, f, 0, true };
            buildIBF_(build_IBF);

        }

    }

    for (std::filesystem::path file : deplete_files_)
    {
        if (!std::filesystem::exists(file))
        {
            // TODO: write message in log file
            throw ConfigReaderException("[Error] The following target file does not exist: " + file.string());
        }

        if (!configReader::filterException(file))
        {

            // TODO: write in log file
            std::cout << "The deplete file is a fasta file, start building ibf ......." << '\n';


            std::filesystem::path deplete = std::filesystem::path(output_fileTOML);
            deplete /= file;
            deplete.replace_extension("ibf");

            build_IBF = { deplete, file, false, false, k, t, f, 0, true };
            buildIBF_(build_IBF);

        }

    }
    

    return build_IBF;
};


/**
 * Parse parameters from toml file for reads classification
 * @param : Toml output file, usage, a list of target and deplete files 
 * @return: Struct with the needed parameters to classify reads
 */

configReader::Classify_Params configReader::classifyReader(std::fstream& tomlOutput, std::string& usage, 
                                                                  std::vector<std::filesystem::path>& target_files_, 
                                                                  std::vector<std::filesystem::path>& deplete_files_) {
    

    Classify_Params classifyStruct;
    int k, t, f, l, m;
    double e;
    std::vector<std::filesystem::path> read_files{};
    std::filesystem::path output_fileTOML{};
    std::vector<std::string> rf_tmp{};
    try
    {
        std::string out = toml::find<std::string>(this->toml, "output_directory");
        output_fileTOML = std::filesystem::path(out).make_preferred();
        toml::value IBF = toml::find(this->toml, "IBF");

        k = toml::find_or<int>(this->toml, "IBF", "kmer_size", 13);
        t = toml::find_or<int>(this->toml, "IBF","threads", 1);
        f = toml::find_or<int>(this->toml, "IBF", "fragment_size", 100000);
        e = toml::find_or<double>(this->toml, "IBF", "exp_seq_error_rate", 0.1);
        l = toml::find_or<int>(this->toml, "IBF", "chunk_length", 250);
        m = toml::find_or<int>(this->toml, "IBF", "max_chunks", 5);

        rf_tmp = toml::find<std::vector<std::string>>(IBF, "read_files");
    }
    catch (std::out_of_range& e)
    {
        // TODO: write message in log file
        throw ConfigReaderException(e.what());
    }

    for (std::string file : rf_tmp)
    {
        std::filesystem::path rf(file);
	rf = rf.make_preferred();
        if (!std::filesystem::exists(rf))
        {
            // TODO: write message in log file
            throw ConfigReaderException("[Error] The following read file does not exist: " + rf.string());
        }
        else
        {
            read_files.emplace_back(std::move(rf));
        }

    }

    std::vector<std::filesystem::path> target_holder{};
    std::vector<std::filesystem::path> deplete_holder{};

    for (std::filesystem::path file : target_files_)
    {
        if (!std::filesystem::exists(file)) 
        {
            // TODO: write message in log file
            throw ConfigReaderException("[Error] The following target file does not exist: " + file.string());
        }

        if (configReader::filterException(file)) 
        {
            target_holder.emplace_back(std::move(file));// If ibf then the target in the same dir  
        }
        else
        {

            // TODO: write in log file
            std::cout << "The target file is a fasta file, start building ibf ......." << '\n';


            std::filesystem::path target = std::filesystem::path(output_fileTOML.string());
            target /= file.filename();
            target.replace_extension("ibf");

            IBF_Build_Params build_IBF = { target, file, false, false, k, t, f, 0, true };
            buildIBF_(build_IBF);
            target_holder.emplace_back(std::move(target));
        }

    }

    for (std::filesystem::path file : deplete_files_)
    {
        if (!std::filesystem::exists(file))
        {
            // TODO: write message in log file
            throw ConfigReaderException("[Error] The following target file does not exist: " + file.string());
        }

        if (configReader::filterException(file))
        {
            deplete_holder.emplace_back(std::move(file));// If ibf then the target in the same dir  
        }
        else
        {

            // TODO: write in log file
            std::cout << "The deplete file is a fasta file, start building ibf ......." << '\n';


            std::filesystem::path deplete = std::filesystem::path(output_fileTOML);
            deplete /= file.filename();
            deplete.replace_extension("ibf");

            IBF_Build_Params build_IBF = { deplete, file, false, false, k, t, f, 0, true };
            buildIBF_(build_IBF);
            deplete_holder.emplace_back(std::move(deplete));
        }

    }

    classifyStruct = { deplete_holder, target_holder, read_files, output_fileTOML, false, false, 0.95 , e, t, l, m, false };

    return classifyStruct;

 };

    


/**
 * Parse parameters from toml file for live reads targeting
* @param : Toml output file, usage, a list of target and deplete files 
* @return: Struct with the needed parameters for live target ((Targeted Sequencing))
*/
configReader::Target_Params configReader::targetReader(std::fstream& tomlOutput, std::string& usage,
                                                             std::vector<std::filesystem::path>& target_files_, 
                                                             std::vector<std::filesystem::path>& deplete_files_) {

    std::string target, deplete;
    IBF_Build_Params build_IBF;

    int k, classifyThreads, f, l, m, basecallThreads;
    double e;
    double significance = 0.95;
    std::filesystem::path output_fileTOML{};
    std::vector<std::string> rf_tmp{};
    std::vector<int> channels;
    std::string device{};
    std::string MinKNOW_host{};
    std::string MinKNOW_port{};
    std::string weights = "48";
    std::string caller{};
    std::string hostCaller{};
    std::string portBasecaller{};
    std::string guppyConfig{};
    try
    {
        output_fileTOML = std::filesystem::path(toml::find<std::string>(this->toml, "output_directory"));
	output_fileTOML = output_fileTOML.make_preferred();
        toml::value IBF = toml::find(this->toml, "IBF");
        toml::value MinKNOW = toml::find(this->toml, "MinKNOW");
        toml::value basecaller = toml::find(this->toml, "Basecaller");
	
        k = toml::find_or<int>(IBF, "kmer_size", 13);
        classifyThreads = toml::find_or<int>(IBF, "threads", 1);
        f = toml::find_or<int>(IBF, "fragment_size", 100000);
        e = toml::find_or<double>(IBF, "exp_seq_error_rate", 0.1);
	
        device = toml::find<std::string>(MinKNOW, "flowcell");
        MinKNOW_host = toml::find_or<std::string>(MinKNOW, "host", "127.0.0.1");
        MinKNOW_port = toml::find_or<std::string>(MinKNOW, "port", "9501");
        channels = toml::find_or<std::vector<int>>(MinKNOW, "channels", std::vector<int>{});
	
        caller = toml::find_or<std::string>(basecaller, "caller", "DeepNano");
        basecallThreads = toml::find_or<int>(basecaller, "threads", 3);
#if defined(_WIN32)
        if (stricmp(caller.c_str(), "guppy") == 0)
#else
        if (strcasecmp(caller.c_str(), "guppy") == 0)
#endif
        {
            hostCaller = toml::find<std::string>(basecaller, "host");
            portBasecaller = toml::find_or<std::string>(basecaller, "port", "5555");
            guppyConfig = toml::find_or<std::string>(basecaller, "config", "dna_r9.4.1_450bps_fast");
            // TODO: check if guppyConfig is correct configuration file

        }
    }
    catch (std::out_of_range& e)
    {
        // TODO: write message in log file
        throw ConfigReaderException(e.what());
    }

    std::vector<std::filesystem::path> target_holder{};
    std::vector<std::filesystem::path> deplete_holder{};

    for (std::filesystem::path file : target_files_)
    {
        if (!std::filesystem::exists(file))
        {
            // TODO: write message in log file
            throw ConfigReaderException("[Error] The following target file does not exist: " + file.string());
        }

        if (configReader::filterException(file))
        {
            target_holder.emplace_back(std::move(file));// If ibf then the target in the same dir  
        }
        else
        {

            // TODO: write in log file
            std::cout << "The target file is a fasta file, start building ibf ......." << '\n';


            std::filesystem::path target = std::filesystem::path(output_fileTOML);
            target /= file.filename();
            target.replace_extension("ibf");

            build_IBF = { target, file, false, false, k, classifyThreads, f, 0, true };
            buildIBF_(build_IBF);
            target_holder.emplace_back(std::move(target));
        }

    }

    for (std::filesystem::path file : deplete_files_)
    {
        if (!std::filesystem::exists(file))
        {
            // TODO: write message in log file
            throw ConfigReaderException("[Error] The following target file does not exist: " + file.string());
        }

        if (configReader::filterException(file))
        {
            deplete_holder.emplace_back(std::move(file));// If ibf then the target in the same dir  
        }
        else
        {

            // TODO: write in log file
            std::cout << "The deplete file is a fasta file, start building ibf ......." << '\n';

            std::filesystem::path deplete = std::filesystem::path(output_fileTOML);
            deplete /= file.filename();
            deplete.replace_extension("ibf");

            build_IBF = { deplete, file, false, false, k, classifyThreads, f, 0, true };
            buildIBF_(build_IBF);

            deplete_holder.emplace_back(std::move(deplete));
        }

    }

    Target_Params targetStruct;
    if (channels.size() == 2 )
        targetStruct = { MinKNOW_host, device, deplete_holder, target_holder, output_fileTOML, hostCaller, portBasecaller, guppyConfig, caller, MinKNOW_port, basecallThreads, classifyThreads,
                   significance, e, false, false, false, (uint8_t)channels[0], (uint8_t)channels[1] };
    else
        targetStruct = { MinKNOW_host, device, deplete_holder, target_holder, output_fileTOML, hostCaller, portBasecaller, guppyConfig, caller, MinKNOW_port, basecallThreads, classifyThreads,
                    significance, e, false, false, false };   

    return targetStruct;
};

/**
 * Parse parameters from toml file for live reads depletion
* @param : Toml output file, usage, a list of target and deplete files
* @return: Struct with the needed parameters for live deplete
*/
/*
configReader::live_depletion_parser_  configReader::depleteReader(std::fstream& tomlOutput, std::string usage,
                                                                  std::string target_files_, std::string deplete_files_) {


    std::string target, deplete;
    configReader::ibf_build_parser_ build_IBF;
    configReader::live_depletion_parser_ depletionStruct;

    auto output_fileTOML = toml::find<std::string>(this->toml, "output_directory");

    const auto& IBF = toml::find(this->toml, "IBF");
    const auto& MinKNOW = toml::find(this->toml, "MinKNOW");
    const auto& basecaller = toml::find(this->toml, "Basecaller");


    const auto  flowcell = toml::find<std::string>(MinKNOW, "flowcell");
    const auto  hostIP = toml::find<std::string>(MinKNOW, "host");
    toml::value  port_ = toml::get<toml::table >(MinKNOW).at("port");

    std::string device = flowcell;
    std::string MinKNOW_host = hostIP;
    int MinKNOW_port = toml::get<toml::integer>(port_);


    toml::value  threads = toml::get<toml::table >(IBF).at("threads");
    double significance = 0.95;
    double error_rate = toml::find<double     >(IBF, "exp_seq_error_rate");

    const auto  weights_ = "48";
    std::string weights = weights_;
    int classifyThreads = toml::get<toml::integer>(threads);

    const auto  caller_ = toml::find<std::string>(basecaller, "caller");
    toml::value  basecaller_threads = toml::get<toml::table >(basecaller).at("threads");
    const auto  hostCaller_ = toml::find<std::string>(basecaller, "host");
    toml::value  portBasecaller_ = toml::get<toml::table >(basecaller).at("port");

    int basecallThreads = toml::get<toml::integer>(basecaller_threads);
    std::string caller = caller_;
    std::string hostCaller = hostCaller_;
    int portBasecaller = toml::get<toml::integer>(portBasecaller_);


    std::string target_holder, deplete_holder;
    std::stringstream s_stream_1(deplete_files_);

    while (s_stream_1.good()) {

        std::string substr1;
        getline(s_stream_1, substr1, ',');

        if (substr1.length() > 1)
        {
            if (configReader::filterException(substr1)) {

                deplete = deplete_files_ + ",";
                deplete_holder += deplete;// if ibf then the deplete file.ibf is in the same dir

            }
            else if (!configReader::filterException(substr1)) {

                std::cout << "The deplete file is an fasta file, start building ibf ......." << '\n';

                deplete = output_fileTOML + substr1 + "_.ibf";
                deplete_holder = deplete_holder + deplete + ",";

                toml::value  kmerS = toml::get<toml::table >(IBF).at("kmer_size");
                toml::value  fragment_size = toml::get<toml::table >(IBF).at("fragment_size");
                // Load parameters to build IBF from given files
                int f = toml::get<toml::integer>(fragment_size);
                int k = toml::get<toml::integer>(kmerS);

                build_IBF = { deplete, substr1, false, false, k, classifyThreads, f, 0, true };

                buildIBF_(build_IBF);

            }
        }

        else {

            std::cout << "No deplete file found! " << '\n';
            std::cout << "  " << '\n';
            deplete_holder = ",";
        }
    }


    std::stringstream s_stream(target_files_);
    while (s_stream.good()) {

        std::string substr;
        getline(s_stream, substr, ',');

        if (substr.length() > 1)
        {
            if (configReader::filterException(substr)) {

                target = target_files_ + ",";
                target_holder += target;
            }
            else if (!configReader::filterException(substr)) {

                std::cout << "The target file is an fasta file, start building ibf ......." << '\n';

                target = output_fileTOML + substr + "_.ibf";
                target_holder = target_holder + target + ",";

                toml::value  kmerS = toml::get<toml::table >(IBF).at("kmer_size");
                toml::value  fragment_size = toml::get<toml::table >(IBF).at("fragment_size");

                int f = toml::get<toml::integer>(fragment_size);
                int k = toml::get<toml::integer>(kmerS);

                build_IBF = { target, substr, false, false, k, classifyThreads, f, 0, true };

                buildIBF_(build_IBF);

            }
        }

        else {

            std::cout << "No target file found! " << '\n';
            std::cout << "  " << '\n';
            target_holder = ",";
        }
    }


    deplete_holder.pop_back();
    target_holder.pop_back();

    depletionStruct = { MinKNOW_host, device, deplete_holder, target_holder, weights, MinKNOW_port, basecallThreads, classifyThreads, significance, error_rate, false, false, false };

    return depletionStruct;
};
*/
