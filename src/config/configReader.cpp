

#include "configReader.hpp"
#include <string.h>
#include <filesystem>
#include <sstream>
#include <string>
// seqan libraries
#include <seqan/binning_directory.h>
#include <IBF.hpp>


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
    output_fileTOML += "configLog.toml";
    std::fstream tomlOutput(output_fileTOML, std::ios::app | std::ios::out | std::ios::in);

    return tomlOutput;
    
}


/**
 * Check if the target/deplete input files are IBF or not
 * @param : File name
 * @return: Bool due to decision
 */

bool configReader::filterException(std::string file) {

    TIbf_ filter;

    try
    {
        seqan::retrieve(filter, seqan::toCString(file));
    }
    catch (seqan::Exception const& e)
    {
        return false;
    }

    return true;
}

// Copied from buildIBF (src/main/ibfbuild.hpp)

void configReader::buildIBF_(ibf_build_parser_& parser)
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
        interleave::print_build_stats(stats);
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


/**
 * Parse parameters from toml file to build IBF from the reference sequence(s)
 * @param : Toml output file, usage, a list of target and deplete files 
 * @return: Struct with the needed parameters to construct IBF (all IBF will be constructed within ibfReader)
 */

configReader::ibf_build_parser_ configReader::ibfReader(std::fstream& tomlOutput, std::string usage,
                                                        std::string target_files_, std::string deplete_files_) {
    configReader::ibf_build_parser_ build_IBF;
    configReader::ibf_build_parser_ buildStruct;

    auto output_fileTOML = toml::find<std::string>(this->toml, "output_directory");

    // find IBF module 
    const auto& IBF = toml::find(this->toml, "IBF");

    toml::value  kmerS = toml::get<toml::table >(IBF).at("kmer_size");
    int k = toml::get<toml::integer>(kmerS);

    toml::value  threads = toml::get<toml::table >(IBF).at("threads");
    int t = toml::get<toml::integer>(threads);

    toml::value  fragment_size = toml::get<toml::table >(IBF).at("fragment_size");
    int f = toml::get<toml::integer>(fragment_size);


    std::stringstream s_stream(target_files_);// list of target files 

    while (s_stream.good()) {

        std::string substr;
        getline(s_stream, substr, ','); 

    if (substr.length() > 1)
    {
        if (configReader::filterException(substr)) {}

        else if (!configReader::filterException(substr)) {

            std::cout << "The target file is an fasta file, start building ibf ......." << '\n';
            std::string outputIBF = substr + "_.ibf";
            build_IBF = { output_fileTOML+outputIBF, substr, false, false, k, t, f, 0, true };
            buildIBF_(build_IBF);
        }
    }
    
    else {

        std::cout << "No target file found! " << '\n';
        std::cout << "  " << '\n';
    }
    }

    std::stringstream s_stream1(deplete_files_);// list of depletion files 

    while (s_stream1.good()) {

        std::string substr1;
        getline(s_stream1, substr1, ','); 

    if (substr1.length() > 1)
    {
        if (configReader::filterException(substr1)) {}

        else if (!configReader::filterException(substr1)) {

            std::string outputIBF = substr1 + "_.ibf";
            build_IBF = { output_fileTOML+outputIBF, substr1, false, false, k, t, f, 0, true };

            buildIBF_(build_IBF);

            buildStruct = { outputIBF, substr1, false, false, k, t, f, 0, true };
        }
    }
    else {
        std::cout << "No deplete file found! " << '\n';
        std::cout << "  " << '\n';
    }
    }

    return buildStruct;
};


/**
 * Parse parameters from toml file for reads classification
 * @param : Toml output file, usage, a list of target and deplete files 
 * @return: Struct with the needed parameters to classify reads
 */

configReader::read_classify_parser_  configReader::classifyReader(std::fstream& tomlOutput, std::string usage, 
                                                                  std::string target_files_, std::string deplete_files_) {
    

    configReader::ibf_build_parser_ build_IBF;
    configReader::read_classify_parser_ classifyStruct;

    auto output_fileTOML = toml::find<std::string>(this->toml, "output_directory");

    const auto& IBF = toml::find(this->toml, "IBF");

    toml::value  kmerS = toml::get<toml::table >(IBF).at("kmer_size");
    int k = toml::get<toml::integer>(kmerS);

    toml::value  threads = toml::get<toml::table >(IBF).at("threads");
    int t = toml::get<toml::integer>(threads);

    toml::value  fragment_size = toml::get<toml::table >(IBF).at("fragment_size");
    int f = toml::get<toml::integer>(fragment_size);

    double e = toml::find<double     >(IBF, "exp_seq_error_rate");

    toml::value  chunk_length = toml::get<toml::table >(IBF).at("chunk_length");
    int l = toml::get<toml::integer>(chunk_length);

    toml::value  max_chunks = toml::get<toml::table >(IBF).at("max_chunks");
    int m = toml::get<toml::integer>(max_chunks);

    const auto  read_files = toml::find<std::string>(IBF, "read_files");
    std::string read_files_ = read_files;

    std::string target, deplete;

    std::string target_holder, deplete_holder;

    std::stringstream s_stream(target_files_);

    while (s_stream.good()) {

        std::string substr;
        getline(s_stream, substr, ',');

    if (substr.length() > 1)
    {
        if (configReader::filterException(substr)) {

            target = substr + ",";
            target_holder += target;// If ibf then the target in the same dir
            //target.clear();
        }
        
        else if (!configReader::filterException(substr)) {

            std::cout << "The target file is an fasta file, start building ibf ......." << '\n';

            target = output_fileTOML + substr + "_.ibf";
            target_holder = target_holder + target + ",";
            build_IBF = { target, substr, false, false, k, t, f, 0, true };
            buildIBF_(build_IBF);

        }   
        
    }
    
    else {

        std::cout << "No target file found! " << '\n';
        std::cout << "  " << '\n';
        target_holder = ",";

    }
    }

    std::stringstream s_stream_1(deplete_files_);

    while (s_stream_1.good()) {

        std::string substr1;
        getline(s_stream_1, substr1, ','); 

        if (substr1.length() > 1)
        {
            if (configReader::filterException(substr1)) {

                deplete = substr1 + ",";
                deplete_holder += deplete;// if ibf then the deplete file.ibf is in the same dir

            }
            else if (!configReader::filterException(substr1)) {

                std::cout << "The deplete file is an fasta file, start building ibf ......." << '\n';

                deplete = output_fileTOML + substr1 + "_.ibf";
                deplete_holder = deplete_holder + deplete + ",";
                build_IBF = { deplete, substr1, false, false, k, t, f, 0, true };
                buildIBF_(build_IBF);

            }
        }

        else {

            std::cout << "No deplete file found! " << '\n';
            std::cout << "  " << '\n';
            deplete_holder = ",";
        }
    
    }

 

    deplete_holder.pop_back();
    target_holder.pop_back();

    classifyStruct = { deplete_holder, target_holder, read_files_, output_fileTOML, false, false, 0.95 , e, t, l, m, false };

    return classifyStruct;

 };

    


/**
 * Parse parameters from toml file for live reads targeting
* @param : Toml output file, usage, a list of target and deplete files 
* @return: Struct with the needed parameters for live target ((Targeted Sequencing))
*/
configReader::live_target_parser_ configReader::targetReader(std::fstream& tomlOutput, std::string usage,
                                                             std::string target_files_, std::string deplete_files_) {

    std::string target, deplete;
    configReader::ibf_build_parser_ build_IBF;

    auto output_fileTOML = toml::find<std::string>(this->toml, "output_directory");
    std::string output_dir = output_fileTOML;
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

    // Find caller settings 
    const auto  caller_ = toml::find<std::string>(basecaller, "caller");
    toml::value  basecaller_threads = toml::get<toml::table >(basecaller).at("threads");
    const auto  hostCaller_ = toml::find<std::string>(basecaller, "host");
    toml::value  portBasecaller_ = toml::get<toml::table >(basecaller).at("port");
    int basecallThreads = toml::get<toml::integer>(basecaller_threads);

    std::string caller = caller_;
    std::string hostCaller = hostCaller_;
    std::string portBasecaller = toml::get<std::string>(portBasecaller_);

    std::string target_holder, deplete_holder;

    std::stringstream s_stream_1(deplete_files_);

    while (s_stream_1.good()) {

        std::string substr1;
        getline(s_stream_1, substr1, ','); 

        if (substr1.length() > 1)
        {
            if (configReader::filterException(substr1)) {

                deplete = substr1 + ",";
                deplete_holder += deplete;
            }

            else if (!configReader::filterException(substr1)) {

                std::cout << "The deplete file is an fasta file, start building ibf ......." << '\n';

                deplete = output_fileTOML + substr1 + "_.ibf";
                deplete_holder = deplete_holder + deplete + ",";

                toml::value  kmerS = toml::get<toml::table >(IBF).at("kmer_size");
                toml::value  fragment_size = toml::get<toml::table >(IBF).at("fragment_size");
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
 
                target = substr + ",";
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


    configReader::live_target_parser_ targetStruct = { MinKNOW_host, device, deplete_holder, target_holder, output_dir, hostCaller, portBasecaller, caller, MinKNOW_port, basecallThreads, classifyThreads,
                    significance, error_rate, false, false, false };

    

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
