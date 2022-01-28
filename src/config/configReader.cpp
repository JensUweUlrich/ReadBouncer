

#include "configReader.hpp"
#include <string.h>
#include <filesystem>
#include <sstream>
#include <string>
// seqan libraries
#include <seqan/binning_directory.h>
#include <chrono>
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

/**
 * ConfigReader constructor 
 * @param tomlFile /path/to/config.toml 
 */

ConfigReader::ConfigReader(std::string const tomlFile) {

        this->tomlInputFile = tomlFile;
        std::ifstream tomlFileReadBouncer(tomlInputFile, std::ios_base::binary);
        this->configuration_ = toml::parse(tomlFileReadBouncer, /*optional -> */ tomlInputFile);

        if (!tomlFileReadBouncer.is_open()) {
            std::cerr << "Error parsing the toml file: " << tomlInputFile << '\n';
    }
        
};

/**
 * Parse output dir with usage  
 * 
 */

void ConfigReader::parse_general(){

	try
	{
		 this->log_dir = toml::find<std::string>(this->configuration_, "log_directory");
		 this->log_dir = log_dir.make_preferred();

		 if (!std::filesystem::is_directory(this->log_dir) || !std::filesystem::exists(this->log_dir))
		{
			std::filesystem::create_directories(this->log_dir);
		}

		 this->output_dir = toml::find<std::string>(this->configuration_, "output_directory");
		 this->output_dir = this->output_dir.make_preferred();

		 if (!std::filesystem::is_directory(this->output_dir) || !std::filesystem::exists(this->output_dir))
		{
			 std::filesystem::create_directories(this->output_dir); 
		}
         this->usage = toml::find<std::string>(configuration_, "usage");

	}
	catch (const toml::exception& e)
	{
		std::cerr << "Could not parse " << tomlInputFile << std::endl;
		std::cerr << e.what() << std::endl;
	}
	catch (std::out_of_range& e)
	{
		std::cerr << "Error in " << tomlInputFile << std::endl;
		std::cerr << e.what() << std::endl;

	}
}


/**
 * Write log file based on usage path/to/log_dir/configLog.toml
 * @param usage [build, classify, target, test]
 */

void ConfigReader::createLog(std::string& usage){

    std::filesystem::path configLog(this->output_dir); 
    configLog /= "configLog.toml";
    std::fstream outputLog(configLog, std::ios::app | std::ios::out | std::ios::in);


    toml::value tbl;
    toml::value target_files(toml::array{});
		for (std::filesystem::path file : IBF_Parsed.target_files)
			target_files.push_back(file.string());

    toml::value deplete_files(toml::array{});
		for (std::filesystem::path file : IBF_Parsed.deplete_files)
			deplete_files.push_back(file.string());

    toml::value read_files(toml::array{});
		for (std::filesystem::path file : IBF_Parsed.read_files)
			read_files.push_back(file.string());
    

    if (usage == "build"){

         tbl = toml::value{ {
		{usage, toml::table{{
				{ "target_files", target_files},
				{ "deplete_files", deplete_files},
				{ "kmer-size", IBF_Parsed.size_k },
				{ "threads", IBF_Parsed.threads },
				{ "fragment-size", IBF_Parsed.fragment_size}
				 }}
				},
		} };
        
    }

    else if (usage == "classify"){

        tbl = toml::value{ {
		{usage, toml::table{{
				{ "target_files", target_files},
				{ "deplete_files", deplete_files},
                { "read_files", read_files},
				{ "kmer-size", IBF_Parsed.size_k },
				{ "threads", IBF_Parsed.threads },
				{ "fragment-size", IBF_Parsed.fragment_size},
                { "exp_seq_error_rate", IBF_Parsed.error_rate},
                { "chunk_length", IBF_Parsed.chunk_length},
                { "max_chunks", IBF_Parsed.max_chunks},
                
				 }}
				},
		} };
     }

    else if (usage == "target"){

         tbl = toml::value{ {
		{usage, toml::table{{
				{ "target_files", target_files},
				{ "deplete_files", deplete_files},
				{ "kmer-size", IBF_Parsed.size_k },
				{ "threads", IBF_Parsed.threads },
				{ "fragment-size", IBF_Parsed.fragment_size},
                { "exp_seq_error_rate", IBF_Parsed.error_rate},
                { "chunk_length", IBF_Parsed.chunk_length},
                { "max_chunks", IBF_Parsed.max_chunks},
                {"host" , MinKNOW_Parsed.host},
                {"port" , MinKNOW_Parsed.port},
                {"flowcell" , MinKNOW_Parsed.flowcell},
                {"caller", Basecaller_Parsed.caller},
                {"host", Basecaller_Parsed.guppy_host},
                {"port", Basecaller_Parsed.guppy_port},
                {"threads", Basecaller_Parsed.basecall_threads},
                
				 }}
				},
		} };
    }
    else if (usage == "test"){

         tbl = toml::value{ {
		{usage, toml::table{{
                {"host" , MinKNOW_Parsed.host},
                {"port" , MinKNOW_Parsed.port},
                {"flowcell" , MinKNOW_Parsed.flowcell},
				 }}
				},
		} };

    }
        

    auto start = std::chrono::system_clock::now();
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	outputLog << "#Computation time: " << std::ctime(&end_time) << '\n';
	outputLog << toml::format(tbl) << '\n';
	outputLog.close();
}


/**
 * Check if the target/deplete input files are IBF or not
 * @param  file input deplete/target file
 * @return Bool due to decision
 * @throw seqan::Exception if fasta file
 */

bool ConfigReader::filterException(std::filesystem::path& file) {

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


/**
 * Parse all parameters from [IBF] module
 * @throw ConfigReader::Exception
 */

void ConfigReader::readIBF(){

    std::vector<std::string> rf_tmp;

    try
    {
        IBF_Parsed.size_k= toml::find_or<int>(this->configuration_, "IBF", "kmer_size", 13);
        IBF_Parsed.fragment_size = toml::find_or<int>(this->configuration_, "IBF", "fragment_size", 100000);
        IBF_Parsed.threads = toml::find_or<int>(this->configuration_, "IBF", "threads", 1);
        IBF_Parsed.error_rate = toml::find_or<double>(this->configuration_, "IBF", "exp_seq_error_rate", 0.1);
        IBF_Parsed.chunk_length = toml::find_or<int>(this->configuration_, "IBF", "chunk_length", 250);
        IBF_Parsed.max_chunks = toml::find_or<int>(this->configuration_, "IBF", "max_chunks", 5);
    }
    catch (std::out_of_range& e)
    {
        // TODO: write message in log file
        throw ConfigReader(e.what());
    }

    try
	{
		std::vector<std::string> tmp = toml::find<std::vector<std::string>>(this->configuration_, "IBF", "target_files");
		for (std::string s : tmp)
			IBF_Parsed.target_files.emplace_back((std::filesystem::path(s)).make_preferred());
		tmp.clear();
		tmp = toml::find<std::vector<std::string>>(this->configuration_, "IBF", "deplete_files");
		for (std::string s : tmp)
			IBF_Parsed.deplete_files.emplace_back((std::filesystem::path(s)).make_preferred()); 
	}
	catch (toml::exception& e)
	{
		throw ConfigReader(e.what());

	}

    for (std::filesystem::path file : IBF_Parsed.target_files)
		{
			if (!std::filesystem::exists(file))
			{
				// TODO: write message in log file
				throw ConfigReader("[Error] The following target file does not exist: " + file.string());
			}
		}

		for (std::filesystem::path file : IBF_Parsed.deplete_files)
		{
			if (!std::filesystem::exists(file))
			{
				// TODO: write message in log file
				throw ConfigReader("[Error] The following deplete file does not exist: " + file.string());
			}
		}


    try
    {
        rf_tmp = toml::find<std::vector<std::string>>(this->configuration_, "IBF", "read_files");
    }
    catch (toml::exception& e)
	{

		throw ConfigReader(e.what());
	}

    for (std::string file : rf_tmp)
    {
        std::filesystem::path rf(file);
        rf = rf.make_preferred();

        if (!std::filesystem::exists(rf))
        {
            // TODO: write message in log file
            throw ConfigReader("[Error] The following read file does not exist: " + rf.string());
        }
        else
        {
            IBF_Parsed.read_files.emplace_back(std::move(rf));
        }

    }

}


/**
 * Parse all parameters from [MinKNOW] module
 * @throw ConfigReader::Exception
 */

void ConfigReader::readMinKNOW(){
    

    try
    {
        toml::value MinKNOW = toml::find(this->configuration_, "MinKNOW");
        MinKNOW_Parsed.flowcell = toml::find<std::string>(MinKNOW, "flowcell");
        MinKNOW_Parsed.host = toml::find_or<std::string>(MinKNOW, "host", "127.0.0.1");
        MinKNOW_Parsed.port = toml::find_or<std::string>(MinKNOW, "port", "9501");
        //channels = toml::find_or<std::vector<int>>(MinKNOW, "channels", std::vector<int>{});
    
    }
    catch (std::out_of_range& e)
    {
        // TODO: write message in log file
        throw ConfigReaderException(e.what());
    }

}

/**
 * Parse all parameters from [Basecaller] module
 * @throw ConfigReader::Exception
 */

void ConfigReader::readBasecaller(){

    try
    {
        toml::value basecaller = toml::find(this->configuration_, "Basecaller");
        Basecaller_Parsed.caller = toml::find_or<std::string>(basecaller, "caller", "DeepNano");
        Basecaller_Parsed.guppy_host = toml::find<std::string>(basecaller, "host");
        Basecaller_Parsed.guppy_port = toml::find_or<std::string>(basecaller, "port", "5555");
        Basecaller_Parsed.basecall_threads = toml::find_or<int>(basecaller, "threads", 3);
        //Basecaller_Parsed.guppy_config = toml::find_or<std::string>(basecaller, "config", "dna_r9.4.1_450bps_fast");
    }
    catch (std::out_of_range& e)
    {
        // TODO: write message in log file
        throw ConfigReaderException(e.what());
    }

}

/**
 * Call private methods to parse parameters from the different three moduls
 * @throw ConfigReader::Exception
 */

void ConfigReader::parse(){
    

    ConfigReader::readIBF();
    ConfigReader::readMinKNOW();
    ConfigReader::readBasecaller();
    
}
    

//@TODO
/**
 * Parse parameters from toml file for live reads targeting
* @param : Toml output file, usage, a list of target and deplete files 
* @return: Struct with the needed parameters for live target ((Targeted Sequencing))
*/

/*


ConfigReader::Target_Params ConfigReader::targetReader( std::string& usage,
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

        if (ConfigReader::filterException(file))
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

        if (ConfigReader::filterException(file))
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
};*/