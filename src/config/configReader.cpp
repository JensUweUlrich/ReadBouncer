

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
        try
        {
            this->configuration_ = toml::parse(tomlFileReadBouncer, /*optional -> */ tomlInputFile);
        }
        catch (toml::exception& e)
        {
            throw ConfigReaderException(e.what());
        }

        // TODO: throw ConfigReaderException
        if (!tomlFileReadBouncer.is_open()) {
            std::cerr << "Error parsing the toml file: " << tomlInputFile << '\n';
    }
        
};

/**
 * Parse main parameters from toml file [log_directory, output_directory and usage] 
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
        throw ConfigReaderException(e.what());
	}
	catch (std::out_of_range& e)
	{
        throw ConfigReaderException(e.what());

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
                {"host" , MinKNOW_Parsed.host},
                {"port" , MinKNOW_Parsed.port},
                {"flowcell" , MinKNOW_Parsed.flowcell},
                {"MinChannel", MinKNOW_Parsed.minChannel},
                {"MaxChannel", MinKNOW_Parsed.maxChannel},
                {"caller", Basecaller_Parsed.caller},
                {"GuppyConfig", Basecaller_Parsed.guppy_config},
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
 * Check if the target/deplete input files are IBF or fasta files
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
        throw ConfigReaderException(e.what());
    }

    try
	{
		std::vector<std::string> tmp = toml::find<std::vector<std::string>>(this->configuration_, "IBF", "target_files");
		for (std::string s : tmp)
			IBF_Parsed.target_files.emplace_back((std::filesystem::path(s)).make_preferred());
	}
	catch (toml::exception& e)
	{
		throw ConfigReaderException(e.what());
	}
    catch (std::out_of_range& e)
    {
        // Do nothing
        // sometimes we only want to specify deplete files
    }


    try
    {
        std::vector<std::string> tmp = toml::find<std::vector<std::string>>(this->configuration_, "IBF", "deplete_files");
        for (std::string s : tmp)
            IBF_Parsed.deplete_files.emplace_back((std::filesystem::path(s)).make_preferred());
    }
    catch (toml::exception& e)
    {
        throw ConfigReaderException(e.what());
    }
    catch (std::out_of_range& e)
    {

        // Do nothing
        // sometimes we only want to specify target files
    }

    if (!(this->usage.compare("test") == 0))
    {
        if (IBF_Parsed.deplete_files.size() + IBF_Parsed.target_files.size() == 0)
        {
            throw ConfigReaderException("[Error] At least one target or deplete file has to be specified!");
        }
    }

    for (std::filesystem::path file : IBF_Parsed.target_files)
	{
		if (!std::filesystem::exists(file))
		{
			// TODO: write message in log file
			throw ConfigReaderException("[Error] The following target file does not exist: " + file.string());
		}
	}

	for (std::filesystem::path file : IBF_Parsed.deplete_files)
	{
		if (!std::filesystem::exists(file))
		{
			// TODO: write message in log file
			throw ConfigReaderException("[Error] The following deplete file does not exist: " + file.string());
		}
	}

    try
    {
        rf_tmp = toml::find<std::vector<std::string>>(this->configuration_, "IBF", "read_files");
    }
    catch (toml::exception& e)
	{
		throw ConfigReaderException(e.what());
	}
    catch (std::out_of_range& e)
    {
        if (this->usage.compare("classify") == 0)
        {
            throw ConfigReaderException(e.what());
        }
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
            IBF_Parsed.read_files.emplace_back(std::move(rf));
        }

    }

}


/**
 * Parse all parameters from [MinKNOW] module
 * @throw ConfigReader::Exception
 */

void ConfigReader::readMinKNOW(){
    
    toml::value MinKNOW;
    try
    {
        MinKNOW = toml::find(this->configuration_, "MinKNOW");
    }
    catch (std::out_of_range& e)
    {
        // Do nothing and use default values
        return;
    }

    try
    {
        
        MinKNOW_Parsed.flowcell = toml::find<std::string>(MinKNOW, "flowcell");
        MinKNOW_Parsed.host = toml::find_or<std::string>(MinKNOW, "host", "127.0.0.1");
        MinKNOW_Parsed.port = toml::find_or<std::string>(MinKNOW, "port", "9501");
        std::vector<int> channels = toml::find_or<std::vector<int>>(MinKNOW, "channels", std::vector<int>{});
        MinKNOW_Parsed.token_path = toml::find_or<std::string>(MinKNOW, "token_path", "");
        if (channels.size() == 2)
        {
            MinKNOW_Parsed.minChannel = (uint16_t) channels[0];
            MinKNOW_Parsed.maxChannel = (uint16_t) channels[1];
        }
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

    toml::value basecaller;
    try
    {
        basecaller = toml::find(this->configuration_, "Basecaller");
    }
    catch (std::out_of_range& e)
    {
        // Do nothing and use default values
        return;
    }


    try
    {
        Basecaller_Parsed.caller = toml::find_or<std::string>(basecaller, "caller", "DeepNano");
        Basecaller_Parsed.guppy_host = toml::find_or<std::string>(basecaller, "host", "127.0.0.1");
        Basecaller_Parsed.guppy_port = toml::find_or<std::string>(basecaller, "port", "5555");
        Basecaller_Parsed.basecall_threads = toml::find_or<int>(basecaller, "threads", 3);
        Basecaller_Parsed.guppy_config = toml::find_or<std::string>(basecaller, "config", "dna_r9.4.1_450bps_fast");
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
    
    try
    {
        ConfigReader::readIBF();
        ConfigReader::readMinKNOW();
        ConfigReader::readBasecaller();
    }
    catch (ConfigReaderException& e)
    {
        throw;
    }
}
    
