#include "gtest/gtest.h"
//#include "gmock/gmock.h"
#include "configReader.hpp"
#include <algorithm>



class ConfigReaderTest: public ::testing::Test
{
	protected:

		ConfigReader* config;
        
		// prepare the objects for each test.
		void SetUp() override
		{
            tomlF = "./test_toml_file.toml";
			config = new ConfigReader(tomlF.string());

            std::ifstream tomlFileReadBouncer(tomlF, std::ios_base::binary);
            configurationSet = toml::parse(tomlFileReadBouncer, /*optional -> */ tomlF.string());
			
		}
        // release any resources we allocated in SetUp()
		void TearDown() override
		{
			delete config;
		}

    public:
     toml::basic_value<struct toml::discard_comments, std::unordered_map, std::vector> configurationSet;
     std::filesystem::path tomlF, log, output; 
     std::string subcommand;

};

/*
* Test reading config.toml and construction the configuration settings
*/

TEST_F(ConfigReaderTest, ConfigReaderConstructurTest)
{

    EXPECT_TRUE(std::filesystem::exists(tomlF));
	EXPECT_EQ(config->configuration_, configurationSet); // check if the loaded configurations are same

}

/*
* Test parsing [General] from test_toml_file.toml 
*/

TEST_F(ConfigReaderTest, ParseGeneralTest)
{
    config->parse_general(); // parse using the main ConfigReader class 

    log        = toml::find<std::string>(configurationSet, "log_directory"); // testing class
    output     = toml::find<std::string>(configurationSet, "output_directory"); // testing class
    subcommand = toml::find<std::string>(configurationSet, "usage"); // testing class


    // Test if toml::find acutally parses correctly! 
    EXPECT_EQ("log_tests", log);// Non-fatal assertion to check afterwards if ConfigReader throws an exception at this point! 
    EXPECT_EQ("output_tests", output);

    // Breaking point of testing file
    if(subcommand == "build" || subcommand == "test" || subcommand == "classify" || subcommand == "target"){

        SUCCEED();
    }

    else{

        std::cerr << "[ERROR] This test failed due to the changes in the used commands [build, target, classify, test]" << '\n';
        FAIL(); 
    }

    // Testing configReader parsing (general) process
	EXPECT_EQ(config->log_dir, log);
    EXPECT_EQ(config->output_dir, output);
    EXPECT_EQ(config->usage, subcommand);

    // If any test fails! Expect an throw exception from configReader lib as ConfigReaderException
    if (Test::HasFailure()){
        
        EXPECT_THROW(config->parse_general(),ConfigReaderException);// Does ConfigReader throws exception? 
    }

    else{

        EXPECT_NO_THROW(config->parse_general());
    } 
}



/*
* Testparse function from ConfigReader class
*/

TEST_F(ConfigReaderTest, ParseTest)
{
    config->parse();

   EXPECT_EQ(15, config->IBF_Parsed.size_k); 
   EXPECT_EQ(100000, config->IBF_Parsed.fragment_size); 
   EXPECT_EQ(3, config->IBF_Parsed.threads); 
   EXPECT_EQ(0.1, config->IBF_Parsed.error_rate); 
   EXPECT_EQ(360, config->IBF_Parsed.chunk_length); 
   EXPECT_EQ(1, config->IBF_Parsed.max_chunks); 

   std::vector<std::filesystem::path> target_files = { "./target_test1.fasta", "./target_test2.fasta" };
   std::vector<std::filesystem::path> deplete_files = { "./deplete_test1.fasta", "./deplete_test2.fasta" };
   std::vector<std::filesystem::path> reads = {"./reads_test.fastq"};

   EXPECT_EQ(target_files, config->IBF_Parsed.target_files);
   EXPECT_EQ(deplete_files, config->IBF_Parsed.deplete_files);
   EXPECT_EQ(reads, config->IBF_Parsed.read_files);

   EXPECT_EQ("localhost", config->MinKNOW_Parsed.host);
   EXPECT_EQ("9502", config->MinKNOW_Parsed.port);
   EXPECT_EQ("MS00000", config->MinKNOW_Parsed.flowcell);
   EXPECT_EQ("test/tmp/minknow-auth-token.json", config->MinKNOW_Parsed.token_path);
   EXPECT_EQ(1, config->MinKNOW_Parsed.minChannel);
   EXPECT_EQ(512, config->MinKNOW_Parsed.maxChannel);


   EXPECT_EQ("DeepNano", config->Basecaller_Parsed.caller);
   EXPECT_EQ("127.0.0.1", config->Basecaller_Parsed.guppy_host);
   EXPECT_EQ("9502", config->Basecaller_Parsed.guppy_port);
   EXPECT_EQ(3, config->Basecaller_Parsed.basecall_threads);
   EXPECT_EQ("dna_r9.4.1_450bps_fast", config->Basecaller_Parsed.guppy_config);

   if (Test::HasFailure()){
        
        EXPECT_THROW(config->parse_general(),ConfigReaderException);
    }

    else{

        EXPECT_NO_THROW(config->parse_general());
    } 


}

/*
* Testing last struct equality [IBF]
*/

TEST_F(ConfigReaderTest, IBFStructTest){

    config->parse();

    struct IBFTestStruct
    {
        int size_k = 15;
        int fragment_size = 100000;
        int threads = 3;
        std::vector<std::filesystem::path> target_files = { "./target_test1.fasta", "./target_test2.fasta" };
        std::vector<std::filesystem::path> deplete_files = { "./deplete_test1.fasta", "./deplete_test2.fasta" };
        std::vector<std::filesystem::path> reads = { "./reads_test.fastq" };
        double error_rate = 0.1;
        int chunk_length = 360;
        int max_chunks = 1;
    }IBFTestStruct_;

    
    EXPECT_EQ(IBFTestStruct_.size_k, config->IBF_Parsed.size_k);
    EXPECT_EQ(IBFTestStruct_.threads, config->IBF_Parsed.threads);
    EXPECT_EQ(IBFTestStruct_.fragment_size, config->IBF_Parsed.fragment_size);
    EXPECT_EQ(IBFTestStruct_.target_files, config->IBF_Parsed.target_files);
    EXPECT_EQ(IBFTestStruct_.deplete_files, config->IBF_Parsed.deplete_files);
    EXPECT_EQ(IBFTestStruct_.reads, config->IBF_Parsed.read_files);
    EXPECT_EQ(IBFTestStruct_.error_rate, config->IBF_Parsed.error_rate);
    EXPECT_EQ(IBFTestStruct_.chunk_length, config->IBF_Parsed.chunk_length);
    EXPECT_EQ(IBFTestStruct_.max_chunks, config->IBF_Parsed.max_chunks);
}

/*
* Struct equality [MinKNOW]
*/

TEST_F(ConfigReaderTest, MinKNOWStructTest){

    config->parse();

    struct MinknowTestStruct
    {
        std::string host = "localhost";
        std::string port = "9502";
        std::string flowcell = "MS00000";
        std::filesystem::path token_path_s{"test/tmp/minknow-auth-token.json"};
        uint16_t minChannel = 1;
        uint16_t maxChannel = 512;

    }MinknowTestStruct_;

   EXPECT_EQ(MinknowTestStruct_.host, config->MinKNOW_Parsed.host);
   EXPECT_EQ(MinknowTestStruct_.port, config->MinKNOW_Parsed.port);
   EXPECT_EQ(MinknowTestStruct_.flowcell, config->MinKNOW_Parsed.flowcell);
   EXPECT_EQ(MinknowTestStruct_.token_path_s, config->MinKNOW_Parsed.token_path);
   EXPECT_EQ(MinknowTestStruct_.minChannel, config->MinKNOW_Parsed.minChannel);
   EXPECT_EQ(MinknowTestStruct_.maxChannel, config->MinKNOW_Parsed.maxChannel);

}

/*
* Struct equality [IBaseCaller]
*/


TEST_F(ConfigReaderTest, BaseCallerStructTest){

    config->parse();

    struct BasecallerTestStruct
    {
        std::string caller = "DeepNano";
        std::string guppy_host = "127.0.0.1";
        std::string guppy_port = "9502";
        int basecall_threads = 3;
        std::string guppy_config = "dna_r9.4.1_450bps_fast";

    }BasecallerTestStruct_;

   EXPECT_EQ(BasecallerTestStruct_.caller, config->Basecaller_Parsed.caller);
   EXPECT_EQ(BasecallerTestStruct_.guppy_host, config->Basecaller_Parsed.guppy_host);
   EXPECT_EQ(BasecallerTestStruct_.guppy_port, config->Basecaller_Parsed.guppy_port);
   EXPECT_EQ(BasecallerTestStruct_.basecall_threads, config->Basecaller_Parsed.basecall_threads);
   EXPECT_EQ(BasecallerTestStruct_.guppy_config, config->Basecaller_Parsed.guppy_config);

}

/*
* Test writing the correct log file content while using the "build" command
* TODO: Test table length! (if we add any new items inside the main function!) 
*/
TEST_F(ConfigReaderTest, CreateLogIBFTest){

    config->parse_general();
    config->parse();

    // IBF log
    std::string usage{ "build" };
    config->createLog(usage);

    std::filesystem::path logFile = config->output_dir += "/configLog.toml";
    toml::basic_value<struct toml::discard_comments, std::unordered_map, std::vector> configurationSetting;
    std::ifstream tomlLogFile(logFile, std::ios_base::binary);
    configurationSetting = toml::parse(tomlLogFile, /*optional -> */ logFile.string());

    if (std::filesystem::exists(logFile)) {

        SUCCEED();
    }
    else {

        FAIL();
    }

    int kmer_size = toml::find<int>(configurationSetting, "build", "kmer_size");
    int fragment_size = toml::find<int>(configurationSetting, "build", "fragment_size");
    int threads =  toml::find<int>(configurationSetting, "build", "threads");

    std::vector<std::filesystem::path> target_files, deplete_files;
    std::vector<std::string> target_files_tmp = toml::find<std::vector<std::string>>(configurationSetting, "build", "target_files");
    std::vector<std::string> deplete_files_tmp = toml::find<std::vector<std::string>>(configurationSetting, "build", "deplete_files");

    for (std::string s : target_files_tmp) target_files.emplace_back((std::filesystem::path(s)).make_preferred());
    for (std::string s : deplete_files_tmp) deplete_files.emplace_back((std::filesystem::path(s)).make_preferred());
    
    std::vector<std::filesystem::path> target_files_test = { "./target_test1.fasta", "./target_test2.fasta" };
    std::vector<std::filesystem::path> deplete_files_test = { "./deplete_test1.fasta", "./deplete_test2.fasta" };
    
    EXPECT_EQ(kmer_size, 15);
    EXPECT_EQ(fragment_size, 100000);
    EXPECT_EQ(threads, 3);
    EXPECT_EQ(target_files, target_files_test);
    EXPECT_EQ(deplete_files, deplete_files_test);

}


/*
* Test writing the correct log file content while using the "classify" command
* TODO: Test table length! (if we add any new items inside the main function!) 
*/
TEST_F(ConfigReaderTest, CreateLogClassifyTest) {

    config->parse_general();
    config->parse();

    std::string usage{ "classify" };
    config->createLog(usage);

    std::filesystem::path logFile = config->output_dir += "/configLog.toml";
    toml::basic_value<struct toml::discard_comments, std::unordered_map, std::vector> configurationSetting;
    std::ifstream tomlLogFile(logFile, std::ios_base::binary);
    configurationSetting = toml::parse(tomlLogFile, /*optional -> */ logFile.string());

    if (std::filesystem::exists(logFile)) {

        SUCCEED();
    }
    else {

        FAIL();
    }

    int kmer_size = toml::find<int>(configurationSetting, "classify", "kmer_size");
    int fragment_size = toml::find<int>(configurationSetting, "classify", "fragment_size");
    int threads = toml::find<int>(configurationSetting, "classify", "threads");
    double exp_seq_error_rate = toml::find<double>(configurationSetting, "classify", "exp_seq_error_rate");
    int chunk_length = toml::find<int>(configurationSetting, "classify", "chunk_length");
    int max_chunks = toml::find<int>(configurationSetting, "classify", "max_chunks");

    std::vector<std::filesystem::path> target_files, deplete_files, read_files;
    std::vector<std::string> target_files_tmp = toml::find<std::vector<std::string>>(configurationSetting, "classify", "target_files");
    std::vector<std::string> deplete_files_tmp = toml::find<std::vector<std::string>>(configurationSetting, "classify", "deplete_files");
    std::vector<std::string> read_files_tmp = toml::find<std::vector<std::string>>(configurationSetting, "classify", "read_files");

    for (std::string s : target_files_tmp) target_files.emplace_back((std::filesystem::path(s)).make_preferred());
    for (std::string s : deplete_files_tmp) deplete_files.emplace_back((std::filesystem::path(s)).make_preferred());
    for (std::string s : read_files_tmp) read_files.emplace_back((std::filesystem::path(s)).make_preferred());

    std::vector<std::filesystem::path> target_files_test = { "./target_test1.fasta", "./target_test2.fasta" };
    std::vector<std::filesystem::path> deplete_files_test = { "./deplete_test1.fasta", "./deplete_test2.fasta" };
    std::vector<std::filesystem::path> read_files_test = { "./reads_test.fastq" };

    EXPECT_EQ(kmer_size, 15);
    EXPECT_EQ(fragment_size, 100000);
    EXPECT_EQ(threads, 3);
    EXPECT_EQ(exp_seq_error_rate, 0.1);// not fixed value! 
    EXPECT_EQ(chunk_length, 360);
    EXPECT_EQ(max_chunks, 1);
    EXPECT_EQ(target_files, target_files_test);
    EXPECT_EQ(deplete_files, deplete_files_test);
    EXPECT_EQ(read_files, read_files_test);

}


/*
* Test writing the correct log file content while using the "target" command
* TODO: Test table length! (if we add any new items inside the main function!) 
*/
TEST_F(ConfigReaderTest, CreateLogTargetTest) {

    config->parse_general();
    config->parse();

    std::string usage{ "target" };
    config->createLog(usage);

    std::filesystem::path logFile = config->output_dir += "/configLog.toml";
    toml::basic_value<struct toml::discard_comments, std::unordered_map, std::vector> configurationSetting;
    std::ifstream tomlLogFile(logFile, std::ios_base::binary);
    configurationSetting = toml::parse(tomlLogFile, /*optional -> */ logFile.string());

    if (std::filesystem::exists(logFile)) {

        SUCCEED();
    }
    else {

        FAIL();
    }

    int kmer_size = toml::find<int>(configurationSetting, "target", "kmer_size");
    int fragment_size = toml::find<int>(configurationSetting, "target", "fragment_size");
    int threads = toml::find<int>(configurationSetting, "target", "threads");
    int basecaller_threads = toml::find<int>(configurationSetting, "target", "Basecaller_threads");
    double exp_seq_error_rate = toml::find<double>(configurationSetting, "target", "exp_seq_error_rate");

    uint16_t MinChannel = toml::find<uint16_t>(configurationSetting, "target", "MinChannel");
    uint16_t MaxChannel = toml::find<uint16_t>(configurationSetting, "target", "MaxChannel");

    std::string host = toml::find<std::string>(configurationSetting, "target", "host");
    std::string port = toml::find<std::string>(configurationSetting, "target", "port");
    std::string flowcell = toml::find<std::string>(configurationSetting, "target", "flowcell");
    std::string token_path = toml::find<std::string>(configurationSetting, "target", "token_path");

    std::string caller = toml::find<std::string>(configurationSetting, "target", "caller");
    std::string caller_host = toml::find<std::string>(configurationSetting, "target", "caller_host");
    std::string caller_port = toml::find<std::string>(configurationSetting, "target", "caller_port");
    std::string config = toml::find<std::string>(configurationSetting, "target", "GuppyConfig");

    std::vector<std::filesystem::path> target_files, deplete_files, read_files;
    std::vector<std::string> target_files_tmp = toml::find<std::vector<std::string>>(configurationSetting, "target", "target_files");
    std::vector<std::string> deplete_files_tmp = toml::find<std::vector<std::string>>(configurationSetting, "target", "deplete_files");

    for (std::string s : target_files_tmp) target_files.emplace_back((std::filesystem::path(s)).make_preferred());
    for (std::string s : deplete_files_tmp) deplete_files.emplace_back((std::filesystem::path(s)).make_preferred());

    std::vector<std::filesystem::path> target_files_test = { "./target_test1.fasta", "./target_test2.fasta" };
    std::vector<std::filesystem::path> deplete_files_test = { "./deplete_test1.fasta", "./deplete_test2.fasta" };

    
    EXPECT_EQ(kmer_size, 15);
    EXPECT_EQ(fragment_size, 100000);
    EXPECT_EQ(threads, 3);
    EXPECT_EQ(basecaller_threads, 3);
    EXPECT_EQ(exp_seq_error_rate, 0.1);
    EXPECT_EQ(MinChannel, 1);
    EXPECT_EQ(MaxChannel, 512);
    EXPECT_EQ(host, "localhost");
    EXPECT_EQ(port, "9502");
    EXPECT_EQ(flowcell, "MS00000");
    EXPECT_EQ(token_path, "test/tmp/minknow-auth-token.json");
    EXPECT_EQ(caller, "DeepNano");
    EXPECT_EQ(caller_host, "127.0.0.1");
    EXPECT_EQ(caller_port, "9502");
    EXPECT_EQ(config, "dna_r9.4.1_450bps_fast");
    EXPECT_EQ(basecaller_threads, 3);
    EXPECT_EQ(target_files, target_files_test);
    EXPECT_EQ(deplete_files, deplete_files_test);



}


/*
* Test writing the correct log file content while using the "test" command
* TODO: Test table length! (if we add any new items inside the main function!) 
*/
TEST_F(ConfigReaderTest, CreateLogConnectionTestTest) {

    config->parse_general();
    config->parse();

    std::string usage{ "test" };
    config->createLog(usage);

    std::filesystem::path logFile = config->output_dir += "/configLog.toml";
    toml::basic_value<struct toml::discard_comments, std::unordered_map, std::vector> configurationSetting;
    std::ifstream tomlLogFile(logFile, std::ios_base::binary);
    configurationSetting = toml::parse(tomlLogFile, /*optional -> */ logFile.string());

    //std::cout << "Size of setting: " << configurationSetting.size() << std::endl;
    if (std::filesystem::exists(logFile)) {

        SUCCEED();
    }
    else {

        FAIL();
    }

    std::string host = toml::find<std::string>(configurationSetting, "test", "host");
    std::string port = toml::find<std::string>(configurationSetting, "test", "port");
    std::string flowcell = toml::find<std::string>(configurationSetting, "test", "flowcell");
    std::string token_path = toml::find<std::string>(configurationSetting, "test", "token_path");

    EXPECT_EQ(host, "localhost");
    EXPECT_EQ(port, "9502");
    EXPECT_EQ(flowcell, "MS00000");
    EXPECT_EQ(token_path, "test/tmp/minknow-auth-token.json");

    // Creating new empty file in order to avoid having the same table twice 
    std::ofstream file(logFile.string());
}


int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}