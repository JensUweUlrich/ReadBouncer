#include "gtest/gtest.h"
#include "configReader.hpp"

//https://chromium.googlesource.com/external/github.com/google/googletest/+/refs/tags/release-1.8.0/googletest/docs/Primer.md#:~:text=Google%20Test%20assertions%20are%20macros,along%20with%20a%20failure%20message.


class ConfigReaderTest: public ::testing::Test
{
	protected:

		ConfigReader* config;
		// prepare the objects for each test.
		void SetUp() override
		{
			config = new ConfigReader("/mnt/c/bug29/readbouncer/config.toml");

            std::filesystem::path path = "/mnt/c/bug29/readbouncer/config.toml";
            tomlF /= path;
            std::ifstream tomlFileReadBouncer(path, std::ios_base::binary);
            configurationSet = toml::parse(tomlFileReadBouncer, /*optional -> */ path.string());
			
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


// TEST_F:  writing two or more tests that operate on similar data
// 
TEST_F(ConfigReaderTest, TestConfigReaderConstructur)
{
    std::cout<<'\n';
    std::cout << "Testing ConfigReaderConstructur......................................" << '\n';
    std::cout<<'\n';
    // check if file exists
    EXPECT_TRUE(std::filesystem::exists(tomlF));

    // fatal failure if both configurations are unequal
    //ASSERT_STREQ only for C strings // ASSERT_EQ any two values! could be two IBF
	EXPECT_EQ(config->configuration_, configurationSet);
    

}


TEST_F(ConfigReaderTest, TestParseGeneral)
{
    config->parse_general();
    //EXPECT_ANY_THROW(config->parse_general(), ConfigReaderException);
    //EXPECT_NO_THROW(config->parse_general());

    log        = toml::find<std::string>(configurationSet, "log_directory");
    output     = toml::find<std::string>(configurationSet, "output_directory");
    subcommand = toml::find<std::string>(configurationSet, "usage");

    std::cout << "Testing parse_general......................................" << '\n';
    std::cout<<'\n';

    // Test if toml::find acutally parses correctly! 
    EXPECT_EQ("RB_out/logs", log);// Nonfatal assertion to check afterwards if ConfigReader throws an exception at this point! 
    EXPECT_EQ("RB_out", output);
    if(subcommand == "build" || subcommand == "test" || subcommand == "classify" || subcommand == "target"){

        SUCCEED();
    }

    else{

        FAIL(); 
    }


    // ConfigReader Test
	EXPECT_EQ(config->log_dir, log);// Nonfatal assertion to check afterwards if ConfigReader throws an exception at this point! 
    EXPECT_EQ(config->output_dir, output);
    EXPECT_EQ(config->usage, subcommand);

    if (Test::HasFailure()){
        
        EXPECT_THROW(config->parse_general(),ConfigReaderException);
    }

    else{

        EXPECT_NO_THROW(config->parse_general());
    } 
}



TEST_F(ConfigReaderTest, TestParse)
{
    config->parse();

   std::vector<std::filesystem::path> target_files{};
   std::vector<std::filesystem::path> deplete_files{};
   std::vector<std::filesystem::path> read_files{};


    // readIBF()

   std::cout << "Testing parse......................................" << '\n';
   std::cout<<'\n';
   
   //[IBF]
   int    k  = toml::find_or<int>(configurationSet, "IBF", "kmer_size", 13);
   int    f  = toml::find_or<int>(configurationSet, "IBF", "fragment_size", 100000);
   int    t  = toml::find_or<int>(configurationSet, "IBF", "threads", 1);
   double e  = toml::find_or<double>(configurationSet, "IBF", "exp_seq_error_rate", 0.1);
   int    cl = toml::find_or<int>(configurationSet, "IBF", "chunk_length", 250);
   int    mc = toml::find_or<int>(configurationSet, "IBF", "max_chunks", 5);

   std::vector<std::string> target_tmp = toml::find<std::vector<std::string>>(configurationSet, "IBF", "target_files");
   for (std::string s : target_tmp) target_files.emplace_back((std::filesystem::path(s)).make_preferred());
  

   std::vector<std::string> deplet_tmp = toml::find<std::vector<std::string>>(configurationSet, "IBF", "deplete_files");
   for (std::string s : deplet_tmp) deplete_files.emplace_back((std::filesystem::path(s)).make_preferred());
   

   std::vector<std::string>  reads_tmp = toml::find<std::vector<std::string>>(configurationSet, "IBF", "read_files");
   for (std::string s : reads_tmp) read_files.emplace_back((std::filesystem::path(s)).make_preferred());

   // Test if toml::find acutally parses correctly! 
   EXPECT_EQ(k, 15); 
   EXPECT_EQ(f, 100000); 
   EXPECT_EQ(t, 3); 
   EXPECT_EQ(e, 0.1); 
   EXPECT_EQ(cl, 350); 
   EXPECT_EQ(mc, 1); 

   std::vector<std::filesystem::path> target  = {"/mnt/c/ReadBouncerToml/build/main/Release/Listeria_monocytogenes_ATCC_19115_.fasta", "/mnt/c/ReadBouncerToml/build/main/Release/Pseudomonas_aeruginosa_complete_genome.fasta"};
   std::vector<std::filesystem::path> deplete = {"/mnt/c/ReadBouncerToml/build/main/Release/Bacillus_subtilis_complete_genome.ibf", "/mnt/c/ReadBouncerToml/build/main/Release/Enterococcus_faecalis_complete_genome.fasta"};
   std::vector<std::filesystem::path> reads   = {"/mnt/c/ReadBouncerToml/build/main/Release/Listeria.fastq","/mnt/c/ReadBouncerToml/build/main/Release/SaccharomycesReal.fasta"};

   EXPECT_EQ(target_files, target);
   EXPECT_EQ(deplete_files, deplete);
   EXPECT_EQ(read_files, reads);

   

   // test all
   EXPECT_EQ(k, config->IBF_Parsed.size_k); 
   EXPECT_EQ(f, config->IBF_Parsed.fragment_size); 
   EXPECT_EQ(t, config->IBF_Parsed.threads); 
   EXPECT_EQ(e, config->IBF_Parsed.error_rate); 
   EXPECT_EQ(cl, config->IBF_Parsed.chunk_length); 
   EXPECT_EQ(mc, config->IBF_Parsed.max_chunks); 

   EXPECT_EQ(target_files, config->IBF_Parsed.target_files);
   EXPECT_EQ(deplete_files, config->IBF_Parsed.deplete_files);
   EXPECT_EQ(read_files, config->IBF_Parsed.read_files);

   EXPECT_EQ("localhost", config->MinKNOW_Parsed.host);
   EXPECT_EQ("9501", config->MinKNOW_Parsed.port);
   EXPECT_EQ("MS00000", config->MinKNOW_Parsed.flowcell);

   EXPECT_EQ("DeepNano", config->Basecaller_Parsed.caller);
   EXPECT_EQ("127.0.0.1", config->Basecaller_Parsed.guppy_host);
   EXPECT_EQ("9501", config->Basecaller_Parsed.guppy_port);
   EXPECT_EQ(3, config->Basecaller_Parsed.basecall_threads);

   if (Test::HasFailure()){
        
        EXPECT_THROW(config->parse_general(),ConfigReaderException);
    }

    else{

        EXPECT_NO_THROW(config->parse_general());
    } 


}

/*
TEST(ConfigReaderTest, TestIBFStruct)// what ConfigReader sends to main! 
{
    config->parse();

    struct IBF_Struct_Test
    {
        int size_k = 15;
        int fragment_size = 100000;
        int threads = 3;
        std::vector<std::filesystem::path> target_files = {"/mnt/c/ReadBouncerToml/build/main/Release/Listeria_monocytogenes_ATCC_19115_.fasta", "/mnt/c/ReadBouncerToml/build/main/Release/Pseudomonas_aeruginosa_complete_genome.fasta"};
	    std::vector<std::filesystem::path> deplete_files = {"/mnt/c/ReadBouncerToml/build/main/Release/Bacillus_subtilis_complete_genome.ibf", "/mnt/c/ReadBouncerToml/build/main/Release/Enterococcus_faecalis_complete_genome.fasta"};
        std::vector<std::filesystem::path> reads = {"/mnt/c/ReadBouncerToml/build/main/Release/Listeria.fastq","/mnt/c/ReadBouncerToml/build/main/Release/SaccharomycesReal.fasta"};
        double error_rate = 0.1;
        int chunk_length = 360;
        int max_chunks = 1;
    };

    EXPECT_EQ(IBF_Struct_Test, config->IBF_Parsed); 

    

}

TEST(ConfigReaderTest, TestMinKNOWStruct)// what ConfigReader sends to main! 
{
    struct MinKNOW_Struct_Test
    {
        std::string host = "localhost";
        std::string port = "9501";
        std::string flowcell = "MS00000";
        uint16_t minChannel = 1;
        uint16_t maxChannel = 512;
    };

    //ASSERT_EQ(MinKNOW_Struct_Test, config->MinKNOW_Parsed); 

}

TEST(ConfigReaderTest, TestBasecallerStruct)// what ConfigReader sends to main! 
{
  struct Basecaller_Struct_Test
    {
        std::string caller = "DeepNano";
        std::string guppy_host = "127.0.0.1";
        std::string guppy_port = "9501";
        int basecall_threads = 3;
        std::string guppy_config = "dna_r9.4.1_450bps_fast";
    };

    //EXPECT_EQ(Basecaller_Struct_Test, config->Basecaller_Parsed); 

}*/

int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}