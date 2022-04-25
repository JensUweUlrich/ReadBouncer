




class ClassifyTest: public ::testing::Test
{
	//MockIBF mock_ibf_;
	
	protected:
        
		// prepare the objects for each test.
		void SetUp() override
		{
			//ibf = new interleave::IBF();
		}
        // release any resources we allocated in SetUp()
		void TearDown() override
		{
			//delete ibf;
		}

    public:// mocking
	  
       int testCounter
	   {0};
	    

};


typedef std::vector<interleave::IBFMeta> IBFMetaVector;
 
IBFMetaVector getIBF_(ConfigReader config, bool targetFilter, bool depleteFilter){

	IBFMetaVector DepletionFilters{};
	IBFMetaVector TargetFilters{};

	std::shared_ptr<spdlog::logger> nanolive_logger;

	if(depleteFilter){
		// parse depletion IBF if given as parameter
		for (std::filesystem::path deplete_file : config.IBF_Parsed.deplete_files)
		{
			interleave::IBFMeta filter{};
			filter.name = deplete_file.stem().string();
			interleave::IBF tf{};
			interleave::IBFConfig DepleteIBFconfig{};

			if (config.filterException(deplete_file)){
				try
				{
					DepleteIBFconfig.input_filter_file = deplete_file.string();
					interleave::FilterStats stats = tf.load_filter(DepleteIBFconfig);
					filter.filter = std::move(tf.getFilter());
					interleave::print_load_stats(stats);
				}
				catch (interleave::ParseIBFFileException& e)
				{
					nanolive_logger->error("Error parsing depletion IBF using the following parameters");
					nanolive_logger->error("Depletion IBF file                : " + deplete_file.string());
					nanolive_logger->error("Error message : " + std::string(e.what()));
					nanolive_logger->flush();
					throw;
				}

				DepletionFilters.emplace_back(std::move(filter));
			}
		
		    else
			{
				try
				{
					//ibf_build_parser params;
					std::filesystem::path out = std::filesystem::path(config.output_dir);
					out /= deplete_file.filename();
					out.replace_extension("ibf");
					ibf_build_parser params = { out.string(), deplete_file.string(), false, false, config.IBF_Parsed.size_k, config.IBF_Parsed.threads, config.IBF_Parsed.fragment_size, 0, true };
					filter.filter = buildIBF(params);
					}

				catch (std::out_of_range& e)
				{
					throw ConfigReaderException(e.what());
				}
			DepletionFilters.emplace_back(std::move(filter));
			}
		}
		return DepletionFilters;
	}

	if(targetFilter)
	{
		for (std::filesystem::path target_file : config.IBF_Parsed.target_files)
		{
			interleave::IBFMeta filter{};
			filter.name = target_file.stem().string();
			interleave::IBF tf{};
			interleave::IBFConfig TargetIBFconfig{};
			if (config.filterException(target_file))
			{
				try
				{
					TargetIBFconfig.input_filter_file = target_file.string();
					interleave::FilterStats stats = tf.load_filter(TargetIBFconfig);
					filter.filter = std::move(tf.getFilter());
					interleave::print_load_stats(stats);

				}
				catch (interleave::ParseIBFFileException& e)
				{
					nanolive_logger->error("Error building IBF for target file using the following parameters");
					nanolive_logger->error("Depletion IBF file                : " + target_file.string());
					nanolive_logger->error("Error message : " + std::string(e.what()));
					nanolive_logger->flush();
					throw;
				}

				TargetFilters.emplace_back(std::move(filter));
			}
		
			else
			{
				try
				{
					//ibf_build_parser params;
					std::filesystem::path out = std::filesystem::path(config.output_dir);
					out /= target_file.filename();
					out.replace_extension("ibf");
					ibf_build_parser params = { out.string(), target_file.string(), false, false, config.IBF_Parsed.size_k, config.IBF_Parsed.threads, config.IBF_Parsed.fragment_size, 0, true };
					filter.filter = buildIBF(params);
				}

				catch (std::out_of_range& e)
				{
					throw ConfigReaderException(e.what());
				}

				TargetFilters.emplace_back(std::move(filter));
			}
		}

		return TargetFilters;
	}

}

// Goal is to test the methods parse_reads(), split() and classify_reads() from the header src/main/classify.hpp 

TEST_F(ClassifyTest, ParseReadsTest){

	interleave::TReads reads;
	uint64_t prefixLength; 

    //parse_reads("Listeria.fastq", reads, prefixLength);

	std::cout<< "There are no tests for the method parse_reads due to seqan's tests! " << '\n';
}


TEST_F(ClassifyTest, SplitTest){

	std::string s1 {"ACGCGA"};

	/*std::vector<std::string> testEqu {split(s1, 'C')};
	EXPECT_EQ("A", testEqu[0]);  

	std::vector<std::string> testEqu_ {split(s1, 'G')};
	EXPECT_EQ("G", testEqu[1]);  */


    
}


/**
 * @brief Step by step test for classifyTest methods (void method has no return value!)
 * 
 */

TEST_F(ClassifyTest, ClassifyReadsTest){

	// parse parameters struct
	ConfigReader config; 

	config = ConfigReader("/mnt/c/bug29/ReadBouncer/config.toml");
	config.parse_general();
	config.parse();
	
	classify_reads(config,getIBF_(config, true, false),getIBF_(config, false, true));

    std::cout<< "Classify Tests! " << '\n';
}