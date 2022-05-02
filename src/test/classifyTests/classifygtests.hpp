/**
 * ckassifygtests.hpp
 * 
 */

class ClassifyTest: public ::testing::Test
{
	//MockIBF mock_ibf_;
	
	protected:
        
		// prepare the objects for each test.
		ConfigReader config;

		void SetUp() override
		{
			// parse parameters struct 
			config = ConfigReader("/mnt/c/bug29/ReadBouncer/config.toml");

			config.parse_general();
			config.parse();

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


/**
 * @brief Test fragment_start method
 * 
 */

TEST_F(ClassifyTest, fragment_start){
	
	EXPECT_EQ(fragment_start(config.IBF_Parsed.chunk_length, 1), 360);
	EXPECT_EQ(fragment_start(config.IBF_Parsed.chunk_length, 2), 720);
	EXPECT_EQ(fragment_start(450, 2), 900);
}

/**
 * @brief Test fragment_end method
 * 
 */

TEST_F(ClassifyTest, fragment_end){
	
	EXPECT_EQ(fragment_end(config.IBF_Parsed.chunk_length, 1), 720);
	EXPECT_EQ(fragment_end(400, 1), 800);
	EXPECT_EQ(fragment_end(450, 2), 1350);
}

/**
 * @brief Test classification results
 * 
 */

TEST_F(ClassifyTest, ClassifyReadsTest){

	classify_reads(config,getIBF(config, true, false),getIBF(config, false, true));

	ASSERT_EQ(ClassificationResults_.found, 3);
	ASSERT_EQ(ClassificationResults_.failed, 0);
	ASSERT_EQ(ClassificationResults_.too_short, 0);
	ASSERT_EQ(ClassificationResults_.readCounter, 3);

}

/**
 * @brief Test unclassified reads
 * 
 */

TEST_F(ClassifyTest, writeUnclassifed){

	ASSERT_EQ(writeUnclassifed(config).string(), "RB_out/unclassified.fasta");
	
}