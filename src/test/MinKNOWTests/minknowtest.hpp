/**
 * minknowgtests.hpp
 * 
 */

class MinknowTest: public ::testing::Test
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

TEST_F(MinknowTest, fragment_start){
	
	std::cout<< "MinKNOW Tests....."<< std::endl;
	//EXPECT_EQ(fragment_start(config.IBF_Parsed.chunk_length, 1), 360);
	//EXPECT_EQ(fragment_start(config.IBF_Parsed.chunk_length, 2), 720);
	//EXPECT_EQ(fragment_start(450, 2), 900);
}