



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

TEST_F(ClassifyTest, ParseReadsTest){


    std::cout<< "Classify Tests! " << '\n';
    parse_reads("Listeria.fastq");
}

TEST_F(ClassifyTest, SplitTest){


    std::cout<< "Classify Tests! " << '\n';
}

TEST_F(ClassifyTest, ClassifyReadsTest){


    std::cout<< "Classify Tests! " << '\n';
}