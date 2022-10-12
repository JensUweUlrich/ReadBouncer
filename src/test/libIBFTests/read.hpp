



class ReadTest: public ::testing::Test
{
	
	
	protected:

		interleave::Read *read;
		typedef seqan::BinningDirectory< seqan::InterleavedBloomFilter,
                                    seqan::BDConfig< seqan::Dna5, seqan::Normal, seqan::Uncompressed > >
        TIbf;

        std::vector <TIbf> IBFs;
        
		// prepare the objects for each test.
		void SetUp() override
		{
			read = new interleave::Read();
            read->sequence = "AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAGAGCCCCAAAAGAGAGGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAGAGCCCCAAAAGAGAGGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAGAGCCCCAAAAGAGAGGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAGAGCCCCAAAAGAGAGGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAGAGCCCCAAAAGAGAGGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAGAGCCCCAAAAGAGAGGAGA";
            
            
            
		}
        // release any resources we allocated in SetUp()
		void TearDown() override
		{
			delete read;
		}

        template <typename T>
        bool pairEq(std::pair<T, T> , std::pair<T, T>);


    public:// mocking
        int testCounter, TC
	    {0};
    
        
        

};

/**
 * @test
 * test.ibf: AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAGAGCCCCAAAAGAGAGGAGATTTTANNNNNNNNTATATTATA
 * test1.ibf:
 >TestReferenceSequence_1
AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGA
AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGA
AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGA
AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGA
AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGA
AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGA

>TestReferenceSequence_2
AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAA
AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAA
AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAA
AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAA
AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAA
AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAAAAAAAAA
 * 
 */

/*
* Compare two STD pairs. 
* @return true if both pairs are equal
*/

template <typename T>
bool ReadTest::pairEq(std::pair<T, T> p1, std::pair<T, T> p2){

    if((p1.first == p2.first) && (p1.second == p2.second)){

        return true;
    }

    else{

        return false; 
    }

}

/*
* Test the main classification methods, some methods are tested line by line as errors might occur (hidden) in correct results
*/


TEST_F (ReadTest, ClassificationTest){

    std::cout << "\n ########################################## " << '\n';
    std::cout << "Testing Class Interleave::Read! " << '\n';

    const testing::TestInfo* const test_info = testing::UnitTest::GetInstance()->current_test_info();
    printf("We are in test %s of test suite %s.\n", test_info->name(), test_info->test_suite_name());

	std::vector<interleave::IBFMeta> filters {};

    interleave::ClassifyConfig config{};
    config.error_rate = 0.1;
    config.significance = 0.95;
    config.strata_filter = -1;

    testing::NiceMock<MockRead> Mock_Read ;
    
    bool found, found_ = false;
    std::vector<std::filesystem::path> filters_files = {"/mnt/c/bug29/ReadBouncer/build/test/libIBFTests/test.ibf", "/mnt/c/bug29/ReadBouncer/build/test/libIBFTests/test1.ibf"}; 

    Mock_Read.sequence = "AAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAG";
    
    

    for (std::filesystem::path file : filters_files)
		{
			interleave::IBFMeta filter{};
            filter.name = file.stem().string();
			interleave::IBF f{};
			interleave::IBFConfig FilterIBFconfig{};

			FilterIBFconfig.input_filter_file = file.string();
			interleave::FilterStats stats = f.load_filter(FilterIBFconfig);
			filter.filter = std::move(f.getFilter());
			filters.emplace_back(std::move(filter));

            IBFs.emplace_back(std::move(f.getFilter()));
            
        }

    EXPECT_EQ(2, IBFs.size());
    ASSERT_TRUE(std::filesystem::exists("IbfClassificationLog.txt"));
    testCounter += 2;

    for ( TIbf& filter : IBFs )
        {
                std::vector< uint16_t > selectedBins    = seqan::count( filter, Mock_Read.sequence );
                std::vector< uint16_t > selectedBinsRev = seqan::count( filter, TSeqRevComp( Mock_Read.sequence ) );

                /*for(auto i: selectedBins){

                    //std::cout << "selectBins: " << i << std::endl; // found at 23 bins! 
                }*/

                /*
                for(auto i: selectedBinsRev){

                    std::cout << "selectedBinsRev: " << i << std::endl; // found at 0 bins ---> correct! 
                }*/

                
                TInterval ci = calculateCI(config.error_rate, filter.kmerSize, seqan::length(Mock_Read.sequence), config.significance);
                
                EXPECT_EQ(5, ci.first); // see calculation steps below! 
                EXPECT_EQ(30, ci.second);
                testCounter += 2;

                uint16_t readlen = seqan::length(Mock_Read.sequence);
                EXPECT_EQ(35, readlen);

                int16_t threshold = readlen - filter.kmerSize + 1 - ci.second;
                EXPECT_EQ(-7, threshold);
                testCounter++;

                found =  Mock_Read.select_matches ( selectedBins, selectedBinsRev, filter, threshold);
                //found_ = read->select_matches_test( selectedBins, selectedBinsRev, filter, threshold);

                EXPECT_EQ(found, found_);
                testCounter++;
                
                for ( uint32_t binNo = 0; binNo < filter.noOfBins; ++binNo ) // select_matches! 
                {

                    if ( selectedBins[binNo] >= threshold || selectedBinsRev[binNo] >= threshold )
                    {
                        found = true;
                        EXPECT_EQ(found, true);
                        testCounter++;
                    }
                }
        }

   std::vector <TIbf> emptyVector;
   if (emptyVector.empty())// test exception
    {
       EXPECT_THROW(read->classify(emptyVector, config), NullFilterException);
       testCounter++;
    }

    if (read->getReadLength() < IBFs[0].kmerSize)
    {
       EXPECT_THROW(read->classify(emptyVector, config), CountKmerException);
       EXPECT_THROW(read->classify(emptyVector, config), ShortReadException);
       testCounter += 2;
    }

    EXPECT_EQ(13, IBFs[0].kmerSize);
    EXPECT_EQ(354, read->getReadLength());
    //EXPECT_EQ(1,  read->find_matches_test(IBFs, config));// already tested with short read and we had a very big threshold --> fixed
    EXPECT_EQ(true,  read->classify(IBFs, config));
    testCounter += 4;

    std::vector<interleave::IBFMeta> emptyVectorMeta;
    if (emptyVectorMeta.empty())// test exception
    {
        EXPECT_THROW(read->classify(emptyVectorMeta, config), NullFilterException);
        testCounter++;
    }

    if (read->getReadLength() < IBFs[0].kmerSize)// smaller than the k-mer size! 
    {
       EXPECT_THROW(read->classify(emptyVectorMeta, config), CountKmerException);
       EXPECT_THROW(read->classify(emptyVectorMeta, config), ShortReadException);
       testCounter += 2;
    }

    //int classifyReturnValue = read->classify(filters, config);

    std::vector<int> matchingResults= {282, 182};
    ASSERT_EQ(filters.size(), matchingResults.size());
    testCounter++;

    for (uint8_t ind = 0; ind < filters.size(); ++ind)
        {
            //ASSERT_EQ(matchingResults[ind], read->count_matches_test(filters[ind], config));
            testCounter++;
        }
    
    ASSERT_EQ(0, read->classify(filters, config));
    //EXPECT_EQ(282, read->count_matches_test(filters[read->classify(filters, config)], config));
    testCounter += 2;

    //std::pair<int, int> Read::classify(std::vector< IBFMeta >& filt1, std::vector< IBFMeta >& filt2, ClassifyConfig& config)

    std::vector<interleave::IBFMeta> emptyM1;
    std::vector<interleave::IBFMeta> emptyM2;
    if (emptyM1.empty() || emptyM2.empty())// test exception
    {
        EXPECT_THROW(read->classify(emptyM1, emptyM2, config), NullFilterException);
        testCounter++;
    }

    std::vector<interleave::IBFMeta> v1 {};
    std::vector<interleave::IBFMeta> v2 {};
    v1.emplace_back(filters[0]);
    v2.emplace_back(filters[1]);
    
    std::pair <int, int> p1  (282, 182);
    ASSERT_EQ(true, pairEq(p1, read->classify(v1,v2, config)));
    testCounter++;


};


/*
* Test of the count_matches method. Due to the usual use of this method, it is tested line by line alone! 
*/

//uint64_t count_matches(IBFMeta& filter, ClassifyConfig& config);
TEST_F (ReadTest, CountMatchesTest){

    testing::NiceMock<MockRead> Mock_Read ;

    const testing::TestInfo* const test_info = testing::UnitTest::GetInstance()->current_test_info();
    printf("\n We are in test %s of test suite %s.\n", test_info->name(), test_info->test_suite_name());

    std::string file = "/mnt/c/bug29/ReadBouncer/build/test/libIBFTests/test.ibf"; 

    //Mock_Read.sequence = "AAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAG";
    seqan::Dna5String RevCom = "CTCTCCTCTCTCCTCTCTCGGGGGGGGGTTTTTTT"; 

    interleave::IBFMeta filter{};
    interleave::IBF f{};
    interleave::IBFConfig FilterIBFconfig{};

    // fill struct with filter! 
    FilterIBFconfig.input_filter_file = file;
    f.load_filter(FilterIBFconfig);
    filter.filter = std::move(f.getFilter());
    
    // Count_matches:

    std::vector< uint16_t > selectedBins = seqan::count(filter.filter, RevCom);
    std::vector< uint16_t > selectedBinsRev = seqan::count(filter.filter, TSeqRevComp(RevCom));

    /*
    for(auto i: selectedBins){

          //std::cout << "selectBins: " << i << std::endl; // found at 0 bins! 
    }

    for(auto i: selectedBinsRev){

        //std::cout << "selectedBinsRev: " << i << std::endl; // found at 23 bins ---> correct! 
    }*/

    TInterval ci = calculateCI(0.1, 13, seqan::length(RevCom), 0.95);
    EXPECT_EQ(5, ci.first); 
    EXPECT_EQ(30, ci.second);
    testCounter += 2;

    uint16_t readlen = seqan::length(RevCom);
    EXPECT_EQ(35, readlen);
    testCounter++;

    int16_t threshold = readlen - filter.filter.kmerSize + 1 - ci.second;
    EXPECT_EQ(-7, threshold);
    testCounter++;

    uint64_t max_kmer_count = 0;

        for (uint32_t binNo = 0; binNo < filter.filter.noOfBins; ++binNo)
        {

            if (selectedBins[binNo] >= threshold || selectedBinsRev[binNo] >= threshold)
            {
                if (selectedBins[binNo] > max_kmer_count)
                    max_kmer_count = selectedBins[binNo];
                if (selectedBinsRev[binNo] > max_kmer_count)
                    max_kmer_count = selectedBinsRev[binNo];
            }
        }

    EXPECT_EQ(23, max_kmer_count);
    testCounter++;

    
    
}

/*
    { @calculations Example: 

         calculateCI(config.error_rate, filter.kmerSize, seqan::length(Mock_Read.sequence), config.significance);
         r = error_rate = 0.1 
         kmer_size = 13
         readlen = 35
         confidence = 0.95 
         
        double q = 1.0 - pow(1.0 - r, kmer_size) = 1 - pow(1-0.1. 13) = 0.7458134
        double L = ((double)readlen - (double)kmer_size + 1.0) = 35 - 13 + 1 = 23
        double Nmut = L * q = 17.1537
        double varN = L * (1.0 - q) * (q * (2.0 * (double)kmer_size + (2.0 / r) - 1.0) - 2.0 * (double)kmer_size)      --> 44.19227 
                        + (double)kmer_size * ((double)kmer_size - 1.0) * pow((1.0 - q), 2.0)                          --> 10.079
                        + (2.0 * (1.0 - q) / (pow(r, 2.0))) * ((1.0 + ((double)kmer_size - 1.0) * (1.0 - q)) * r - q); --> -17.331 
               varN = 36.9401873 
        double alpha = 1 - confidence = 1 - 0.95 = 0.05
        
        double z = NormalCDFInverse(1.0 - alpha / 2.0) = NormalCDFInverse(0.975) --> z = 1.96

        uint16_t low = (uint16_t) floor(L * q - z * sqrt(varN)) = floor(23*0.7458134 - 1.96*sqrt(36.9401873)) = floor(17.153 - 11.898) = floor(5.255) = 5
        uint16_t high = (uint16_t)ceil(L * q + z * sqrt(varN))  = ceil(23*0.7458134  + 1.96*sqrt(36.9401873)) =  ceil(17.153 + 11.898)  = ceil(29.051) = 30

        #varN 

        23 * (1.0 - 0.7458134) * (0.7458134 * (2.0 * 13 + (2.0 / 0.1) - 1.0) - 2.0 * 13) = 
        23 * (0.2541) * (0.7458134 * (26 + 20 - 1.0) - 26) =
        23 * (0.2541) * (0.7458134 * (26 + 20 - 1.0) - 26) =
        23 * (0.2541) * (0.7458134 * (45) - 26) =
        23 * (0.2541) * (33.5616 - 26) =
        23 * (0.2541) * (7.56160) = 44.19227

        13 * (13 - 1.0) * pow((1.0 - 0.7458134), 2.0) = 
        13 * (12) * pow((1.0 - 0.7458134), 2.0) = 
        13 * (12) * pow((0.25418), 2.0) = 
        13 * (12) * 0.0646 = 10.079

        (2.0 * (1.0 - 0.7458134) / (pow(0.1, 2.0))) * ((1.0 + (13 - 1.0) * (1.0 - 0.7458134)) * 0.1 - 0.7458134) = 
        (2.0 * (0.2541) / (0.01)) * ((1.0 + (12) * (0.2541)) * 0.1 - 0.7458134) = 
        (2.0 * 25.4186) * ((1.0 + 3.049) * 0.1 - 0.7458134) = 
        50.8372 * (4.049 * 0.1 - 0.7458134) = 
        50.8372 * (- 0.34) = -17.331 

        #NormalCDFInverse(0.975)
        RationalApproximation(sqrt(-2.0 * log(1.0 - 0.975))) =                    ** @Note: log() is natural logarithm! **
        RationalApproximation(sqrt(-2.0 * log(0.025))) = 
        RationalApproximation(sqrt(-2.0 *-3.68888)) =
        RationalApproximation(sqrt(7.3776)) =
        RationalApproximation(2.716) =

        2.716 - ((0.010328 * 2.716 + 0.802853) * 2.716 + 2.515517) / (((0.001308 *2.716 + 0.189269) * 2.716 + 1.432788) * 2.716 + 1.0) = 
        2.716 - ((0.8309) * 2.716 + 2.515517) / (((0.19282) * 2.716 + 1.432788) * 2.716 + 1.0) =
        2.716 - (4.7722) / ((1.9564) * 2.716 + 1.0) =
        2.716 - (4.7722) / (6.3138) =
        2.716 - 0.7558 = 1.96
        

    }

        */