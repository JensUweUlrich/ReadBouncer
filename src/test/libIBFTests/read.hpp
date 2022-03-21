



/* Submethods
1- get_threshold_kmers( uint16_t readLen, uint16_t kmerSize, float min_kmers )
2-  uint16_t get_error( uint16_t readLen, uint16_t kmerSize, uint16_t kmer_count)
3- uint16_t get_threshold_errors( uint16_t readLen, uint16_t kmerSize, uint16_t max_error )
4- double RationalApproximation(double t) --> Done
5- double NormalCDFInverse(double p) --> Done
6- TInterval calculateCI(const double r, const uint8_t kmer_size, const uint32_t readlen, const double confidence) --> Done
*/



class ReadTest: public ::testing::Test
{
	//MockIBF mock_ibf_;
	
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
		}
        // release any resources we allocated in SetUp()
		void TearDown() override
		{
			delete read;
		}

    public:// mocking
        

        

};

//bool find_matches(std::vector< TIbf >& filters, ClassifyConfig&  config );
/**
        find matches of kmers in ibfs
        return true if read matches at least one bin of a given IBF
        @matches:   unordered map to store matches 
        @filters:   vector of IBFs
        @config:    configuration settings needed
        @throws:    CountKmerException
    */

/**
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

TEST_F (ReadTest, FindMatchesTest){

    std::cout << "\n ########################################## " << '\n';
    std::cout << "\n Testing Class Interleave::Read! " << '\n';
	std::vector<interleave::IBFMeta> filters {};
    
    ClassifyConfig  config;
    config.error_rate = 0.1;
    config.significance = 0.95;
    config.strata_filter = -1;

    testing::NiceMock<MockRead> Mock_Read ;
    
    bool found = false;
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

            TIbf & holder = filter.filter;
            IBFs.emplace_back(std::move(f.getFilter()));
            
        }

    EXPECT_EQ(2, IBFs.size());
    ASSERT_TRUE(std::filesystem::exists("IbfClassificationLog.txt"));

    //Due to private method we test each step of the method then mock the results! 
    for ( TIbf& filter : IBFs )
        {
                std::vector< uint16_t > selectedBins    = seqan::count( filter, Mock_Read.sequence );
                std::vector< uint16_t > selectedBinsRev = seqan::count( filter, TSeqRevComp( Mock_Read.sequence ) );

                for(auto i: selectedBins){

                    //std::cout << "selectBins: " << i << std::endl; // found at 23 bins! 
                }

                for(auto i: selectedBinsRev){

                    //std::cout << "selectedBinsRev: " << i << std::endl; // found at 0 bins ---> correct! 
                }

                // get calculated threshold for minimum number of kmers needed to report a match
                // this is based on the confidence interval for mutated kmers in a read with expected error rate, kmer size
                TInterval ci = calculateCI(config.error_rate, filter.kmerSize, seqan::length(Mock_Read.sequence), config.significance);
                
                EXPECT_EQ(5, ci.first); // see calculation steps below! 
                EXPECT_EQ(30, ci.second);

                uint16_t readlen = seqan::length(Mock_Read.sequence);
                EXPECT_EQ(35, readlen);

                int16_t threshold = readlen - filter.kmerSize + 1 - ci.second;
                EXPECT_EQ(-7, threshold);

                //found = Mock_Read.select_matches( selectedBins, selectedBinsRev, filter, threshold);
                ASSERT_TRUE(Mock_Read.select_matches( selectedBins, selectedBinsRev, filter, threshold));
                
                //std::cout << "filter.noOfBins: " << filter.noOfBins << std::endl;
                for ( uint32_t binNo = 0; binNo < filter.noOfBins; ++binNo ) 
                {
                    //std::cout << "selectedBins[binNo]: " << selectedBins[binNo] << std::endl;
                    //std::cout << "selectedBinsRev[binNo]: " << selectedBinsRev[binNo] << std::endl;

                    if ( selectedBins[binNo] >= threshold || selectedBinsRev[binNo] >= threshold )
                    {
                        found = true;
                        EXPECT_EQ(found, true);
                    }
                }
        }


    /*
    {
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

        

    
    

}



