



/* Submethods
1- get_threshold_kmers( uint16_t readLen, uint16_t kmerSize, float min_kmers )
2-  uint16_t get_error( uint16_t readLen, uint16_t kmerSize, uint16_t kmer_count)
3- uint16_t get_threshold_errors( uint16_t readLen, uint16_t kmerSize, uint16_t max_error )
4- double RationalApproximation(double t)
5- double NormalCDFInverse(double p)
6- TInterval calculateCI(const double r, const uint8_t kmer_size, const uint32_t readlen, const double confidence)
*/

class ReadTest: public ::testing::Test
{
	//MockIBF mock_ibf_;
	
	protected:

		interleave::Read *read;
		typedef seqan::BinningDirectory< seqan::InterleavedBloomFilter,
                                    seqan::BDConfig< seqan::Dna5, seqan::Normal, seqan::Uncompressed > >
        TIbf;
        
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

TEST_F (ReadTest, FirstTest){

	std::cout << "Test 1" << std::flush; 

}