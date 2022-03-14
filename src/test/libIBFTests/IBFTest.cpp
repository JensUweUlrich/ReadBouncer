#include "IBF.hpp"
#include "gtest/gtest.h"

#include "mockIBF.hpp"


class IBFTest: public ::testing::Test
{
	//MockIBF mock_ibf_;
	
	protected:

		interleave::IBF *ibf;
		typedef seqan::BinningDirectory< seqan::InterleavedBloomFilter,
                                    seqan::BDConfig< seqan::Dna5, seqan::Normal, seqan::Uncompressed > >
        TIbf;
        
		// prepare the objects for each test.
		void SetUp() override
		{
			ibf = new interleave::IBF();
		}
        // release any resources we allocated in SetUp()
		void TearDown() override
		{
			delete ibf;
		}

    public:// mocking
	   // IBFTest(interleave::IBF *_ibf) : ibf(_ibf) {}



};



inline bool operator==(const std::future <void> & future_1, const std::future <void> & future_2)
{
    return future_1 == future_2;
}

inline bool operator!=(const std::future <void> & future_1, const std::future <void> & future_2)
{
    return future_1 != future_2;
}

inline bool operator==(const TIbf & filter_1, const TIbf & filter_2)
{
    return filter_1 == filter_2;
}

// Test create_filter method with all steps from parse_ref_seqs(), cutOutNNNs() and add_sequences_to_filter() by mocking!
TEST_F (IBFTest, CreateFilterTest){

	interleave::IBFConfig config{};
	int testCounter
	{0};

	config.reference_files.emplace_back("test.fasta");
	config.output_filter_file = "test.ibf";
	config.kmer_size = 13;
	config.threads_build = 1;
	config.fragment_length = 100000;

	EXPECT_EQ(true, config.validate());// as we defined a valid IBFConfig config object
	testCounter++;

	
	FilterStats stats = ibf->create_filter(config);// first call

	if(ibf->test_read_task.valid()){// check if the method IBF::parse_ref_seqs returns a valid future

		SUCCEED();
		testCounter++;
	}
	else{

		FAIL();
		EXPECT_THROW(ibf->create_filter(config), InvalidConfigException);
		testCounter++;
	}

	// Test each step in the method IBF::parse_ref_seqs()
	if (config.reference_files.empty())
       {
        EXPECT_THROW(ibf->create_filter(config), MissingReferenceFilesException);
		testCounter++;
       }

	std::stringstream buf;
	std::stringstream buf_test;
	buf_test<<"AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAGAGCCCCAAAAGAGAGGAGATTTTATATATTATA";

	for ( std::string const& reference_file : config.reference_files )
            {
                seqan::SeqFileIn seqFileIn;
                if ( !seqan::open( seqFileIn, seqan::toCString( reference_file ) ) )
                {
					EXPECT_THROW(ibf->create_filter(config), FileParserException);
					testCounter++;
                }
				//EXPECT_EQ(seqan::toCString(reference_file), "0x7f570696c14b");

				 // read current input file
                while ( !seqan::atEnd( seqFileIn ) )
                {

                    seqan::StringSet< seqan::CharString > ids;
                    seqan::StringSet< seqan::CharString > seqs;
					seqan::Dna5String seq; 
					seqan::CharString id; 
					seqan::CharString id_test = "TestReferenceSequence ";
					
                    // read all sequences from the current file
                    seqan::readRecords( ids, seqs, seqFileIn, config.n_refs );
					for(auto i : seqs){

						seq += i;
					}

					for(auto i : ids){

						id += i;
					}
					EXPECT_EQ(id, id_test);
					EXPECT_EQ(seq, "AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAGAGCCCCAAAAGAGAGGAGATTTTANNNNNNNNTATATTATA");
					testCounter += 2;
					//std::cout<<"Sequence length: "<< seqan::length(seqs[0]) << std::endl; 
					if (seqan::length(seqs) < config.kmer_size){

						stats.invalidSeqs += 1;
						EXPECT_TRUE(stats.invalidSeqs > 0);
						testCounter++;
					}
					else{

						EXPECT_EQ(0, stats.invalidSeqs);
						testCounter++;
					}
					
				
                }
				
			}
		//std::cout<<ibf->cutOutNNNsTest<<std::endl;
		EXPECT_EQ(ibf->cutOutNNNsTest, "AAAAAAAACCCCCCCCCGAGAGAGGAGAGAGGAGAGAGAGAGCCCCAAAAGAGAGGAGATTTTATATATTAT");// test cutOutNNNs method! 
		EXPECT_EQ(2, stats.totalBinsBinId); //needed bins: (seq.length()/fragment_length +1)
		EXPECT_EQ(144, stats.sumSeqLen);
		testCounter += 3;

		/* Count the filter size in bits: 
		* Formula: 
		uint64_t optimalNumberOfBins = floor(((double) 2 / 64.0) + 1) * 64; --> 64
		uint64_t BinSizeBits = ceil(-1 / (pow(1 - pow((double) 0.01, 1.0 / (double) 3),
                                                        1.0 / ((double) (3 * 99988))) - 1)); --> 1236269
		uint64_t x = BinSizeBits * optimalNumberOfBins; --> 79121216
		*/

		EXPECT_EQ(79121216, config.filter_size_bits);
		testCounter++;

		TIbf testFilter = TIbf(2, 3, 13, 79121216);

		//EXPECT_EQ(testFilter, ibf->getFilter_());
		EXPECT_EQ(seqan::getNumberOfBins(testFilter), seqan::getNumberOfBins(ibf->getFilter_()));
		testCounter++;

		EXPECT_EQ(3, config.hash_functions);//Number of hash functions
		testCounter++;

		// Test: add_sequences_to_filter(tasks, config, binid, queue_refs); here we use mocking, because this method is a private method! 
		testing::NiceMock<MockIBF> mock_ibf;
		uint64_t binid = 0; // reference to a binid counter that is incremented for every new fragment
		std::vector< std::future< void > > tasks; //empty vector of tasks, which needs to be fuilt 
		SafeQueue< Seqs > queue_refs(config.n_batches ); //thread safe queue that stores reference sequences to put in the ibf
		//config: is the settings used to build the IBF, already defined 
		mock_ibf.add_sequences_to_filter(tasks, config, binid, queue_refs);

		EXPECT_EQ(1, ibf->test_func1.threads_build);
		EXPECT_EQ(1, ibf->test_func1.fragIdx);
		EXPECT_EQ(99988, ibf->test_func1.fragstart); //fragstart = fragIdx * config.fragment_length - config.kmer_size + 1 = 1 * 100000 - 13 + 1 = 99988
		EXPECT_EQ(72, ibf->test_func1.fragend);/* uint64_t fragend = (fragIdx+1) * config.fragment_length = 200000
                                                * make sure that last fragment ends at last position of the reference sequence
                                                * if (fragend > length(val.seq)) fragend = length(val.seq); --> fragend = 72 */
	    testCounter += 4; 

		for ( auto& task : tasks )
        {
            EXPECT_EQ(true, task.valid());// as we defined a valid IBFConfig config object
	        testCounter++;
        }

		std::cout << "'\n' Total number of Tests for creating filter is: "<< testCounter << '\n';
}

// Test create_filter stats after creating the filter! 
TEST_F (IBFTest, FilterStatsTest){

	interleave::IBFConfig config{};
	config.reference_files.emplace_back("test.fasta");
	config.output_filter_file = "test.ibf";
	config.kmer_size = 13;
	config.threads_build = 1;
	config.fragment_length = 100000;
	
	FilterStats stats = ibf->create_filter(config);


}


TEST_F (IBFTest, mockTest){

	// arrange
	interleave::IBFConfig config{};
	MockIBF mock_ibf;

	config.reference_files.emplace_back("test.fasta");
	config.output_filter_file = "test.ibf";
	config.kmer_size = 13;
	config.threads_build = 1;
	config.fragment_length = 100000;

	//mock_ibf.create_filter(config);
	
	/*add_sequences_to_filter(tasks, config, binid, queue_refs);
            for ( auto& task : tasks )
            {
                task.get();
            }*/
	std::vector< std::future< void > > tasks;
	SafeQueue< Seqs > queue_refs(config.n_batches ); 
	uint64_t binid = 0;
	//mock_ibf.add_sequences_to_filter(tasks, config, binid, queue_refs);
	

	//EXPECT_EQ(true, config.validate());// as we defined a valid IBFConfig config object
	//testCounter++;
	//FilterStats stats = ibf->create_filter(config);// first call

	



}

int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
