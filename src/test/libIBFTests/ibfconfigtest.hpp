/*
* Test class  
*
*/
class IBFConfigTest: public ::testing::Test
{	
	protected:

        interleave::IBFConfig *ibfConfig;
        
		// prepare the objects for each test.
		void SetUp() override
		{
			ibfConfig = new interleave::IBFConfig();
		}
        // release any resources we allocated in SetUp()
		void TearDown() override
		{
			delete ibfConfig;
		}

    public:// mocking

};

/**
 * @brief Test the definition variables of the class IBFConfig. IBFConfig is already tested (part of) in creating filter tests (class IBFTest).
 * This test is only for the initial values. 
 * 
 */

TEST_F (IBFConfigTest, VariableTest){

    ASSERT_EQ(8388608, ibfConfig->MBinBits);
    ASSERT_EQ(0, ibfConfig->filter_size);
    ASSERT_EQ(0, ibfConfig->filter_size_bits);
    ASSERT_EQ(0, ibfConfig->fragment_length);
    ASSERT_EQ(1500, ibfConfig->overlap_length);
    ASSERT_EQ(13, ibfConfig->kmer_size);
    ASSERT_EQ(3, ibfConfig->hash_functions);
    ASSERT_EQ(2, ibfConfig->threads);
    ASSERT_EQ(400, ibfConfig->n_refs);
    ASSERT_EQ(500000, ibfConfig->n_batches);
    ASSERT_EQ(0.01, ibfConfig->max_fp);
    ASSERT_EQ(1, ibfConfig->threads_build);
    
    if(ibfConfig->reference_files.empty() && ibfConfig->directory_reference_files.empty() && ibfConfig->extension.empty()
                                          && ibfConfig->output_filter_file.empty() && ibfConfig->input_filter_file.empty()
                                          && ibfConfig->update_filter_file.empty()){

    SUCCEED();

     }

    else{

        FAIL();
    }
                                         
}