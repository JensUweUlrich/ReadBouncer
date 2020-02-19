//#include "mockBloomFilter.hpp"
//#include "gmock/gmock.h"
#include "bloom_filter.hpp"
#include "customBloomFilter.hpp"
#include "gtest/gtest.h"

namespace fs = std::filesystem;

/**
 * Test fixture class used to set up some ressources required for all tests
 */
class BloomFilterTest: public ::testing::Test
{
	protected:

		CustomBloomFilter * bf1;
		fs::path bf_tmp_file;
		std::vector<int> testData {6874, 4654641, 69961, 38945, 685154, 346843149, 6494, 697451, 69851, 6516847, 36541694, 61941, 98515, 194168, 6541154, 68541, 88526, 698524, 9664};

		void SetUp() override
		{
			bloom_parameters parameters;

			// How many elements roughly do we expect to insert?
			parameters.projected_element_count = testData.size();

			// Maximum tolerable false positive probability? (0,1)
			parameters.false_positive_probability = 0.0001; // 1 in 10000

			// Simple randomizer (optional)
			parameters.random_seed = 0xA5A5A5A5;

			// compute optimal bloom filter parameters
			parameters.compute_optimal_parameters();
			
			// initialize bloom filter object
			bf1 = new CustomBloomFilter(parameters, 19);
			
			// insert test data
			bf1->insert(testData.begin(), testData.end());
			
			fs::path tmp_dir = fs::temp_directory_path();
			fs::path f ("tmp_bloom_filter.temp");
			bf_tmp_file = tmp_dir / f;
		}

		void TearDown() override
		{
			delete bf1;
			if (std::filesystem::exists(bf_tmp_file))
			{
				fs::remove(bf_tmp_file);
			}
		}

};

/**
 * Tests
 */


/**
 * Test writing and reading Bloom Filter to a binary output file
 */
TEST_F(BloomFilterTest, TestWriteBloomFilterToFile)
{
	bf1->writeToFile(bf_tmp_file);
	// fatal failure if file does not exist
	ASSERT_TRUE(std::filesystem::exists(bf_tmp_file));

	CustomBloomFilter bf2
	{ };
	// fatal failure if read from file failed
	ASSERT_TRUE(bf2.readFromFile(bf_tmp_file));

	// fatal failure if Kmer size is unequal
	ASSERT_EQ(bf1->getKmerSize(), bf2.getKmerSize());

	// fatal failure number of hash seeds is unequal
	ASSERT_EQ(bf1->hash_count(), bf2.hash_count());

	for (int i = 0; i < bf1->hash_count(); ++i)
	{
		ASSERT_EQ(bf1 -> getHashSeeds().at(i), bf2.getHashSeeds().at(i)) << "hash seeds unequal at position " << i << "\n";
	}

	// fatal failure if bloom filter size is unequal
	ASSERT_EQ(bf1->size(), bf2.size());

	for (int i = 0; i < bf1->getBitTable().size(); ++i)
	{
		ASSERT_EQ(bf1->getBitTable().at(i), bf2.getBitTable().at(i)) << "bit vectors unequal at position " << i << "\n";
	}

}

int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
