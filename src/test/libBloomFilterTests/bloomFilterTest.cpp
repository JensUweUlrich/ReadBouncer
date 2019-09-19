//#include "mockBloomFilter.hpp"
//#include "gmock/gmock.h"
#include "customBloomFilter.hpp"
#include "gtest/gtest.h"

namespace fs = std::experimental::filesystem;

/**
 * Test fixture class used to set up some ressources required for all tests
 */
class BloomFilterTest: public ::testing::Test
{
	protected:

		CustomBloomFilter * bf1;
		fs::path bf_tmp_file;

		void SetUp() override
		{
			bf1 = new CustomBloomFilter(0.05, 600, 19);
			fs::path tmp_dir = fs::temp_directory_path();
			fs::path f ("tmp_bloom_filter.temp");
			bf_tmp_file = tmp_dir / f;
		}

		void TearDown() override
		{
			delete bf1;
			if (std::experimental::filesystem::exists(bf_tmp_file))
			{
				fs::remove(bf_tmp_file);
			}
		}

};

/**
 * Tests
 */

/**
 * Test for correct computation of optimal Bloom Filter size
 */
TEST_F(BloomFilterTest, TestBloomFilterSize)
{
	EXPECT_EQ(3742, bf1->bits.size());
}

/**
 * Test for correct computation of optimal hash function number
 */
TEST_F(BloomFilterTest, TestHashFunctionNumber)
{
	EXPECT_EQ(5, bf1->hashes.size());
}

/**
 * Test for correctly setting Kmer size in Bloom Filter object
 */
TEST_F(BloomFilterTest, TestKmerSize)
{
	EXPECT_EQ(19, bf1->kMerSize);
}

/**
 * Test writing and reading Bloom Filter to a binary output file
 */
TEST_F(BloomFilterTest, TestWriteBloomFilterToFile)
{
	bf1->writeToFile(bf_tmp_file);
	// fatal failure if file does not exist
	ASSERT_TRUE(std::experimental::filesystem::exists(bf_tmp_file));

	CustomBloomFilter bf2
	{ };
	// fatal failure if read from file failed
	ASSERT_TRUE(bf2.readFromFile(bf_tmp_file));

	// fatal failure if Kmer size is unequal
	ASSERT_EQ(bf1->kMerSize, bf2.kMerSize);

	for (int i = 0; i < bf1->hashes.size(); ++i)
	{
		ASSERT_EQ(bf1 -> hashes.at(i), bf2.hashes.at(i)) << "hash functions unequal at position " << i << "\n";
	}

	for (int i = 0; i < bf1->bits.size(); ++i)
	{
		ASSERT_EQ(bf1->bits.at(i), bf2.bits.at(i)) << "bit vectors unequal at position " << i << "\n";
	}

}

int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
