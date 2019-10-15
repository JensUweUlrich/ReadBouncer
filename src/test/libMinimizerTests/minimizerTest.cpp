//#include "mockBloomFilter.hpp"
//#include "gmock/gmock.h"
#include "minimizer3.hpp"
#include "gtest/gtest.h"

/**
 * Test fixture class used to set up some ressources required for all tests
 */
class MinimizerTest: public ::testing::Test
{
	protected:

		Minimizer* mini;
		
		void SetUp() override
		{
			mini = new Minimizer();
			mini->setKmerSize(15);
			mini->setWindowSize(20);
			mini->setMaskedPositions(1);
		}

		void TearDown() override
		{
			delete mini;
		}

};

/**
 * Tests
 */


/**
 * Test writing and reading Bloom Filter to a binary output file
 */
TEST_F(MinimizerTest, TestSpacedSeedComputation)
{
	seqan3::dynamic_bitset t1{"111111101111111"};
	std::string st1 = t1.to_string();
	
	seqan3::shape t2 = mini->computeSpacedSeed();
	std::string st2 = t2.to_string();

	// fatal failure if seeds are unequal
	ASSERT_EQ(st1, st2);

}

TEST_F(MinimizerTest, TestGetMinimizer)
{
	seqan3::shape t2 = mini->computeSpacedSeed();
	seqan3::dna4_vector window{"GCATTATCGTGAAACGCTTT"_dna4};
	std::vector<uint64_t> kmerHashes = window | views::kmer_hash(t2) | ranges::to<std::vector<uint64_t>>();
	seqan3::debug_stream << kmerHashes << std::endl;
	
	
	std::vector<dna4_vector> minimizer = mini->getMinimizer(window, false);
	std::vector<char> min_char = minimizer[0] | seqan3::views::to_char | ranges::to<std::vector<char>>();
	std::string min(min_char.begin(), min_char.end());
	
	ASSERT_EQ(min, "TTATCGTGAAACGCT");
	
}

TEST_F(MinimizerTest, TestGetMinimizerWithSNP1)
{
	seqan3::shape t2 = mini->computeSpacedSeed();
	seqan3::dna4_vector window2{"GCATTATCGTGACACGCTTT"_dna4};
	std::vector<uint64_t> kmerHashes2 = window2 | views::kmer_hash(t2) | ranges::to<std::vector<uint64_t>>();
	seqan3::debug_stream << kmerHashes2 << std::endl;
}

TEST_F(MinimizerTest, TestGetMinimizerWithSNP2)
{
	mini = new Minimizer();
	mini->setKmerSize(18);
	mini->setWindowSize(25);
	seqan3::shape t2{0b111010010100110111_shape};
	mini->setGappedShape(t2);
	//seqan3::shape t2 = mini->computeSpacedSeed();
	seqan3::dna4_vector window1{"GCATTATCGTGAAACGCTTTCGCGTTTTCGTGCGCCGCTTGATAACAAATCCCTGCTTCAGGAAAAGCCTCTAAGCAT"_dna4};
	std::vector<uint64_t> kmerHashes1 = window1 | views::kmer_hash(t2) | ranges::to<std::vector<uint64_t>>();
	//seqan3::debug_stream << kmerHashes1 << std::endl;
	
	//mini->setMaskedPositions(15);
	std::vector<uint64_t> minimizer = mini->getMinimizerHashValues(window1);
	for (int i = 0;i<minimizer.size();++i)
	{
		minimizer[i] ^= mini->seed;
	}
	debug_stream << minimizer << std::endl;
	
	
	seqan3::dna4_vector window2{"GCATTATCGTGAAACGCTTTCGCGTTTCGTGCGCCGCTTGATAACAAAATCCCTGCTTCAGGAAAAAGCCTCTAAGCAT"_dna4};
	std::vector<uint64_t> kmerHashes2 = window2 | views::kmer_hash(t2) | ranges::to<std::vector<uint64_t>>();
	//seqan3::debug_stream << kmerHashes2 << std::endl;
	minimizer = mini->getMinimizerHashValues(window2) ;
	for (int i = 0;i < minimizer.size();++i)
	{
		minimizer[i] ^= mini->seed;
	}
		
	debug_stream << minimizer << std::endl;
	
}



int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
