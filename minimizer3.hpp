#include <deque>

#include <seqan3/core/debug_stream.hpp>
#include "MurmurHash3.h"
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/container/dynamic_bitset.hpp>
#include <seqan3/range/views/all.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/std/ranges>
using namespace seqan3;

class Minimizer
{
	public:


		// Random, but static value for xor for hashes. Counteracts consecutive minimizers.
		// E.g., without it, the next minimizer after a poly-A region AAAAA would be most likely something like AAAAC.
		// Setting seed to 0, will lead to lexicographically smallest k-mer
		uint64_t seed
		{ 0x8F3F73B5CF1C9ADE };
		// k-mer size -> set default value to 31 according to kraken2 paper
		uint8_t k
		{ 5 };
		// window size != w in minimizer paper -> set default value to 35 according to kraken2 paper
		uint8_t w
		{ 15 };
		// number of masked positions
		uint8_t s {0};
		// gapped shape
		seqan3::shape shape{};
		// start positions of minimizers
		std::vector<uint64_t> minBegin;
		// end positions of minimizers
		std::vector<uint64_t> minEnd;

		std::vector<dna4_vector> getMinimizer(dna4_vector & text, const bool ungapped);
		
		std::vector<uint64_t> getMinimizerHashValues(dna4_vector & text);


	   double NumberOfMinimizer;

	   Minimizer()
	   {
		   NumberOfMinimizer = 0.0;
	   }

		

		inline void resize(uint16_t newKmerSize, uint32_t neww, uint64_t newSeed)
		{
			k = newKmerSize;
			w = neww;
			seed = newSeed;
		}

		inline void setKmerSize(uint8_t newKmerSize)
		{
			k = newKmerSize;
		}
		
		inline void setWindowSize(uint8_t newWindowSize)
		{
			w = newWindowSize;
		}
		
		inline void setMaskedPositions(uint8_t newS)
		{
			s = newS;
			shape = computeSpacedSeed();
		}
		///////////////////////////////////////////////////// Changed //////////////////////////////////////////////
		//old 
		/*inline void setGappedShape(std::vector<int>)
		{
			std::vector<int> t2;
			
		}*/
		//new 
		
		inline void setGappedShape(seqan3::shape& newShape)
		{
			shape = newShape;
		}
		///////////////////////////////////////////////////// Changed //////////////////////////////////////////////
		inline seqan3::shape computeSpacedSeed()
		{
			// TODO: exchange with e.g. seqan3::shape s1{seqan3::bin_literal{0b101}}; // represents "101", i.e. gapped 3-mer
			// give parameter for number of masked positions -> every other position beginning at rightmost position -> see kraken2 paper with default of 7 masked positions
			seqan3::dynamic_bitset bs{};
			bs.resize(k);
			bs.set();
			
			uint8_t stepSize = ((uint8_t)(k / (s + 1)));
			
			for (int i = 1 ; i <= s; i++)
			{
				int pos =  stepSize * i;
				bs.flip(pos);
			}
			
			seqan3::shape shape1{};
			shape1.assign(bs.begin(), bs.end());
			return shape1;
		}


};
// double Minimizer::sum = 0.0; Arch error  deswegen da oben , wo cpp ist 

