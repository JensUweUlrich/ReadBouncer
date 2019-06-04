#include <deque>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/view/all.hpp>
#include <seqan3/range/view/kmer_hash.hpp>
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
		// k-mer size
		uint16_t k
		{ 19 };
		// window size != w in minimizer paper
		uint32_t w
		{ 25 };
		// start positions of minimizers
		std::vector<uint64_t> minBegin;
		// end positions of minimizers
		std::vector<uint64_t> minEnd;

		std::vector<uint64_t> getMinimizer(dna5_vector & text);

		inline void resize(uint16_t newKmerSize, uint32_t neww, uint64_t newSeed)
		{
			k = newKmerSize;
			w = neww;
			seed = newSeed;
		}

		inline void setKmerSize(uint16_t newKmerSize)
		{
			k = newKmerSize;
		}

};
