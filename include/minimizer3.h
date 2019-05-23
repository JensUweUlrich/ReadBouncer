#include <deque>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/view/all.hpp>
#include <seqan3/range/view/kmer_hash.hpp>
#include <seqan3/std/ranges>
using namespace seqan3;

struct Minimizer
{
public:

    // Random, but static value for xor for hashes. Counteracts consecutive minimizers.
    // E.g., without it, the next minimizer after a poly-A region AAAAA would be most likely something like AAAAC.
    // Setting seed to 0, will lead to lexicographically smallest k-mer
    uint64_t seed{0x8F3F73B5CF1C9ADE};
    // k-mer size
    uint16_t k{19};
    // window size != w in minimizer paper
    uint32_t w{25};
    // start positions of minimizers
    std::vector<uint64_t> minBegin;
    // end positions of minimizers
    std::vector<uint64_t> minEnd;


    inline void resize(uint16_t newKmerSize, uint32_t neww, uint64_t newSeed = 0x8F3F73B5CF1C9ADE)
    {
        k = newKmerSize;
        w = neww;
        seed = newSeed;
    }

    inline void setKmerSize(uint16_t newKmerSize)
    {
    	k = newKmerSize;
    }

    std::vector<uint64_t> getMinimizer(dna5_vector & text)
    {
        if (k > text.size())
            return std::vector<uint64_t> {};

        // Reverse complement without copying/modifying the original string
        dna5_vector revComp = text | std::view::reverse | view::complement;

        uint64_t possible = text.size() > w ? text.size() - w + 1 : 1; 	// number of all possible windows in text
        uint32_t windowKmers = w - k + 1;								// number of kmers in a window
        																// corresponds to w in minimizer paper
        uint64_t kmerSize = text.size() - k + 1;						// number of all possible k-mers in text

        std::vector<uint64_t> kmerHashes{};
        kmerHashes.reserve(possible); 									// maybe rather reserve to expected?

        // Stores hash, begin and end for all k-mers in the window
        std::deque<uint64_t> windowValues;

        std::vector<uint64_t> kmerHashIt = text | view::kmer_hash(k);
        std::vector<uint64_t> revcHashIt = revComp | view::kmer_hash(k);

        // Initialisation. We need to compute all hashes for the first window.
        for (uint32_t i = 0; i < windowKmers; ++i)
        {
            // Get smallest canonical k-mer
            uint64_t kmerHash = kmerHashIt[i] ^ seed;
            uint64_t revcHash = revcHashIt[kmerSize-i-1] ^ seed;
            if (kmerHash <= revcHash)
            {
                windowValues.push_back(kmerHash);

            }
            else
            {
                windowValues.push_back(revcHash);
            }
        }

        auto min = std::min_element(std::begin(windowValues), std::end(windowValues));
        kmerHashes.push_back(*min);

        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting
        bool minimizer_changed{false};
        for (uint64_t i = 1; i < possible; ++i)
        {
            if (min == std::begin(windowValues))
            {
                windowValues.pop_front();
                min = std::min_element(std::begin(windowValues), std::end(windowValues));
                minimizer_changed = true;
            }
            else
                windowValues.pop_front();

            uint64_t kmerHash = kmerHashIt[windowKmers + i - 1] ^ seed;
            uint64_t revcHash = revcHashIt[kmerSize - windowKmers - i] ^ seed;
            if (kmerHash <= revcHash)
                windowValues.push_back(kmerHash);
            else
                windowValues.push_back(revcHash);

            if (windowValues.back() < *min)
            {
                min = std::end(windowValues) - 1;
                minimizer_changed = true;
            }

            if (minimizer_changed)
            {
                kmerHashes.push_back(*min);
                minimizer_changed = false;
            }
        }

        return kmerHashes;
    }
};
