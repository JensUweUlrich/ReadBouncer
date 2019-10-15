#include "minimizer3.hpp"
#include <math.h>
using namespace seqan3;

std::vector<dna4_vector> Minimizer::getMinimizer(dna4_vector & text, const bool ungapped)
{
	// no minimizers if text to short
	if (k > text.size())
		return std::vector<dna4_vector>
		{ };

	// Reverse complement without copying/modifying the original string
	dna4_vector revComp = text | std::views::reverse | views::complement | ranges::to<std::vector<seqan3::dna4>>();

	// computation of some parameters
	uint64_t possible = text.size() > w ? text.size() - w + 1 : 1; 	// number of all possible windows in text
	uint32_t windowKmers = w - k + 1;								// number of kmers in a window
	                                                                // corresponds to w in minimizer paper
	uint64_t kmerSize = text.size() - k + 1;						// number of all possible k-mers in text

	// resulting container with minimizer sequences that will be returned
	std::vector<dna4_vector> minimizers
	{ };
	// reserve space for maximum number of possible minimizer
	minimizers.reserve(possible); 									// maybe rather reserve to expected?

	// Stores hash values for all k-mers in the current window
	std::deque<uint64_t> windowValues;
	// stores all kmers corresponding to stored hash values in the current window
	std::deque<dna4_vector> windowValueKmers;

	// ungapped shape or masked shape used for kmer hashing?
	if (shape.size() != k)
	{
		seqan3::shape shape = computeSpacedSeed();
		if (ungapped)
		{
			shape = seqan3::ungapped
			{ k };
		}
	}
	std::vector<uint64_t> kmerHashIt = text | views::kmer_hash(shape) | ranges::to<std::vector<uint64_t>>();
	std::vector<uint64_t> revcHashIt = revComp | views::kmer_hash(shape) | ranges::to<std::vector<uint64_t>>();

	// stores the minimizer for the current window
	dna4_vector curMinimizer
	{ };
	// stores minimum kmer hash value in the current window
	uint64_t minHashValue = std::numeric_limits<uint64_t>::max();

	// Initialisation. We need to compute all hashes for the first window.
	for (uint32_t i = 0; i < windowKmers; ++i)
	{
		// Get smallest canonical k-mer (for reverse complement of kmer as well)
		uint64_t kmerHash = kmerHashIt[i] ^ seed;
		uint64_t revcHash = revcHashIt[kmerSize - i - 1] ^ seed;

		// get current kmer
		dna4_vector curKmer = text | views::slice(i, i + k) | ranges::to<std::vector<seqan3::dna4>>();
		dna4_vector curRevComKmer = revComp | views::slice(revComp.size() - i - k, revComp.size() - i) | ranges::to<std::vector<seqan3::dna4>>();

		// if hash value of current kmer is minimum hash value in the window -> set as minimizer of the window
		if (kmerHash < minHashValue)
		{
			minHashValue = kmerHash;
			curMinimizer = curKmer;
		}
		// hash value of reverse complement smaller? -> corresponding kmer is minimizer
		if (revcHash < minHashValue)
		{
			minHashValue = revcHash;
			curMinimizer = curRevComKmer;
		}

		// store sequence and hash value for current kmer for next window computation
		if (kmerHash <= revcHash)
		{
			windowValues.push_back(kmerHash);
			windowValueKmers.push_back(curKmer);
		}
		// use reverse complement if its hash value is smaller
		else
		{
			windowValues.push_back(revcHash);
			windowValueKmers.push_back(curRevComKmer);
		}
	}

	// store minimizer sequence of first window
	minimizers.push_back(curMinimizer);

	// For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
	// that results from the window shifting
	bool minimizer_changed
	{ false };
	for (uint64_t i = 1; i < possible; ++i)
	{
		// was first kmer in last window the minimizer?
		if (minHashValue == *std::begin(windowValues))
		{
			// delete first kmer from the containers since it's not part of the current window
			windowValues.pop_front();
			windowValueKmers.pop_front();
			// compute minimum hash value and kmer sequence as current minimizer for all but the last kmer in window
			minHashValue = std::numeric_limits<uint64_t>::max();
			for (unsigned int j = 0; j < windowValues.size(); ++j)
			{
				if (windowValues[j] < minHashValue)
				{
					minHashValue = windowValues[j];
					curMinimizer = windowValueKmers[j];
				}
			}
			minimizer_changed = true;
		}
		// just delete first kmer from the containers since it's not part of the current window
		else
		{
			windowValues.pop_front();
			windowValueKmers.pop_front();
		}

		// get hash values for latest kmer in current window
		uint64_t kmerHash = kmerHashIt[windowKmers + i - 1] ^ seed;
		uint64_t revcHash = revcHashIt[kmerSize - windowKmers - i] ^ seed;

		// get sequence of latest kmer in current window
		dna4_vector curKmer = text | views::slice(windowKmers + i - 1, windowKmers + i + k - 1) | ranges::to<std::vector<seqan3::dna4>>();
		dna4_vector curRevComKmer = revComp | views::slice(revComp.size() - windowKmers - i - k + 1, revComp.size() - windowKmers - i + 1)
		                            | ranges::to<std::vector<seqan3::dna4>>();
		auto test2 = curRevComKmer | std::views::reverse | views::complement;

		// store latest kmer hash value and minimizer sequence
		if (kmerHash <= revcHash)
		{
			windowValues.push_back(kmerHash);
			windowValueKmers.push_back(curKmer);
		}
		else
		{
			windowValues.push_back(revcHash);
			windowValueKmers.push_back(curRevComKmer);
		}

		// if latest kmer hash value is minmum in current window -> use as minimizer
		if (windowValues.back() < minHashValue)
		{
			minHashValue = windowValues.back();
			curMinimizer = windowValueKmers.back();
			minimizer_changed = true;
		}

		// if current window does not have the same minimizer as the window before -> add minimizer
		if (minimizer_changed)
		{
			minimizers.push_back(curMinimizer);
			minimizer_changed = false;
		}
	}

	return minimizers;
}

uint64_t hashFunction(dna4_vector& seq, uint64_t& seed)
{
	char data[16];
	std::vector<char> c = seq | seqan3::views::to_char | ranges::to<std::vector<char>>();
    MurmurHash3_x64_128(reinterpret_cast<char*>(c.data()), seq.size(), seed, data);
	return *((uint64_t *)data);
}

std::vector<uint64_t> Minimizer::getMinimizerHashValues(dna4_vector & text)
{
	//seqan3::debug_stream << text << std::endl;
	// no minimizers if text to short
	if (k > text.size())
		return std::vector<uint64_t>
		{ };

	// Reverse complement without copying/modifying the original string
	dna4_vector revComp = text | std::views::reverse | views::complement | ranges::to<std::vector<seqan3::dna4>>();;

	// computation of some parameters
	uint64_t possible = text.size() > w ? text.size() - w + 1 : 1; 	// number of all possible windows in text
	uint32_t windowKmers = w - k + 1;								// number of kmers in a window
	                                                                // corresponds to w in minimizer paper
	uint64_t kmerSize = text.size() - k + 1;						// number of all possible k-mers in text

	// resulting container with minimizer sequences that will be returned
	std::vector<uint64_t> minimizers
	{ };
	// reserve space for maximum number of possible minimizer
	minimizers.reserve(possible); 									// maybe rather reserve to expected?

	// Stores hash values for all k-mers in the current window
	std::deque<uint64_t> windowValues;
	// stores all kmers corresponding to stored hash values in the current window
	std::deque<dna4_vector> windowValueKmers;

	// ungapped shape or masked shape used for kmer hashing?
	if (shape.size() != k)
	{
		shape = computeSpacedSeed();
	}
	std::vector<uint64_t> kmerHashIt = text | views::kmer_hash(shape) | ranges::to<std::vector<uint64_t>>();
	std::vector<uint64_t> revcHashIt = revComp | views::kmer_hash(shape) | ranges::to<std::vector<uint64_t>>();

	// stores the minimizer for the current window
	dna4_vector curMinimizer
	{ };
	// stores minimum kmer hash value in the current window
	uint64_t minHashValue = std::numeric_limits<uint64_t>::max();

	// Initialisation. We need to compute all hashes for the first window.
	for (uint32_t i = 0; i < windowKmers; ++i)
	{
		// Get smallest canonical k-mer (for reverse complement of kmer as well)
		uint64_t kmerHash = kmerHashIt[i] ^ seed;
		uint64_t revcHash = revcHashIt[kmerSize - i - 1] ^ seed;

		// get current kmer
		dna4_vector curKmer = text | views::slice(i, i + k) | ranges::to<std::vector<seqan3::dna4>>();
		dna4_vector curRevComKmer = revComp | views::slice(revComp.size() - i - k, revComp.size() - i) | ranges::to<std::vector<seqan3::dna4>>();

		//uint64_t kmerHash = hashFunction(curKmer, seed);
		//uint64_t revcHash = hashFunction(curRevComKmer, seed);

		// if hash value of current kmer is minimum hash value in the window -> set as minimizer of the window
		if (kmerHash < minHashValue)
		{
			minHashValue = kmerHash;
			curMinimizer = curKmer;
		}
		// hash value of reverse complement smaller? -> corresponding kmer is minimizer
		if (revcHash < minHashValue)
		{
			minHashValue = revcHash;
			curMinimizer = curRevComKmer;
		}

		// store sequence and hash value for current kmer for next window computation
		if (kmerHash <= revcHash)
		{
			windowValues.push_back(kmerHash);
			windowValueKmers.push_back(curKmer);
		}
		// use reverse complement if its hash value is smaller
		else
		{
			windowValues.push_back(revcHash);
			windowValueKmers.push_back(curRevComKmer);
		}
	}

	// store minimizer sequence of first window
	auto mini = std::min_element(std::begin(windowValues), std::end(windowValues));
	minimizers.push_back(*mini);
	
	dna4_vector window = text | views::slice(0, w) | ranges::to<std::vector<seqan3::dna4>>();
	//seqan3::debug_stream << window << " " << curMinimizer << " " << *mini << std::endl;

	// For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
	// that results from the window shifting
	bool minimizer_changed
	{ false };
	for (uint64_t i = 1; i < possible; ++i)
	{
		// was first kmer in last window the minimizer?
		if (minHashValue == *std::begin(windowValues))
		{
			// delete first kmer from the containers since it's not part of the current window
			windowValues.pop_front();
			windowValueKmers.pop_front();
			// compute minimum hash value and kmer sequence as current minimizer for all but the last kmer in window
			minHashValue = std::numeric_limits<uint64_t>::max();
			for (unsigned int j = 0; j < windowValues.size(); ++j)
			{
				if (windowValues[j] < minHashValue)
				{
					minHashValue = windowValues[j];
					curMinimizer = windowValueKmers[j];
				}
			}
			minimizer_changed = true;
		}
		// just delete first kmer from the containers since it's not part of the current window
		else
		{
			windowValues.pop_front();
			windowValueKmers.pop_front();
		}

		// get hash values for latest kmer in current window
		uint64_t kmerHash = kmerHashIt[windowKmers + i - 1] ^ seed;
		uint64_t revcHash = revcHashIt[kmerSize - windowKmers - i] ^ seed;

		// get sequence of latest kmer in current window
		dna4_vector curKmer = text | views::slice(windowKmers + i - 1, windowKmers + i + k - 1) | ranges::to<std::vector<seqan3::dna4>>();
		dna4_vector curRevComKmer = revComp | views::slice(revComp.size() - windowKmers - i - k + 1, revComp.size() - windowKmers - i + 1)
		                            | ranges::to<std::vector<seqan3::dna4>>();
		auto test2 = curRevComKmer | std::views::reverse | views::complement;
		
		//uint64_t kmerHash = hashFunction(curKmer, seed);
		//uint64_t revcHash = hashFunction(curRevComKmer, seed);
		

		// store latest kmer hash value and minimizer sequence
		if (kmerHash <= revcHash)
		{
			windowValues.push_back(kmerHash);
			windowValueKmers.push_back(curKmer);
		}
		else
		{
			windowValues.push_back(revcHash);
			windowValueKmers.push_back(curRevComKmer);
		}

		// if latest kmer hash value is minmum in current window -> use as minimizer
		if (windowValues.back() < minHashValue)
		{
			minHashValue = windowValues.back();
			curMinimizer = windowValueKmers.back();
			minimizer_changed = true;
		}
		dna4_vector win = text | views::slice(windowKmers + i - 1, windowKmers + i + w - 1) | ranges::to<std::vector<seqan3::dna4>>();
		
		// if current window does not have the same minimizer as the window before -> add minimizer
		if (minimizer_changed)
		{
			auto mini = std::min_element(std::begin(windowValues), std::end(windowValues));
			minimizers.push_back(*mini);
			minimizer_changed = false;
		}
	//	seqan3::debug_stream << win << " " << curMinimizer << " " << minimizers.back() << std::endl;
	}

	return minimizers;
}