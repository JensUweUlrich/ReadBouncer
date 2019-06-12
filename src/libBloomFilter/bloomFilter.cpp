#include "bloomFilter.hpp"

// default constructor
BloomFilter::BloomFilter()
{
	// do nothing
}

/**
 * Bloom Filter constructor
 * initializing bloom filter with optimal filter size and hash function number
 *
 * @param error_rate : allowed false positive rate for bloom filter size calculation
 * @param minimizer_number : number of minimizers to be stored in bloom filter
 * @param k : k-mer size
 */
BloomFilter::BloomFilter(const float error_rate, const uint64_t minimizer_number, const uint16_t k)
{
	kMerSize = k;
	initialize_bloom_filter(error_rate, minimizer_number);
}

/**
 * Bloom Filter Constructor
 * for fixed filter size and hash function number
 *
 * @param size : default bloom filter size
 * @param numHashes : number of hash functions
 */
BloomFilter::BloomFilter(uint64_t size, uint8_t numHashes)
{
	bits.resize(size);
	initialize_hash_functions(numHashes);
}

BloomFilter::~BloomFilter()
{
	// do nothing
}

/**
 * stores the bloom filter and all necessary information for reconstruction
 * in a binary file with the following setup:
 * 1. k-mer size
 * 2. number of hash functions/values used
 * 3. hash values
 * 4. bloom filter size (number of bits)
 * 5. bloom filter as bytes
 *
 * @param file : filesystem path to the output file
 *
 */
void BloomFilter::writeToFile(const std::experimental::filesystem::path file)
{
	std::ofstream fout(file, std::ofstream::out);

	// store KmerSize in file
	fout.write((const char*) &kMerSize, sizeof(uint16_t));

	// first get number of hash functions/values to store to file
	std::vector<uint16_t>::size_type m = hashes.size();
	fout.write((const char*) &m, sizeof(std::vector<bool>::size_type));

	// then convert hash functions to char and save in file
	for (uint16_t v : hashes)
	{
		fout.write((const char*) &v, sizeof(uint16_t));
	}

	// first get number of bits in bloom filter and write it to file
	std::vector<bool>::size_type n = bits.size();
	fout.write((const char*) &n, sizeof(std::vector<bool>::size_type));

	// then convert every bit to a char and write it to file
	for (std::vector<bool>::size_type i = 0; i < n;)
	{
		unsigned char aggr = 0;
		for (unsigned char mask = 1; mask > 0 && i < n; ++i, mask <<= 1)
			if (bits.at(i))
				aggr |= mask;
		fout.write((const char*) &aggr, sizeof(unsigned char));
	}

	fout.close();

}

/**
 * parse bloom filter from an input binary file
 *
 * @param file : filesystem path to bloom filter input file
 *
 */
bool BloomFilter::readFromFile(const std::experimental::filesystem::path file)
{
	std::ifstream fin(file, std::ifstream::in);

	// read K-mer size used for minimizer computation
	fin.read((char*) &kMerSize, sizeof(uint16_t));

	// read number of hash functions and resize array to store hash functions
	std::vector<uint16_t>::size_type m;
	fin.read((char*) &m, sizeof(std::vector<uint16_t>::size_type));
	hashes.resize(m);

	// read hash functions from file
	for (int i = 0; i < m; ++i)
	{
		fin.read((char*) &hashes.at(i), sizeof(uint16_t));
	}

	// read bloom filter bitwise from file
	std::vector<bool>::size_type n;
	fin.read((char*) &n, sizeof(std::vector<bool>::size_type));
	bits.resize(n);
	for (std::vector<bool>::size_type i = 0; i < n;)
	{
		unsigned char aggr;
		fin.read((char*) &aggr, sizeof(unsigned char));
		for (unsigned char mask = 1; mask > 0 && i < n; ++i, mask <<= 1)
			bits.at(i) = aggr & mask;
	}
	fin.close();

	return true;
}

/**
 * calculates the bloom filter size and the number of hash functions needed for
 * given false positive rate and the number of minimizers to be stored in the bloom filter;
 * resizes bit vector and hash value vector accordingly
 *
 * @param p : maximum false positive rate
 * @param minimizer number : number of minimizers to be stored in the bloom filter
 */
void BloomFilter::initialize_bloom_filter(const float p, const uint64_t minimizer_number)
{
	// initialize optimal bloom filter size
	uint64_t size = calculate_bloom_filter_size(p, minimizer_number);
	bits.resize(size);

	// initialize optimal number of hash functions
	uint64_t hash_func = calculate_hash_function_number(size, minimizer_number);
	initialize_hash_functions(hash_func);
}

/**
 * initialize a vector with random numbers used as hash function values
 *
 * @param func_number : number of hash function values to initialize
 */
void BloomFilter::initialize_hash_functions(const uint64_t func_number)
{
	srand(time(NULL));
	// initialize vector with random numbers to XOR given hash value
	for (uint64_t i = 0; i < func_number; ++i)
	{
		uint64_t randNr = (uint64_t) rand();
		hashes.push_back(randNr);
	}
}

