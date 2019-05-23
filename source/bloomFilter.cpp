#include "bloomFilter.hpp"


// default constructor
BloomFilter::BloomFilter()
{
	// do nothing
}

/**
 * Bloom Filter constructor
 * initializing bloom filter with optimal filter size and hash function number
 */
BloomFilter::BloomFilter(const float error_rate, const uint64_t minimizer_number)
{
	initialize_bloom_filter(error_rate, minimizer_number);
}

/**
 * Bloom Filter Constructor
 * for fixed filter size and hash function number
 */
BloomFilter::BloomFilter(uint64_t size, uint8_t numHashes)
{
	m_bits.resize(size);
	initialize_hash_functions(numHashes);
}

/**
 * adds a k-mer to the Bloom Filter
 * @param value : hash value of the k-mer
 */
void BloomFilter::addHashValue(const uint64_t value)
{
	m_bits[value % m_bits.size()] = true;
	for (uint16_t f : m_hashes)
	{
		m_bits[(value ^ f) % m_bits.size()] = true;
	}
}

/**
 * check if a k-mer is contained in the Bloom Filter
 * @param value : hash value of the k-mer
 */
bool BloomFilter::possiblyContains(const uint64_t value) const
{
	bool result = m_bits[value];
	for (uint16_t f : m_hashes)
	{
		result &= m_bits[value ^ f];
	}
	return result;

}

void BloomFilter::writeToFile(const std::experimental::filesystem::path file)
{
	std::ofstream fout(file, std::ofstream::out);

	// first get number of hash functions/values to store to file
	std::vector<uint16_t>::size_type m = m_hashes.size();
	fout.write((const char*) &m, sizeof(std::vector<bool>::size_type));

	// then convert hash functions to char and save in file
	for (uint16_t v : m_hashes)
	{
		fout.write((const char*) &v, sizeof(uint16_t));
	}

	// first get number of bits in bloom filter and write it to file
	std::vector<bool>::size_type n = m_bits.size();
	fout.write((const char*) &n, sizeof(std::vector<bool>::size_type));

	// then convert every bit to a char and write it to file
	for (std::vector<bool>::size_type i = 0; i < n;)
	{
		unsigned char aggr = 0;
		for (unsigned char mask = 1; mask > 0 && i < n; ++i, mask <<= 1)
			if (m_bits.at(i))
				aggr |= mask;
		fout.write((const char*) &aggr, sizeof(unsigned char));
	}

	fout.close();

}

void BloomFilter::readFromFile(const std::experimental::filesystem::path file)
{
	std::ifstream fin(file, std::ifstream::in);

	std::vector<uint16_t>::size_type m;
	fin.read((char*) &m, sizeof(std::vector<uint16_t>::size_type));
	m_hashes.resize(m);

	for (int i = 0; i < m; ++i)
	{
		fin.read((char*) &m_hashes.at(i), sizeof(uint16_t));
	}

	std::vector<bool>::size_type n;
	fin.read((char*) &n, sizeof(std::vector<bool>::size_type));
	m_bits.resize(n);
	for (std::vector<bool>::size_type i = 0; i < n;)
	{
		unsigned char aggr;
		fin.read((char*) &aggr, sizeof(unsigned char));
		for (unsigned char mask = 1; mask > 0 && i < n; ++i, mask <<= 1)
			m_bits.at(i) = aggr & mask;
	}
	fin.close();
}

void BloomFilter::initialize_bloom_filter(const float p, const uint64_t minimizer_number)
{
	// initialize optimal bloom filter size
	uint64_t size = calculate_bloom_filter_size(p, minimizer_number);
	m_bits.resize(size);

	seqan3::debug_stream << "Bloom Filter size: " << size << "\n";

	// initialize optimal number of hash functions
	uint64_t hash_func = calculate_hash_function_number(size, minimizer_number);
	initialize_hash_functions(hash_func);

	seqan3::debug_stream << "Number of Hash Functions: " << hash_func << "\n";
}

void BloomFilter::initialize_hash_functions(const uint64_t func_number)
{
	srand(time(NULL));
	// initialize vector with random numbers to XOR given hash value
	for (uint64_t i = 0; i < func_number; ++i)
	{
		uint64_t randNr = (uint64_t) rand();
		m_hashes.push_back(randNr);
	}
}

uint64_t BloomFilter::calculate_bloom_filter_size(const float p, const uint64_t minimizer_number)
{
	return ceil((-1) * (minimizer_number * log(p)) / (pow(log(2.0), 2)));
}

uint64_t BloomFilter::calculate_hash_function_number(const uint64_t filter_size, uint64_t minimizer_number)
{
	return ceil(filter_size / minimizer_number * log(2.0));
}

