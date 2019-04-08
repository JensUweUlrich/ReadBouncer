#include "bloomFilter.hpp"

#include <vector>
#include <ctime>
#include <cstdint>
#include <cstdlib>


// Implementation of bloom filter constructor
BloomFilter::BloomFilter(uint64_t size, uint8_t numHashes)
{
	m_bits.resize(size);
	srand(time(NULL));
	// initialize vector with random numbers to XOR given hash value
	for (uint8_t i = 0; i < numHashes; ++i)
	{
		uint16_t randNr = (uint16_t) rand();
		m_hashes.push_back(randNr);
	}

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

