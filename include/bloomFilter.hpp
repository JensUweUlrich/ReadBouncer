#include <vector>
#include <cstdint>

// Bloom Filter class
class BloomFilter
{
	// class variables
	private:
		std::vector<uint16_t> m_hashes;
		std::vector<bool> m_bits;


	// Constructor
	BloomFilter(const uint64_t size, const uint8_t numHashes);

	public:
		// add new  to bloom filter
		void addHashValue(const uint64_t value);

		// check if the data is contained in the bloom filter
		bool possiblyContains(const uint64_t value) const;


};
