#include <vector>
#include <cstdint>
#include <experimental/filesystem>
#include <fstream>
#include <vector>
#include <ctime>
#include <cstdint>
#include <cstdlib>
#include <math.h>

#include <seqan3/core/debug_stream.hpp>

// Bloom Filter class
class BloomFilter
{
	// class variables
	private:
		std::vector<uint16_t> m_hashes;
		std::vector<bool> m_bits;
		uint16_t kMerSize;

		void initialize_bloom_filter(const float p, const uint64_t minimizer_number);
		void initialize_hash_functions(const uint64_t func_number);
		uint64_t calculate_bloom_filter_size(const float p, const uint64_t minimizer_number);
		uint64_t calculate_hash_function_number(const uint64_t filter_size, uint64_t minimizer_number);


	// Constructor
	BloomFilter();
	BloomFilter(const uint64_t size, const uint8_t numHashes);
	BloomFilter(const float error_rate, const uint64_t minimizer_number, const uint16_t newKmerSize);

	public:
		// add new  to bloom filter
		void addHashValue(const uint64_t value);

		// check if the data is contained in the bloom filter
		bool possiblyContains(const uint64_t value) const;

		void writeToFile(const std::experimental::filesystem::path file);

		void readFromFile(const std::experimental::filesystem::path file);

};
