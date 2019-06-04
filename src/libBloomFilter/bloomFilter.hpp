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

		void initialize_bloom_filter(const float p, const uint64_t minimizer_number);
		void initialize_hash_functions(const uint64_t func_number);

		inline uint64_t calculate_bloom_filter_size(const float p, const uint64_t minimizer_number)
		{
			return ceil((-1) * (minimizer_number * log(p)) / (pow(log(2.0), 2)));
		}

		inline uint64_t calculate_hash_function_number(const uint64_t filter_size, uint64_t minimizer_number)
		{
			return ceil(filter_size / minimizer_number * log(2.0));
		}

	// Constructor
	BloomFilter();
	BloomFilter(const uint64_t size, const uint8_t numHashes);
	BloomFilter(const float error_rate, const uint64_t minimizer_number, const uint16_t newKmerSize);

	public:

		std::vector<uint16_t> hashes;
		std::vector<bool> bits;
		uint16_t kMerSize;

		/**
		 * adds a k-mer to the Bloom Filter
		 * @param value : hash value of the k-mer
		 */
		inline void addHashValue(const uint64_t value)
		{
			bits[value % bits.size()] = true;
			for (uint16_t f : hashes)
			{
				bits[(value ^ f) % bits.size()] = true;
			}
		}

		/**
		 * check if a k-mer is contained in the Bloom Filter
		 * @param value : hash value of the k-mer
		 */
		inline bool possiblyContains(const uint64_t value) const
		{
			bool result = bits[value % bits.size()];
			for (uint16_t f : hashes)
			{
				result &= bits[(value ^ f) % bits.size()];
			}
			return result;

		}

		void writeToFile(const std::experimental::filesystem::path file);

		void readFromFile(const std::experimental::filesystem::path file);

};
