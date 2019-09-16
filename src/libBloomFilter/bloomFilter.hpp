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

		virtual void initialize_bloom_filter(const float p, const uint64_t minimizer_number);
		virtual void initialize_hash_functions(const uint64_t func_number);

		virtual inline uint64_t calculate_bloom_filter_size(const float p, const uint64_t minimizer_number)
		{
			return ceil(minimizer_number * (log(1/p)) / (pow(log(2.0), 2)));
		}

		virtual inline uint64_t calculate_hash_function_number(const uint64_t filter_size, uint64_t minimizer_number)
		{
			return ceil(filter_size * log(2.0)/ minimizer_number);
		}

	// Constructor
	BloomFilter();
	BloomFilter(const uint64_t size, const uint8_t numHashes);
	BloomFilter(const float error_rate, const uint64_t minimizer_number, const uint16_t newKmerSize);

	// Destructor
	~BloomFilter();

	public:

		std::vector<uint16_t> hashes;
		std::vector<bool> bits;
		uint16_t kMerSize;

		/**
		 * adds a k-mer to the Bloom Filter
		 * @param value : hash value of the k-mer
		 */
		virtual inline void addHashValue(const uint64_t value)
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
		virtual inline bool possiblyContains(const uint64_t value)
		{
			bool result = bits[value % bits.size()];
			for (uint16_t f : hashes)
			{
				result &= bits[(value ^ f) % bits.size()];
			}
			return result;

		}

		virtual void writeToFile(const std::experimental::filesystem::path file);

		virtual bool readFromFile(const std::experimental::filesystem::path file);

};
