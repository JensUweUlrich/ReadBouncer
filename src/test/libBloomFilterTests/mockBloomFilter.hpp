#include "bloomFilter.hpp"
#include "gmock/gmock.h"

class MockBloomFilter : public BloomFilter
{

	public:

		MOCK_METHOD1(initialize_hash_functions, void(const uint64_t func_number));
		MOCK_METHOD1(addHashValue, void(const uint64_t value));
		MOCK_METHOD1(possiblyContains, bool(const uint64_t value));
		MOCK_METHOD1(writeToFile, void(const std::experimental::filesystem::path file));
		MOCK_METHOD1(readFromFile, bool(const std::experimental::filesystem::path file));

		MOCK_METHOD2(initialize_bloom_filter, void(float p, uint64_t minimizer_number));
		MOCK_METHOD2(calculate_hash_function_number, uint64_t(const uint64_t filter_size, uint64_t minimizer_number));
		MOCK_METHOD2(calculate_bloom_filter_size, uint64_t(const float p, const uint64_t minimizer_number));

};
