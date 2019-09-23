#include <vector>
#include <cstdint>
#include <experimental/filesystem>
#include <fstream>
#include <vector>
#include <ctime>
#include <cstdint>
#include <cstdlib>
#include <math.h>

#include "bloom_filter.hpp"

#include <seqan3/core/debug_stream.hpp>

// Bloom Filter class
class CustomBloomFilter : public bloom_filter
{
	// class variables
	private:
		uint16_t kMerSize;
		

	// Constructor
	CustomBloomFilter();
	CustomBloomFilter(const bloom_parameters& parameters, const uint16_t& newKmerSize);

	// Destructor
	~CustomBloomFilter();

	public:

		virtual inline uint16_t getKmerSize()
		{
			return kMerSize;
		}
		
		virtual inline std::vector<bloom_type> getHashSeeds()
		{
			std::vector<bloom_type> f(salt_);
			return f;
		}
		
		virtual inline table_type getBitTable()
		{
			table_type b(bit_table_);
			return b;
		}

		virtual void writeToFile(const std::experimental::filesystem::path& file);

		virtual bool readFromFile(const std::experimental::filesystem::path& file);

};
