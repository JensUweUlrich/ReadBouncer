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


	public:

		uint16_t kMerSize;

		// Constructor
		CustomBloomFilter();
		CustomBloomFilter(const bloom_parameters& parameters, const uint16_t& newKmerSize);

		// Destructor
		~CustomBloomFilter();

		inline uint16_t getKmerSize()
		{
			return kMerSize;
		}

		
		inline std::vector<bloom_type> getHashSeeds()
		{
			std::vector<bloom_type> f(salt_);
			return f;
		}
		
		inline table_type getBitTable()
		{
			table_type b(bit_table_);
			return b;
		}

		void writeToFile(const std::experimental::filesystem::path& file);

		bool readFromFile(const std::experimental::filesystem::path& file);

};
