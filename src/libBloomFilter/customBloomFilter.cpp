#include "customBloomFilter.hpp"
#include "bloomFilterException.hpp"


// default constructor
CustomBloomFilter::CustomBloomFilter()
{
	// do nothing
}

/**
 * Bloom Filter constructor
 * initializing bloom filter with optimal parameters and hash function number
 *
 * @param parameters : parameters used for creating the bloom filter
 * @param k : k-mer size
 */
CustomBloomFilter::CustomBloomFilter(const bloom_parameters& parameters, const uint16_t& newKmerSize) : bloom_filter(parameters)
{
	// just call base class constructor and set k-mer size
	kMerSize = newKmerSize;
}

CustomBloomFilter::~CustomBloomFilter()
{
	// do nothing
}

/**
 * stores the bloom filter and all necessary information for reconstruction
 * in a binary file with the following setup:
 * 1. k-mer size
 * 2. number of hash seeds used
 * 3. hash seeds
 * 4. bloom filter size (number of bits)
 * 5. bloom filter as bytes
 *
 * @param file : filesystem path to the output file
 *
 */
void CustomBloomFilter::writeToFile(const std::filesystem::path& file)
{
	::std::ofstream fout;
	fout.exceptions(std::ofstream::failbit | std::ofstream::badbit);

	try
	{
		fout.open(file, std::ofstream::out);
	}
	catch (std::ofstream::failure &e)
	{
		throw BloomFilterException("Could not open " + file.u8string());
	}

	// store KmerSize in file
	try
	{
		fout.write((const char*) &kMerSize, sizeof(uint16_t));
	}
	catch (std::ofstream::failure &e)
	{
		throw BloomFilterException("Could not write Kmer size to output file!");
	}


	// store computed hash seeds
	try
	{
		// first get number of hash seeds to store to file
		std::size_t s = hash_count();
		fout.write((const char*) &s, sizeof(s));

		// then convert hash seeds to char and save in file
		for (bloom_type v : salt_)
		{
			std::cout << v << std::endl;
			fout.write((const char*) &v, sizeof(bloom_type));
		}
	}
	catch (std::ofstream::failure &e)
	{
		throw BloomFilterException("Could not write hash function values to output file!");
	}

	// first get number of bits in bloom filter and write it to file
	try
	{
		unsigned long long int n = size();
		fout.write((const char*) &n, sizeof(unsigned long long int));
		
		// then convert every bit to a char and write it to file
		int i = 0;
		for (unsigned char aggr : bit_table_)
		{
			if (i < 10 && int(aggr) != 0)
			{
				//std::cout << int(aggr) << std::endl;
				++i;
			}
			fout.write((const char*) &aggr, sizeof(unsigned char));
		}
	}
	catch (std::ofstream::failure &e)
	{
		throw BloomFilterException("Could not write bloom filter bit vector to output file!");
	}

	fout.close();

}


/**
 * parse bloom filter from an input binary file
 *
 * @param file : filesystem path to bloom filter input file
 *
 */
bool CustomBloomFilter::readFromFile(const std::filesystem::path& file)
{
	std::ifstream fin;
	fin.exceptions(std::ofstream::failbit | std::ofstream::badbit);

	try
	{
		fin.open(file, std::ifstream::in);
	}
	catch(std::ifstream::failure &e)
	{
		throw BloomFilterException("Could not open " + file.u8string());
	}

	// read K-mer size used for minimizer computation
	try
	{
		fin.read((char*) &kMerSize, sizeof(uint16_t));
	}
	catch (std::ifstream::failure &e)
	{
		throw BloomFilterException("Could not read Kmer size from "  + file.u8string());
	}

	std::cout << "K-mer size set to : " << kMerSize << std::endl;
	// read number of hash functions and resize array to store hash functions
	try
	{
		std::size_t m;
		fin.read((char*) &m, sizeof(size_t));
		std::cout << "Number of hash functions set to : " << m << std::endl;
		salt_count_ = m;
		salt_.resize(salt_count_);

		// read hash functions from file
		for (int i = 0; i < salt_count_; ++i)
		{
			fin.read((char*) &salt_.at(i), sizeof(bloom_type));
		}
		
		/*for (bloom_type seed : salt_)
		{
			std::cout << seed << std::endl;
		}*/
	}
	catch (std::ifstream::failure &e)
	{
		throw BloomFilterException("Could not read hash function values from " + file.u8string());
	}
	
	// read bloom filter from file
	try
	{
		unsigned long long int n;
		fin.read((char*) &n, sizeof(unsigned long long int));
		std::cout << "Bloom filter size : " << n << std::endl;
		table_size_ = n;
		bit_table_.resize(table_size_ / bits_per_char, static_cast<unsigned char>(0x00));
		
		int k = 0;
		for (unsigned long long int i = 0; i < bit_table_.size(); ++i)
		{
			fin.read((char*) &bit_table_.at(i), sizeof(unsigned char));
			
			if (k < 10 && int(bit_table_.at(i)) != 0)
			{
				//std::cout << int(bit_table_.at(i)) << std::endl;
				++k;
			}
		}
	}
	catch (std::ifstream::failure &e)
	{
		throw BloomFilterException("Could not read bit vector from " + file.u8string());
	}
	fin.close();

	return true;
}


