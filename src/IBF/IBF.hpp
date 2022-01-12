#pragma once

#include "IBFConfig.hpp"
#include "SafeQueue.hpp"
#include "StopClock.hpp"
#include "IBFExceptions.hpp"
#include <ont_read.hpp>

// seqan libraries
#include <seqan/binning_directory.h>

// spdlog
#include "spdlog/spdlog.h"
#include "spdlog/sinks/rotating_file_sink.h"

#include <cinttypes>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace interfaces;

#ifndef INTERLEAVE_IBFBUILD_HPP_
#define INTERLEAVE_IBFBUILD_HPP_

namespace interleave
{

    struct ReadMatch
    {
        //std::string target;
        uint32_t    binNo;
        uint16_t    kmer_count;
    };

    struct Seqs
    {
        std::string       seqid;
        seqan::Dna5String seq;
    };

    struct FilterStats
    {
        FilterStats()
        : sumSeqLen{ 0 }
        , totalSeqsBinId{ 0 }
        , totalBinsBinId{ 0 }
        , totalSeqsFile{ 0 }
        , totalBinsFile{ 0 }
        , invalidSeqs{ 0 }
        , newBins{ 0 }
        {
        }
        

        uint64_t sumSeqLen;
        uint64_t totalSeqsBinId;
        uint32_t totalBinsBinId;
        uint64_t totalSeqsFile;
        uint32_t totalBinsFile;
        uint64_t invalidSeqs;
        uint32_t newBins;

        StopClock timeLoadSeq;
        StopClock timeLoadFilter;
        StopClock timeSaveFilter;
        StopClock timeBuild;
        StopClock timeLoadFiles;
        StopClock timeIBF;
    };

    struct FragmentBin
    {
        uint64_t start;
        uint64_t end;
        uint64_t bin_id;
    };

     

    typedef std::map< std::string, std::vector< FragmentBin > > TSeqBin;

    typedef seqan::BinningDirectory< seqan::InterleavedBloomFilter,
                                    seqan::BDConfig< seqan::Dna5, seqan::Normal, seqan::Uncompressed > >
        TIbf;

    typedef seqan::ModifiedString< seqan::ModifiedString< seqan::Dna5String, seqan::ModComplementDna >, seqan::ModReverse >
        TSeqRevComp;

    typedef std::unordered_map< uint32_t, uint16_t > TMatches;

    typedef std::pair<uint16_t, uint16_t> TInterval;

    class IBF
    {

        private:
            TIbf filter{};

            // parse reference sequences
            std::future< void > parse_ref_seqs(SafeQueue< Seqs > &queue_refs, interleave::IBFConfig &config, FilterStats &stats);
            void add_sequences_to_filter(std::vector< std::future< void > >& tasks, IBFConfig &config, uint64_t &binid, SafeQueue< Seqs > &queue_refs);
            std::vector<std::string> cutOutNNNs(std::string& seq, uint64_t seqlen);
            uint64_t calculate_filter_size_bits(IBFConfig& config, uint64_t numberOfBins);
            std::shared_ptr<spdlog::logger> ibf_logger;

        public:
            IBF()
            {
                if ((ibf_logger = spdlog::get("IbfLog")) == nullptr)
                {
                    try
                    {
                        ibf_logger = spdlog::rotating_logger_mt("IbfLog", InterleavedBloomFilterLog + "logs/InterleavedBloomFilterLog.txt", 1048576 * 5, 100);
                    }
                    catch (const spdlog::spdlog_ex& e)
                    {
                        std::cerr << "IBF Log initialization failed: " << e.what() << std::endl;
                    }
                    ibf_logger->set_level(spdlog::level::debug);
                    //ibf_logger->flush_on(spdlog::level::debug);
                }
            }
			~IBF(){}
            FilterStats create_filter(IBFConfig& config);
            FilterStats load_filter( IBFConfig& config);
            FilterStats update_filter(IBFConfig& config);
            inline TIbf getFilter()
            {
                return filter;
            }
        

    };

     struct IBFMeta
    {
        TIbf filter;
        std::string name;
        seqan::SeqFileOut outfile;
        uint64_t classified;
    };

    class Read
    {
        private:
            
            std::vector< ReadMatch > matches;
            uint16_t max_kmer_count = 0;
          
            bool find_matches(std::vector< TIbf >& filters, ClassifyConfig&  config );
            uint64_t count_matches(IBFMeta& filter, ClassifyConfig& config);
            bool select_matches(std::vector< uint16_t >& selectedBins,
                                std::vector< uint16_t >& selectedBinsRev,
                                TIbf&                    filter,
                                uint16_t                 threshold);
            uint64_t max_matches(std::vector< uint16_t >& selectedBins, std::vector< uint16_t >& selectedBinsRev,
                TIbf& filter, uint16_t threshold);

        public:
            seqan::Dna5String sequence{};
            std::string id{};

            Read() {}

            Read( std::string& id, seqan::Dna5String& seq)
            {
                sequence = seq;
                this->id = id;
            }

            Read(seqan::CharString& id, seqan::Dna5String& seq)
            {
                sequence = seq;
                this->id = toCString(id);
            }
            
            ~Read(){}

            inline uint32_t getReadLength()
            {
                return seqan::length(this->sequence);
            }

            bool classify(std::vector< TIbf >& filters, ClassifyConfig& config);
            int classify(std::vector< IBFMeta >& filters, ClassifyConfig& config);
            std::pair<int, int> classify(std::vector< IBFMeta >& filt1, std::vector< IBFMeta >& filt2, ClassifyConfig& config);
        
    };
    
    typedef std::vector<Read> TReads;

     /**
         Return threshold (number of kmers) based on an percentage of kmers. 0 for anything with at least 1 k-mer
         @DEPRECATED
     */
    inline uint16_t get_threshold_kmers( uint16_t readLen, uint16_t kmerSize, float min_kmers )
    {
        // ceil -> round-up min # k-mers

        return min_kmers > 0 ? std::ceil( ( readLen - kmerSize + 1) * min_kmers ) : 1u;
    }

     /**
         Return the optimal number of errors for a certain sequence based on the kmer_count
         @DEPRECATED
    */
    inline uint16_t get_error( uint16_t readLen, uint16_t kmerSize, uint16_t kmer_count)
    {
        return std::ceil( ( readLen - kmerSize - kmer_count + 1 )
                            / static_cast< float >( kmerSize ) );
    }

    /**
        Return threshold (number of kmers) based on an optimal number of errors
        minimum number of matching kmers of a read given max_error are allowed
        1 instead of 0 - meaning that if a higher number of errors are allowed the threshold here is
        just one kmer match (0 would match every read everywhere)
        @DEPRECATED
    */
    inline uint16_t get_threshold_errors( uint16_t readLen, uint16_t kmerSize, uint16_t max_error )
    {
        return readLen + 1u > kmerSize * ( 1u + max_error )
                ? readLen - kmerSize - max_error * kmerSize  + 1 
                : 1u;
    }

    /**
    *   Abramowitz-Stegun-Approximation for the inverse normal CDF
    */
    inline double RationalApproximation(double t)
    {
        // Abramowitz and Stegun formula 26.2.23.
        // The absolute value of the error should be less than 4.5 e-4.
        double c[] = { 2.515517, 0.802853, 0.010328 };
        double d[] = { 1.432788, 0.189269, 0.001308 };
        return t - ((c[2] * t + c[1]) * t + c[0]) /
            (((d[2] * t + d[1]) * t + d[0]) * t + 1.0);
    }

    /**
    *   approximates the value of the inverse normal cumulative distribution function
    *   @p      : probability, has to be between 0 and 1
    *   @return : z score
    */
    inline double NormalCDFInverse(double p)
    {
        if (p <= 0.0 || p >= 1.0)
        {
            std::stringstream os;
            os << "Invalid input argument (" << p
                << "); must be larger than 0 but less than 1.";
            throw std::invalid_argument(os.str());
        }

        // See article above for explanation of this section.
        if (p < 0.5)
        {
            // F^-1(p) = - G^-1(p)
            return -RationalApproximation(sqrt(-2.0 * log(p)));
        }
        else
        {
            // F^-1(p) = G^-1(1-p)
            return RationalApproximation(sqrt(-2.0 * log(1.0 - p)));
        }
    }

    /**
    *   calculate die confidence interval for the number of errorneous kmers based on read length, kmer size and error rate
    *   based on "statistics of kmers from a sequence undergoing a simple mutation process without spurious matches" 
    *   by Blanca, A., Harris, R., Koslicki, D. and Medvedev, P.
    *   @r          : assumed sequencing error rate 
    *   @kmer_size  : size of kmers used
    *   @readlen    : length of the current read
    *   @confidence : significance level, e.g. 0.95 for a 95% confidence interval
    *   @return     : pair of integer values representing the boundaries of the confidence interval 
    */
    inline TInterval calculateCI(const double r, const uint8_t kmer_size, const uint32_t readlen, const double confidence)
    {
        double q = 1.0 - pow(1.0 - r, kmer_size);
        // number of kmers in sequence of length readlen
        double L = ((double)readlen - (double)kmer_size + 1.0);
        // expected number of errorneous/mutated kmers in sequence of length readlen
        double Nmut = L * q;
        // compute variance
        double varN = L * (1.0 - q) * (q * (2.0 * (double)kmer_size + (2.0 / r) - 1.0) - 2.0 * (double)kmer_size)
                        + (double)kmer_size * ((double)kmer_size - 1.0) * pow((1.0 - q), 2.0)
                        + (2.0 * (1.0 - q) / (pow(r, 2.0))) * ((1.0 + ((double)kmer_size - 1.0) * (1.0 - q)) * r - q);
        double alpha = 1 - confidence;
        
        double z = NormalCDFInverse(1.0 - alpha / 2.0);
        uint16_t low = (uint16_t) floor(L * q - z * sqrt(varN));
        uint16_t high = (uint16_t)ceil(L * q + z * sqrt(varN));
        TInterval ci{ low , high };
        return ci;
    }

    void print_build_stats( interleave::FilterStats& stats);
    void print_load_stats(interleave::FilterStats& stats);
    void print_time( const interleave::IBFConfig& config, interleave::FilterStats& stats );

} // namespace interleave

#endif /* INTERLEAVE_IBFBUILD_HPP_ */
