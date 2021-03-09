#pragma once

#include "IBFConfig.hpp"
#include "SafeQueue.hpp"
#include "StopClock.hpp"
#include "IBFExceptions.hpp"

#include <seqan/binning_directory.h>

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

     struct DepleteConfig
    {
        float       min_kmers;
        uint16_t    max_error;
        uint16_t    strata_filter;
    };

    typedef std::map< std::string, std::vector< FragmentBin > > TSeqBin;

    typedef seqan::BinningDirectory< seqan::InterleavedBloomFilter,
                                    seqan::BDConfig< seqan::Dna5, seqan::Normal, seqan::Uncompressed > >
        TIbf;

    typedef seqan::ModifiedString< seqan::ModifiedString< seqan::Dna5String, seqan::ModComplementDna >, seqan::ModReverse >
        TSeqRevComp;

    typedef std::unordered_map< uint32_t, uint16_t > TMatches;

    class IBF
    {

        private:
            TIbf filter{};

            // parse reference sequences
            std::future< void > parse_ref_seqs(SafeQueue< Seqs > &queue_refs, interleave::IBFConfig &config, FilterStats &stats);
            void add_sequences_to_filter(std::vector< std::future< void > >& tasks, IBFConfig &config, uint64_t &binid, SafeQueue< Seqs > &queue_refs);
            std::vector<std::string> cutOutNNNs(std::string& seq, uint64_t seqlen);
            uint64_t calculate_filter_size_bits(IBFConfig& config, uint64_t numberOfBins);

        public:
            IBF(){}
			~IBF(){}
            FilterStats create_filter(IBFConfig& config);
            FilterStats load_filter( IBFConfig& config);
            FilterStats update_filter(IBFConfig& config);
            //std::vector<std::string> cutOutNNNs(std::string& seq);
            inline TIbf getFilter()
            {
                return filter;
            }
        

    };

    class Read
    {
        private:
            seqan::CharString id;
            seqan::Dna5String seq;
            std::vector< ReadMatch > matches;
            uint16_t max_kmer_count = 0;
            uint16_t channelNr = 0;
            uint16_t readNr = 0;
            TimeMeasures processingTimes{};

            uint32_t filter_matches(TMatches& matches,
                                    uint16_t  len,
                                    uint16_t  kmer_size,
                                    int16_t   strata_filter );
            void find_matches( TMatches& matches, 
                                    std::vector< TIbf >& filters, 
                                    DepleteConfig&  config );
            void select_matches(TMatches&                matches,
                                std::vector< uint16_t >& selectedBins,
                                std::vector< uint16_t >& selectedBinsRev,
                                TIbf&                    filter,
                                uint16_t                 threshold);

        public:
            Read(){}
            Read(  seqan::CharString id, seqan::Dna5String seq)
            {
                this->id    = id;
                this->seq   = seq;
            }
            Read(seqan::CharString id, seqan::Dna5String seq, 
                 uint16_t channelNr, uint16_t readNr, TimeMeasures processTimes)
            {
                this->id = id;
                this->seq = seq;
                this->channelNr = channelNr;
                this->readNr = readNr;
                this->processingTimes = processTimes;
            }
            ~Read(){}
            inline seqan::CharString getID()
            {
                return id;
            }
            inline std::vector< ReadMatch > getMatches()
            {
                return matches;
            }
            inline uint16_t getMaxKmerCount()
            {
                return max_kmer_count;
            }
            inline uint32_t getSeqLength()
            {
                return seqan::length(seq);
            }
            inline uint16_t getChannelNr()
            {
                return channelNr;
            }
            inline uint16_t getReadNr()
            {
                return readNr;
            }
            inline TimeMeasures getProcessingTimes()
            {
                return processingTimes;
            }
            bool classify(std::vector< TIbf >& filters, DepleteConfig& config);
        
        
    };
    
    typedef std::vector<Read> TReads;

     /**
         Return threshold (number of kmers) based on an percentage of kmers. 0 for anything with at least 1 k-mer
     */
    inline uint16_t get_threshold_kmers( uint16_t readLen, uint16_t kmerSize, float min_kmers )
    {
        // ceil -> round-up min # k-mers

        return min_kmers > 0 ? std::ceil( ( readLen - kmerSize + 1) * min_kmers ) : 1u;
    }

     /**
         Return the optimal number of errors for a certain sequence based on the kmer_count
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
    */
    inline uint16_t get_threshold_errors( uint16_t readLen, uint16_t kmerSize, uint16_t max_error )
    {
        return readLen + 1u > kmerSize * ( 1u + max_error )
                ? readLen - kmerSize - max_error * kmerSize  + 1 
                : 1u;
    }

    void print_stats( interleave::FilterStats& stats);
    void print_time( const interleave::IBFConfig& config, interleave::FilterStats& stats );

} // namespace interleave

#endif /* INTERLEAVE_IBFBUILD_HPP_ */
