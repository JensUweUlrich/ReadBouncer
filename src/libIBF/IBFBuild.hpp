#pragma once

#include "Config.hpp"
#include "SafeQueue.hpp"
#include "StopClock.hpp"

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

    struct Seqs
    {
        std::string       seqid;
        seqan::Dna5String seq;
    };

    struct Stats
    {
        Stats()
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
        Tfilter;

    class IBF
    {

        private:
            Tfilter filter{};

            // parse reference sequences
            std::future< void > parse_ref_seqs(SafeQueue< Seqs > &queue_refs, std::mutex &mtx, interleave::Config &config, Stats &stats);
            void parse_seqid_bin( const std::string& seqid_bin_file, TSeqBin& seq_bin, std::set< uint64_t >& bin_ids );
            std::vector< std::future< void > > add_sequences_to_filter(Config &config, uint64_t &binid, SafeQueue< Seqs > &queue_refs);

        public:
            IBF(){}
			~IBF(){}
            void load_filter( Config& config, const std::set< uint64_t >& bin_ids, Stats& stats );
            bool build( Config config );
            void store_filter();
        

    };


} // namespace GanonBuild

#endif /* INTERLEAVE_IBFBUILD_HPP_ */
