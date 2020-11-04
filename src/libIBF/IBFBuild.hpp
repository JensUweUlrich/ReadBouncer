#pragma once

#include "Config.hpp"
#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>

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

namespace ibf
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
            Tfilter filter;

            void parse_seqid_bin( const std::string& seqid_bin_file, TSeqBin& seq_bin, std::set< uint64_t >& bin_ids );

        public:
            void load_filter( Config& config, const std::set< uint64_t >& bin_ids, Stats& stats );
            bool build( Config config );
            void store_filer();
        

    }


} // namespace GanonBuild
