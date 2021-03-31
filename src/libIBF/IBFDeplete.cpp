#include "IBF.hpp"

namespace interleave
{

    /**
        store all targets that share > threshold kmers with the read in matches map
        store maximum number of kmers that match against best matching target
        @matches:           unordered map to store matches (matches[binNo] = maxKmerCountBin)
        @selectedBins:      number of kmer matches for every bin in the filter
        @selectedBinsRev:   number of reverse complement kmer matches for every bin in the filter
        @filter:            actual IBF
        @threshold:         minimum number of kmer matches needed to add the bin number with its kmer match count to the matches map 
    */
    void Read::select_matches(TMatches&                matches,
                        std::vector< uint16_t >& selectedBins,
                        std::vector< uint16_t >& selectedBinsRev,
                        TIbf&                    filter,
                        uint16_t                 threshold )
    {
        // for each bin
        // for ( uint32_t binNo = 0; binNo < filter.ibf.noOfBins; ++binNo )
        // loop in map structure to avoid extra validations when map.size() < filter.ibf.noOfBins when ibf is updated and
        // sequences removed also avoid the error of spurius results from empty bins (bug reported)
       
        for ( uint32_t binNo = 0; binNo < filter.noOfBins; ++binNo ) 
        {
            //std::cout<<selectedBins[binNo] << " , " << selectedBinsRev[binNo] <<std::endl;
            // if kmer count is higher than threshold
            if ( selectedBins[binNo] >= threshold || selectedBinsRev[binNo] >= threshold )
            {
                // get best matching strand
                uint16_t maxKmerCountBin = std::max( selectedBins[binNo], selectedBinsRev[binNo] );
                
                // keep only the best match target/read when same targets are split in several
                // bins
                // only needed if we have different ibfs for different ref seqs
                //if ( matches.count( binNo ) == 0 || maxKmerCountBin > matches[target] )
                //{
                // store match to target
                    matches[binNo] = maxKmerCountBin;
                    if ( maxKmerCountBin > this->max_kmer_count )
                        this->max_kmer_count = maxKmerCountBin;
                //}
            }
        }
        //if (this->max_kmer_count > 0)
        //    std::cout << "threshold : " << threshold << ", Max Kmer Count: " << this->max_kmer_count <<std::endl;
    }

    /**
        find matches of kmers in ibfs
        report only targets where at least n > threshold kmers match target in the ibf
        match forward and reverse complement of kmers
        @matches:   unordered map to store matches 
        @filters:   vector of IBFs
        @config:    configuration settings needed
        @throws:    CountKmerException
    */
    void Read::find_matches(TMatches& matches, 
                                std::vector< TIbf >& filters, 
                                ClassifyConfig& config )
    {
        std::shared_ptr<spdlog::logger> logger = config.classification_logger;
        // iterate over all ibfs
        for ( TIbf& filter : filters )
        {
            try
            {
                // IBF count
                // find all bins that have at least 1 kmer in common with read_seq
                // selectedBins[binNo] = # of kmer matches for this bin
                std::vector< uint16_t > selectedBins    = seqan::count( filter, this->seq );
                std::vector< uint16_t > selectedBinsRev = seqan::count( filter, TSeqRevComp( this->seq ) );

                // get calculated threshold for minimum number of kmers needed to report a match
                // this is based on the confidence interval for mutated kmers in a read with expected error rate, kmer size
                TInterval ci = calculateCI(config.error_rate, filter.kmerSize, seqan::length(this->seq), config.significance);
                //std::cout << "CI = [" << ci.first << " , " << ci.second << "]" << std::endl;
                uint16_t readlen = seqan::length(this->seq);
                // minimum number of kmers = max number of kmer in read - upper bound of the CI
                //std::cout << seqan::length(this->seq) << " " << filter.kmerSize << std::endl;
                int16_t threshold = readlen - filter.kmerSize + 1 - ci.second;
                //uint16_t threshold = seqan::length(this->seq) - filter.kmerSize + 1 - (floor((ci.second - ci.first) / 2) + ci.first);
                // select matches above chosen threshold
                select_matches( matches, selectedBins, selectedBinsRev, filter, threshold);
            }
            catch (seqan::Exception const& e)
            {
                std::stringstream sstr;
                sstr << "Error counting kmers of sequence from " << this->id << " in IBF bins: " << e.what();
                logger->error(sstr.str());
                throw CountKmerException(sstr.str());
            }
        }
    }

    /**
        Filter matches based on strata_filter settings

        @matches:       unordered map to store matches (matches[binNo] = maxKmerCountBin)
        @len:           length of the current read
        @kmer_size:     kmer size used for building the IBF
        @strata_filter: -1 => no strata filter active, report everything up to min_kmers/max-error
                        n >= 0 => report all matches with maximum of n errors
        @return:        number of matches between reads and targets in ibfs
    */
    uint32_t Read::filter_matches(  TMatches& matches,
                                    uint16_t  len,
                                    uint16_t  kmer_size,
                                    int16_t   strata_filter )
    {

        uint16_t threshold_strata = 1; // minimum threshold (when strata_filter == -1)
        if ( strata_filter > -1 )
        {
            // get maximum possible number of error for this read
            uint16_t max_error = get_error( len, kmer_size, this->max_kmer_count );
            // get min kmer count neccessary to achieve the calculated number of errors
            threshold_strata = get_threshold_errors( len, kmer_size, max_error + strata_filter );
        }

        for ( auto const& v : matches )
        { // matches[binNo] = kmerCount
            if ( v.second >= threshold_strata )
            { // apply strata filter
                this->matches.push_back( ReadMatch{ v.first, v.second } );
            }
        }

        return this->matches.size();
    }

    /**
        find matches of kmers from the read within the bins of the IBFs
        if certain number of kmers within one bin were found, read is classified as host read
        @filters: vector of IBFs used to filter out reads
        @config: configuration settings needed
        @return: true, if at least one specific match was found, false otherwise 
        @throws: NullFilterException, ShortReadException, CountKmerException
    */
    bool Read::classify(  std::vector< TIbf >& filters, ClassifyConfig& config)
    {
        std::shared_ptr<spdlog::logger> logger = config.classification_logger;
        if (filters.empty())
        {
            logger->error("No IBF provided to classify the read!");
            throw NullFilterException("No IBF provided to classify the read!");
        }
        // k-mer sizes should be the same among filters
        uint16_t kmer_size = filters[0].kmerSize;
        TMatches matches;
        if ( seqan::length(this->seq) >= kmer_size )
        {
            // try to find kmer matches between read and ibfs
            try
            {
                find_matches( matches, filters, config );
            }
            catch (const CountKmerException& e)
            {
                throw;
            }

            if ( this->max_kmer_count > 0 )
            {
                //std::cout << "kmer count : " << this->max_kmer_count << std::endl;
                // filter matches based on strata filter settings
                return (filter_matches( matches, seqan::length(this->seq), kmer_size, config.strata_filter ) > 0) ;
            }
        }
        else
        {
            std::stringstream sstr;
            sstr << "Read " << this->id << " shorter than kmer size (" << kmer_size << ")";
            logger->error(sstr.str());
            throw ShortReadException("Read " + std::string(seqan::toCString(this->id)) + " shorter than kmer size");
        }
        return false;
    }
} // end namespace interleave