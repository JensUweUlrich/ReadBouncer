#include "IBF.hpp"

using namespace seqan;

namespace interleave
{

    /**
        return true if there is at least one bin with number of 
        matching kmers above the threshold
        @selectedBins:      number of kmer matches for every bin in the filter
        @selectedBinsRev:   number of reverse complement kmer matches for every bin in the filter
        @filter:            actual IBF
        @threshold:         minimum number of kmer matches needed to add the bin number with its kmer match count to the matches map 
    */
    bool Read::select_matches(std::vector< uint16_t >& selectedBins,
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
                return true;
            }
        }
        return false;
    }

    /**
        find matches of kmers in ibfs
        return true if read matches at least one bin of a given IBF
        @matches:   unordered map to store matches 
        @filters:   vector of IBFs
        @config:    configuration settings needed
        @throws:    CountKmerException
    */
    bool Read::find_matches(std::vector< TIbf >& filters, 
                                ClassifyConfig& config )
    {
        std::shared_ptr<spdlog::logger> logger = config.classification_logger;
        // iterate over all ibfs
        bool found = false;
        for ( TIbf& filter : filters )
        {
            try
            {
                // IBF count
                // find all bins that have at least 1 kmer in common with read_seq
                // selectedBins[binNo] = # of kmer matches for this bin
                std::vector< uint16_t > selectedBins    = seqan::count( filter, this->sequence );
                std::vector< uint16_t > selectedBinsRev = seqan::count( filter, TSeqRevComp( this->sequence ) );

                // get calculated threshold for minimum number of kmers needed to report a match
                // this is based on the confidence interval for mutated kmers in a read with expected error rate, kmer size
                TInterval ci = calculateCI(config.error_rate, filter.kmerSize, seqan::length(this->sequence), config.significance);
                //std::cout << "CI = [" << ci.first << " , " << ci.second << "]" << std::endl;
                uint16_t readlen = seqan::length(this->sequence);
                // minimum number of kmers = max number of kmer in read - upper bound of the CI
                //std::cout << seqan::length(this->seq) << " " << filter.kmerSize << std::endl;
                int16_t threshold = readlen - filter.kmerSize + 1 - ci.second;
                //uint16_t threshold = seqan::length(this->seq) - filter.kmerSize + 1 - (floor((ci.second - ci.first) / 2) + ci.first);
                // select matches above chosen threshold
                found = select_matches( selectedBins, selectedBinsRev, filter, threshold);
                if (found)
                    break;
            }
            catch (seqan::Exception const& e)
            {
                std::stringstream sstr;
                sstr << "Error counting kmers of sequence from " <<this->id << " in IBF bins: " << e.what();
                logger->error(sstr.str());
                throw CountKmerException(sstr.str());
            }
        }
        return found;
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
        if ( getReadLength() >= kmer_size )
        {
            // try to find kmer matches between read and ibfs
            try
            {
                return find_matches(filters, config );
            }
            catch (const CountKmerException& e)
            {
                throw;
            }
            // not needed anymore
            /* 
            if ( this->max_kmer_count > 0 )
            {
                //std::cout << "kmer count : " << this->max_kmer_count << std::endl;
                // filter matches based on strata filter settings
                return (filter_matches( matches, seqan::length(this->seq), kmer_size, config.strata_filter ) > 0) ;
            }
            */
        }
        else
        {
            std::stringstream sstr;
            sstr << "Read " << this->id << " shorter than kmer size (" << kmer_size << ")";
            logger->error(sstr.str());
            throw ShortReadException("Read " + this->id + " shorter than kmer size");
        }
        return false;
    }
} // end namespace interleave