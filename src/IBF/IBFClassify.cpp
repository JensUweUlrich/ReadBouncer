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
        return true if there is at least one bin with number of
        matching kmers above the threshold
        @selectedBins:      number of kmer matches for every bin in the filter
        @selectedBinsRev:   number of reverse complement kmer matches for every bin in the filter
        @filter:            actual IBF
        @threshold:         minimum number of kmer matches needed to add the bin number with its kmer match count to the matches map
    */
    uint64_t Read::max_matches(std::vector< uint16_t >& selectedBins,
        std::vector< uint16_t >& selectedBinsRev,
        TIbf& filter,
        uint16_t                 threshold)
    {
        // for each bin
        // for ( uint32_t binNo = 0; binNo < filter.ibf.noOfBins; ++binNo )
        // loop in map structure to avoid extra validations when map.size() < filter.ibf.noOfBins when ibf is updated and
        // sequences removed also avoid the error of spurius results from empty bins (bug reported)
        uint64_t max_kmer_count = 0;
        for (uint32_t binNo = 0; binNo < filter.noOfBins; ++binNo)
        {
            //std::cout<<selectedBins[binNo] << " , " << selectedBinsRev[binNo] <<std::endl;
            // if kmer count is higher than threshold
            if (selectedBins[binNo] >= threshold || selectedBinsRev[binNo] >= threshold)
            {
                if (selectedBins[binNo] > max_kmer_count)
                    max_kmer_count = selectedBins[binNo];
                if (selectedBinsRev[binNo] > max_kmer_count)
                    max_kmer_count = selectedBinsRev[binNo];
            }
        }
        return max_kmer_count;
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
                //QMessageBox::critical(NULL, "Error", QString::fromStdString(sstr.str()));
                throw CountKmerException(sstr.str());
            }
        }
        return found;
    }

    /**
        find matches of kmers in ibfs
        return true if read matches at least one bin of a given IBF
        @matches:   unordered map to store matches
        @filters:   vector of IBFs
        @config:    configuration settings needed
        @throws:    CountKmerException
    */
    uint64_t Read::count_matches(IBFMeta& filter, ClassifyConfig& config)
    {
         std::shared_ptr<spdlog::logger> logger = config.classification_logger;
        // iterate over all ibfs

        try
        {
             // IBF count
            // find all bins that have at least 1 kmer in common with read_seq
            // selectedBins[binNo] = # of kmer matches for this bin
            std::vector< uint16_t > selectedBins = seqan::count(filter.filter, this->sequence);
            std::vector< uint16_t > selectedBinsRev = seqan::count(filter.filter, TSeqRevComp(this->sequence));

            // get calculated threshold for minimum number of kmers needed to report a match
            // this is based on the confidence interval for mutated kmers in a read with expected error rate, kmer size
            TInterval ci = calculateCI(config.error_rate, filter.filter.kmerSize, seqan::length(this->sequence), config.significance);
            //std::cout << "CI = [" << ci.first << " , " << ci.second << "]" << std::endl;
            uint16_t readlen = seqan::length(this->sequence);
            // minimum number of kmers = max number of kmer in read - upper bound of the CI
            //std::cout << seqan::length(this->seq) << " " << filter.kmerSize << std::endl;
            int16_t threshold = readlen - filter.filter.kmerSize + 1 - ci.second;
            //uint16_t threshold = seqan::length(this->seq) - filter.kmerSize + 1 - (floor((ci.second - ci.first) / 2) + ci.first);
            // select matches above chosen threshold
            return max_matches(selectedBins, selectedBinsRev, filter.filter, threshold);
        }
        catch (seqan::Exception const& e)
        {
            std::stringstream sstr;
            sstr << "Error counting kmers of sequence from " << this->id << " in IBF bins: " << e.what();
            logger->error(sstr.str());
            //QMessageBox::critical(NULL, "Error", QString::fromStdString(sstr.str()));
            throw CountKmerException(sstr.str());
        }
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
            //QMessageBox::critical(NULL, "Error", "No IBF provided to classify the read!");
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
                //QMessageBox::critical(NULL, "Error", e.what());
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

            //QMessageBox::critical(NULL, "Error", QString::fromStdString(sstr.str()));

            throw ShortReadException("Read " + this->id + " shorter than kmer size");
        }
        return false;
    }

    

    /**
        find matches of kmers from the read within the bins of the IBFs
        if certain number of kmers within one bin were found, read is considered classified
        account read to the best matching filter
        @filters: vector of IBFs used to filter out reads
        @config: configuration settings needed
        @return: index of best matching IBF in filters vector
        @throws: NullFilterException, ShortReadException, CountKmerException
    */
    int Read::classify(std::vector< IBFMeta >& filters, ClassifyConfig& config)
    {
        std::shared_ptr<spdlog::logger> logger = config.classification_logger;
        if (filters.empty())
        {
            logger->error("No IBF provided to classify the read!");
            //QMessageBox::critical(NULL, "Error", "No IBF provided to classify the read!");
            throw NullFilterException("No IBF provided to classify the read!");
        }
        // k-mer sizes should be the same among filters
        uint16_t kmer_size = filters[0].filter.kmerSize;
        TMatches matches;

        if (seqan::length(this->sequence) >= kmer_size)
        {
            // try to find kmer matches between read and ibfs
            try
            {
                std::vector< std::future< uint64_t > > tasks;
                for (uint8_t ind = 0; ind < filters.size(); ++ind)
                {
                    tasks.emplace_back(std::async(std::launch::async, &interleave::Read::count_matches, this, std::ref(filters[ind]), std::ref(config)));
                }

                uint64_t max_kmer_count = 0;
                int best_index = -1;
                for (int ind = 0; ind < tasks.size();++ind)
                {
                    uint64_t count = tasks[ind].get();
                    if (count > max_kmer_count)
                    {
                        best_index = ind;
                        max_kmer_count = count;
                    }
                }
                return best_index;
            }
            catch (const CountKmerException& e)
            {
                //QMessageBox::critical(NULL, "Error", e.what());
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
            logger->error(sstr.str());//

            //QMessageBox::critical(NULL, "Error", QString::fromStdString(sstr.str()));
            throw ShortReadException("Read " + std::string(seqan::toCString(this->id)) + " shorter than kmer size");
        }
        return -1;
    }

    std::pair<int, int> Read::classify(std::vector< IBFMeta >& filt1, std::vector< IBFMeta >& filt2, ClassifyConfig& config)
    {
        std::shared_ptr<spdlog::logger> logger = config.classification_logger;
        if (filt1.empty() || filt2.empty())
        {
            logger->error("No IBF provided to classify the read!");
            throw NullFilterException("No IBF provided to classify the read!");
            //QMessageBox::critical(NULL, "Error", "No IBF provided to classify the read!");
        }
        // k-mer sizes should be the same among filters

        TMatches matches;
        std::pair<uint64_t, uint64_t> result = std::make_pair(0,0);

        // try to find kmer matches between read and ibfs
        try
        {
            std::vector< std::future< uint64_t > > tasks;
            for (uint8_t ind = 0; ind < filt1.size(); ++ind)
            {
                if (seqan::length(this->sequence) >= filt1[ind].filter.kmerSize)
                    tasks.emplace_back(std::async(std::launch::async, &interleave::Read::count_matches, this, std::ref(filt1[ind]), std::ref(config)));
            }

            uint64_t max_kmer_count = 0;
            int best_index = -1;
            for (int ind = 0; ind < tasks.size(); ++ind)
            {
                uint64_t count = tasks[ind].get();
                if (count > max_kmer_count)
                {
                    best_index = ind;
                    max_kmer_count = count;
                }
            }

            tasks.clear();

            result.first = max_kmer_count;
            //std::cerr << result.first << std::endl;
            for (uint8_t ind = 0; ind < filt2.size(); ++ind)
            {
                if (seqan::length(this->sequence) >= filt2[ind].filter.kmerSize)
                    tasks.emplace_back(std::async(std::launch::async, &interleave::Read::count_matches, this, std::ref(filt2[ind]), std::ref(config)));
            }

            max_kmer_count = 0;
            best_index = -1;
            for (int ind = 0; ind < tasks.size(); ++ind)
            {
                uint64_t count = tasks[ind].get();
                if (count > max_kmer_count)
                {
                    best_index = ind;
                    max_kmer_count = count;
                }
            }
            result.second = max_kmer_count;
            //std::cerr << result.second << std::endl;

        }
        catch (const CountKmerException& e)
        {
            //QMessageBox::critical(NULL, "Error", e.what());
            throw;
        }

        return result;
    }

} // end namespace interleave
