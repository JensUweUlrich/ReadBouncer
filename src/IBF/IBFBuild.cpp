
#include "IBF.hpp"

namespace interleave
{
    std::string IbfClassificationLog, InterleavedBloomFilterLog = "";

    /**
        read reference sequences from files and store them in reference sequence queue
        @queue_refs:    thread safe queue that stores reference sequences to put in the ibf
        @config:        settings used to build the IBF (reference file path, number of references, kmer size, fragment length)
        @stats:         statistics (like runtime for specific tasks) for building the ibf
        @throws:        FileParserException  
    */
    std::future< void > IBF::parse_ref_seqs(SafeQueue< Seqs > &queue_refs, interleave::IBFConfig &config, FilterStats &stats)
    {
        
        this->ibf_logger->info("Starting to parse reference sequence files");
        this->ibf_logger->flush();
        if (config.reference_files.empty())
        {
            this->ibf_logger->error("No refernce file specified!");
            this->ibf_logger->flush();
            throw MissingReferenceFilesException("There were no reference files specified!");
        }

        std::future< void > read_task( std::async( std::launch::async, [=, &queue_refs, &stats] {
            // iterate over all reference sequence input files stated in config file
            std::shared_ptr<spdlog::logger> logger = spdlog::get("IbfLog");
            for ( std::string const& reference_file : config.reference_files )
            {
                logger->info("Start parsing sequences from " + reference_file);
                logger->flush();
                seqan::SeqFileIn seqFileIn;
                // open input file in first iteration
                if ( !seqan::open( seqFileIn, seqan::toCString( reference_file ) ) )
                {
                    logger->error("Could not open " + reference_file + ": Check existence or permissions!");
                    logger->flush();
                    throw FileParserException("Unable to open the file: " + reference_file );
                }
                // read current input file
                while ( !seqan::atEnd( seqFileIn ) )
                {

                    seqan::StringSet< seqan::CharString > ids;
                    seqan::StringSet< seqan::CharString > seqs;
                    // read all sequences from the current file
                    try
                    {
                        seqan::readRecords( ids, seqs, seqFileIn, config.n_refs );
                    }
                    catch ( seqan::Exception const& e )
                    {
                        logger->error("Problems parsing the file : " + reference_file + "[" + e.what() + "]");
                        logger->flush();
                        throw FileParserException( "ERROR: Problems parsing the file: " + reference_file + "[" + e.what() + "]");
                    }

                    //assert(seqan::length( ids ) == seqan::length( seq ) )

                    // iterate over all sequences loaded into the sequence bin
                    logger->info("Fragment reference sequences based on NNNs");
                    logger->flush();
                    for ( uint64_t i = 0; i < seqan::length( ids ); ++i )
                    {
                        stats.totalSeqsFile += 1;
                        // sequence with length < kmer size is invalid and will not be added to the queue
                        if ( seqan::length( seqs[i] ) < config.kmer_size )
                        {
                            stats.invalidSeqs += 1;
                            continue;
                        }
                        std::string cid   = seqan::toCString( ids[i] );
                        // delete everything after first space in seq identifier
                        std::string seqid = cid.substr( 0, cid.find( ' ' ) );
                        // add reference sequences to the queue
                        int counter = 1;
                        std::string seq = std::string(seqan::toCString(seqs[i]));
			// remove all Ns from the sequence
			std::stringstream buf; 
                        for (std::string seq : cutOutNNNs(seq, seqan::length(seqs[i])))
                        {
                            buf << seq;			   
                        }
			std::string newseq = buf.str();
			queue_refs.push(Seqs{ seqid, ((seqan::Dna5String) newseq) });
			// calculate bins needed for that sequence
	                stats.totalBinsBinId += ( newseq.length() / config.fragment_length) + 1;
                        stats.sumSeqLen += newseq.length();
                    }
                }
                seqan::close( seqFileIn );
                logger->info("Finished parsing sequences from " + reference_file);
                logger->flush();
            }
            logger->info("Finished to parse all reference sequence files");
            logger->flush();
            // notify all threads that the last ref sequences were pushed into the queue
            queue_refs.notify_push_over();
        } ) );
        return read_task;
    }

    /**
    * split underlying sequence based on stretches of N's
    * @seq      : reference sequence as a string
    * @seqlen   : length of the sequence
    * @return   : vector of subsequences, that results from splitting the original sequence
    */
    std::vector<std::string> IBF::cutOutNNNs(std::string& seq, uint64_t seqlen)
    {
        std::vector<std::string> splittedStrings;
        size_t start = 0;
        size_t end = 0;

        while ((start = seq.find_first_not_of("N", end)) != std::string::npos)
        {
            end = seq.find("N", start);
            if (end > seqlen)
            {
                std::string s = seq.substr(start, seqlen - start - 1);
                splittedStrings.push_back(s);
                break;
            }
            std::string s = seq.substr(start, end - start);
            
            splittedStrings.push_back(s);
        }
        return splittedStrings;
    }

    /**
        split reference sequence into fragments of predefined length and add the fragments to the IBF
        @tasks:         vector of asynchronous tasks executed in parallel, where every task splits one reference into fragments and adds them to the IBF
        @config:        settings used to build the IBF (number of threads, kmer size, fragment length)
        @binid:         reference to a binid counter that is incremented for every new fragment
        @queue_refs:    thread safe queue that stores reference sequences to put in the ibf
        @throws:        IBFBuildException
        TODO: is binid thread safe? mutex needed?  
    */
    void IBF::add_sequences_to_filter( std::vector< std::future< void > >& tasks, 
                                        IBFConfig &config, 
                                        uint64_t &binid,
                                        SafeQueue< Seqs > &queue_refs)
    {
        this->ibf_logger->info("Start adding fragmented sequences to Interleaved Bloom Filter");
        std::stringstream sstr;
        sstr << "Maximum fragment size stored in each bin : " << config.fragment_length;
        this->ibf_logger->info(sstr.str());
        this->ibf_logger->flush();
        for ( uint16_t taskNo = 0; taskNo < config.threads_build; ++taskNo )
        {
            tasks.emplace_back( std::async( std::launch::async, [=, &queue_refs, &binid] {
                std::shared_ptr<spdlog::logger> logger = spdlog::get("IbfLog");
                while ( true )
                {
                    // take ref seq from the queue
                    Seqs val = queue_refs.pop();
                    if ( val.seqid != "" )
                    {
                        // add all fragments of the current reference to the IBF
                        int64_t fragIdx = 0;
                        int64_t fragstart = fragIdx * config.fragment_length - config.overlap_length + 1;
                        if (fragstart < 0) fragstart = 0;


                        int64_t seqlen = (int64_t) (length(val.seq));
                        while ( fragstart < (seqlen - 1))
                        {
                            // consecutive fragments need to overlap by kmer_size-1 nts
                            // otherwise we would miss to include kmers spanning the borders of the fragments
                                
                            // first fragment starts at sequence index 0
                            uint64_t fragend = (fragIdx+1) * config.fragment_length;
                            // make sure that last fragment ends at last position of the reference sequence
                            if (fragend > length(val.seq)) fragend = length(val.seq);
                            //auto [fragstart, fragend, binid] = seq_bin.at( val.seqid )[i];
                            // fragment of the reference is defined as infix
                            // For infixes, we have to provide both the including start and the excluding end position.
                            // fragstart -1 to fix offset
                            // fragend -1+1 to fix offset and not exclude last position

                            
                            try
                            {
                                seqan::Infix< seqan::Dna5String >::Type fragment = seqan::infix( val.seq, fragstart, fragend );
                                seqan::insertKmer( this->filter, fragment, binid++ );
                            }
                            catch (seqan::Exception const& e )
                            {
                                logger->error("Error inserting the sequence of " + val.seqid + " to the IBF!");
                                std::stringstream sstr;
                                sstr << "Sequence Length=" << seqan::length(val.seq) << " , Fragment start=" << fragstart << " , Fragment end=" << fragend;
                                logger->error(sstr.str());
                                logger->flush();
                                throw InsertSequenceException("Error inserting the sequence of " + val.seqid + " to the IBF!");
                            }
                            fragIdx++;
                            fragstart = fragIdx * config.fragment_length - config.kmer_size + 1;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            } ) );
        }
    }

    /**
        update existing IBF by adding additional sequences to it
        @config: settings used to load and update the IBF (update_filter_file, reference_sequence_files etc.)
        @return: all statistics describing the process of updating the IBF
        @throws: FileParserException, StoreFilterException, IBFBuildException
    */
    FilterStats IBF::update_filter(IBFConfig& config)
    {
        this->ibf_logger->info("Start updating Interleaved Bloom Filter from file");
        this->ibf_logger->flush();
        if (!config.validate())
        {
            ibf_logger->error("Given IBFConfig is not valid");
            logIBFConfig(config);
            throw InvalidConfigException("Config not valid!");
        }

        FilterStats stats;
        stats.timeIBF.start();

        // load existing filter from file
        try
        {
            stats = load_filter(config);
        }
        catch (const FileParserException &e)
        {
            throw;
        }

        // config.n_batches*config.n_refs = max. amount of references in memory
        SafeQueue< Seqs > queue_refs(config.n_batches ); 

        // Start extra thread for reading the input
        try{
            stats.timeLoadSeq.start();
            std::future< void > read_task = parse_ref_seqs(queue_refs, config, stats);
            read_task.get();
            stats.timeLoadSeq.stop();
        }
        catch(const FileParserException &e)
        {
            throw;
        }    

        this->ibf_logger->info("Re-sizing Interleaved Bloom Filter ");
        this->ibf_logger->flush();

        stats.timeBuild.start();
        // add new bins only if reference file is provided
        // number of bins calculated by length of input ref sequences and fragment length
        // number of new bins = number of old bins + number of bins from new ref seqs
        uint32_t number_new_bins = stats.totalBinsBinId + stats.totalBinsFile;
        if ( number_new_bins > stats.totalBinsFile )
        {
            // just resize if number of bins is bigger than amount on IBF
            // when updating an IBF with empty bins or removing the last bins, this will not be true
            this->filter.resizeBins( number_new_bins );
            // number of new bins from ref seq file
            stats.newBins = stats.totalBinsBinId;
            // total number of all bins in the filter
            stats.totalBinsBinId = number_new_bins;
        } // if new bins are smaller (less bins, sequences removed) IBF still keep all bins but empty

        // Start execution threads to add kmers
        uint64_t binid = stats.totalBinsFile;
        std::vector< std::future< void > > tasks;
        try
        {
            add_sequences_to_filter(tasks, config, binid, queue_refs);
            for ( auto& task : tasks )
            {
                task.get();
            }
        }
        catch(const IBFBuildException &e)
        {
            throw;
        }
        stats.timeBuild.stop();

        this->ibf_logger->info("Finished adding fragmented sequences to Interleaved Bloom Filter");
        ibf_logger->flush();

        this->ibf_logger->info("Writing updated IBF to file : " + config.update_filter_file);
        this->ibf_logger->flush();
        // Store filter
        stats.timeSaveFilter.start();
        try
        {
            seqan::store( this->filter, seqan::toCString( config.update_filter_file ) );
        }
        catch (seqan::Exception const& e )
        {
            std::stringstream sstr;
            sstr << "Could not store updated IBF in " << config.update_filter_file << ":" << e.what();
            this->ibf_logger->error(sstr.str());
            this->ibf_logger->flush();
            throw StoreFilterException("Could not store IBF to " + config.update_filter_file + ":" + e.what());
        }
        stats.timeSaveFilter.stop();
        stats.timeIBF.stop();

        return stats;
    }

    /**
        load an IBF from given update or input filter file#
        @config: settings used to load the IBF (update_filter_file, input_filter_file)
        @return: all statistics describing the process of loading the IBF
        @throws: ParseIBFFileException, MissingIBFFileException
    */
    FilterStats IBF::load_filter( IBFConfig& config)
    {
        this->ibf_logger->info("Start loading Interleaved Bloom Filter from file");
        this->ibf_logger->flush();
        FilterStats stats;
        // take time needed to load the filter
        stats.timeLoadFilter.start();
        
        // load from disk in case of update
        if (!config.update_filter_file.empty())
        {
            // load filter for updating
            try
            {
                seqan::retrieve( this->filter, seqan::toCString( config.update_filter_file ) );
            }
            catch(seqan::Exception const& e)
            {
                std::stringstream sstr;
                sstr << "Error parsing IBF update file " << config.update_filter_file << ": " << e.what();
                this->ibf_logger->error(sstr.str());
                this->ibf_logger->flush();
                throw ParseIBFFileException("Error parsing IBF update file " + config.update_filter_file + ": " + e.what());
            }
            
        }
        else if(!config.input_filter_file.empty())
        {
            // load filter without updating
            try
            {
                seqan::retrieve( this->filter, seqan::toCString( config.input_filter_file ) );
            }
            catch(seqan::Exception const& e)
            {
                std::stringstream sstr;
                sstr << "Error parsing IBF input file " << config.input_filter_file << ": " << e.what();
                this->ibf_logger->error(sstr.str());
                this->ibf_logger->flush();
                throw ParseIBFFileException("Error parsing IBF input file " + config.input_filter_file + ": " + e.what());
            }
        }
        else
        {
            this->ibf_logger->error("Error: Either update_filter_file or input_filter_file have to be specified.");
            this->ibf_logger->flush();
            throw MissingIBFFileException("Error: Either update_filter_file or input_filter_file have to be specified.");
        }


        // totalBinsFile account for all bins, even empty
        stats.totalBinsFile = seqan::getNumberOfBins( this->filter );
        config.kmer_size    = seqan::getKmerSize( this->filter );
        // config.hash_functions = seqan::get...( filter ); // not avail.
        // config.filter_size_bits = seqan::get...( filter ); // not avail.

        this->ibf_logger->info("Successfully loaded IBF from file");
        std::stringstream sstr;
        sstr << "Kmer size      : " << config.kmer_size;
        this->ibf_logger->info(sstr.str());
        sstr.str("");
        sstr << "Number of Bins : " << stats.totalBinsFile;
        this->ibf_logger->info(sstr.str());
        this->ibf_logger->flush();
        stats.timeLoadFilter.stop();
        
        return stats;
    }

    /**
    *   caalculate optimal IBF size based on number of bins and maximum false positive rate
    *   @config         : settings used to build the IBF (kmer size, fragment length etc.)
    *   @numberOfBins   : Actual number of bins stored in the IBF
    *   @return         : optimal size of the IBF in Bits
    */
    uint64_t IBF::calculate_filter_size_bits(IBFConfig& config, uint64_t numberOfBins)
    {
        uint64_t max_kmer_count = config.fragment_length - config.kmer_size + 1;
        uint64_t optimalNumberOfBins = floor(((double) numberOfBins / 64.0) + 1) * 64;
        uint64_t BinSizeBits = ceil(-1 / (pow(1 - pow((double) config.max_fp, 1.0 / (double) config.hash_functions),
                                                        1.0 / ((double) (config.hash_functions * max_kmer_count))) - 1));
        //std::cout << "OptimalNumberOfBins : " << optimalNumberOfBins << " , " << "BinSizeBits : " << BinSizeBits << std::endl;
        
        return BinSizeBits * optimalNumberOfBins;
    }

    /**
        create IBF based on given reference sequence files and other configuration settings
        @config:    settings used to build the IBF (reference files, number of threads, kmer size, fragment length etc.)
        @return:    all statistics describing the process of building this IBF
        @throws:    FileParserException, StoreFilterException, InvalidConfigException, IBFBuildException
    */
    FilterStats IBF::create_filter(IBFConfig& config)
    {
        this->ibf_logger->info("Start creating Interleaved Bloom Filter");
        this->ibf_logger->flush();
        if (!config.validate())
        {
            this->ibf_logger->error("Given IBFConfig is not valid");
            logIBFConfig(config);
            throw InvalidConfigException("Config not valid!");
        }
        FilterStats stats;
        
        stats.timeIBF.start();

        // config.n_batches*config.n_refs = max. amount of references in memory
        SafeQueue< Seqs > queue_refs(config.n_batches ); 

        // Start extra thread for reading the input
        
        try{
            stats.timeLoadSeq.start();
            std::future< void > read_task = parse_ref_seqs(queue_refs, config, stats);
            read_task.get();
            stats.timeLoadSeq.stop();
        }
        catch(const FileParserException &e)
        {
            throw;
        }

        this->ibf_logger->info("Calculating Interleaved Bloom Filter size");
        config.filter_size_bits = calculate_filter_size_bits(config, stats.totalBinsBinId);
        std::stringstream sstr;
        sstr << "Calculated size of the IBF : " << (double)config.filter_size_bits / (double)config.MBinBits << " MBytes";
        if (config.verbose)
            std::cout << sstr.str()  << std::endl;

        this->ibf_logger->info(sstr.str());
        this->ibf_logger->flush();
        stats.timeBuild.start();
        try
        {
            this->filter = TIbf(stats.totalBinsBinId, config.hash_functions, config.kmer_size, config.filter_size_bits);
            stats.totalBinsFile = seqan::getNumberOfBins(this->filter);
        }
        catch (seqan::Exception const& e)
        {
            this->ibf_logger->error("Could not instantiate IBF Filter");
            sstr.str("");
            sstr << "Number of Bins : " << stats.totalBinsBinId << ", Number of hash functions : " << config.hash_functions;
            this->ibf_logger->error(sstr.str());
            sstr.str("");
            sstr << "Kmer size : " << config.kmer_size << ", IBF size in bits : " << config.filter_size_bits;
            this->ibf_logger->error(sstr.str());
            this->ibf_logger->flush();
            throw NullFilterException("Could not instantiate IBF Filter");
        }
        // Start execution threads to add kmers
        uint64_t binid = 0;
        std::vector< std::future< void > > tasks;
        try
        {
            add_sequences_to_filter(tasks, config, binid, queue_refs);
            for ( auto& task : tasks )
            {
                task.get();
            }
        }
        catch(const IBFBuildException &e)
        {
            throw;
        }
        stats.timeBuild.stop();
        this->ibf_logger->info("Finished adding fragmented sequences to Interleaved Bloom Filter");

        this->ibf_logger->info("Writing IBF to file : " + config.output_filter_file);
        this->ibf_logger->flush();
        // Store filter
        stats.timeSaveFilter.start();
        try{
            seqan::store( this->filter, seqan::toCString( config.output_filter_file ) );
        }
        catch (seqan::Exception const& e )
        {
            sstr.str("");
            sstr << "Could not store IBF to " << config.output_filter_file << ":" << e.what();
            this->ibf_logger->error(sstr.str());
            this->ibf_logger->flush();
            throw StoreFilterException("Could not store IBF to " + config.output_filter_file + ":" + e.what());
        }
        stats.timeSaveFilter.stop();

        stats.timeIBF.stop();
        this->ibf_logger->info("Successfully stored IBF to file.");
        this->ibf_logger->flush();
        return stats;
    }

    /**
        print time statistics for loading, building and updating an IBF
        @config: settings used to build, update or load the IBF (reference files, number of threads, kmer size, fragment length etc.)
        @stats:  all statistics describing the process
    */
    void print_time( const interleave::IBFConfig& config, interleave::FilterStats& stats )
    {
        using ::operator<<;

        std::cerr << "IBF-build         start time: " << stats.timeIBF.begin() << std::endl;
        std::cerr << "Loading files     start time: " << stats.timeLoadFiles.begin() << std::endl;
        std::cerr << "Loading files       end time: " << stats.timeLoadFiles.end() << std::endl;
        std::cerr << "Loading sequences start time: " << stats.timeLoadSeq.begin() << std::endl;
        std::cerr << "Loading filter    start time: " << stats.timeLoadFilter.begin() << std::endl;
        std::cerr << "Loading filter      end time: " << stats.timeLoadFilter.end() << std::endl;
        std::cerr << "Building filter   start time: " << stats.timeBuild.begin() << std::endl;
        std::cerr << "Loading sequences   end time: " << stats.timeLoadSeq.end() << std::endl;
        std::cerr << "Building filter     end time: " << stats.timeBuild.end() << std::endl;
        std::cerr << "Saving filter     start time: " << stats.timeSaveFilter.begin() << std::endl;
        std::cerr << "Saving filter       end time: " << stats.timeSaveFilter.end() << std::endl;
        std::cerr << "IBF-build         end time: " << stats.timeIBF.end() << std::endl;
        std::cerr << std::endl;
        std::cerr << " - loading files: " << stats.timeLoadFiles.elapsed() << std::endl;
        std::cerr << " - loading filter: " << stats.timeLoadFilter.elapsed() << std::endl;
        std::cerr << " - loading sequences (1t): " << stats.timeLoadSeq.elapsed() << std::endl;
        std::cerr << " - building filter (" << config.threads_build << "t): " << stats.timeBuild.elapsed() << std::endl;
        std::cerr << " - saving filter: " << stats.timeSaveFilter.elapsed() << std::endl;
        std::cerr << " - total: " << stats.timeIBF.elapsed() << std::endl;
        std::cerr << std::endl;
    }

    /**
        print all statistics describing the creation of the IBF
        @stats: all statistics describing the build process of the IBF
    */
    void print_build_stats( interleave::FilterStats& stats)
    {
        double   elapsed_build = stats.timeBuild.elapsed();
        uint64_t validSeqs     = stats.totalSeqsFile - stats.invalidSeqs;
        std::cerr << "IBF-build processed " << validSeqs << " sequences (" << stats.sumSeqLen / 1000000.0 << " Mbp) in "
                << elapsed_build << " seconds (" << ( validSeqs / 1000.0 ) / ( elapsed_build / 60.0 ) << " Kseq/m, "
                << ( stats.sumSeqLen / 1000000.0 ) / ( elapsed_build / 60.0 ) << " Mbp/m)" << std::endl;
        if ( stats.invalidSeqs > 0 )
            std::cerr << " - " << stats.invalidSeqs << " invalid sequences were skipped" << std::endl;
        if ( stats.newBins > 0 )
            std::cerr << " - " << stats.newBins << " bins were added to the IBF" << std::endl;
        std::cerr << " - " << validSeqs << " sequences in " << stats.totalBinsFile + stats.newBins
                << " bins were written to the IBF" << std::endl;
    }

    void print_load_stats(interleave::FilterStats& stats)
    {
        double   elapsed_load = stats.timeLoadFilter.elapsed();
        std::cerr << stats.totalBinsFile << " bins were loaded in " << elapsed_load << " seconds from the IBF" << std::endl;
    }

} // namespace interleave
