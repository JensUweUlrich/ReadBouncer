
#include "IBFBuild.hpp"

namespace interleave
{


/*  @DEPRECATED
    parse tsv file with information about reference sequences
    seqid <tab> seqstart <tab> seqend <tab> binid
    sequences with same binid will be put in same bloom filter ???
    sometimes only specific parts of the reference should be added to the IBF
    these fragments can be specified in seq_id_bin file
    eg
    NC12345 1       1938    1
    NC12345 2561    7145    1
*/
void IBF::parse_seqid_bin( const std::string& seqid_bin_file, TSeqBin& seq_bin, std::set< uint64_t >& bin_ids )
{
    std::string   line;
    std::ifstream infile( seqid_bin_file );
    while ( std::getline( infile, line, '\n' ) )
    {
        std::istringstream         stream_line( line );
        std::vector< std::string > fields;
        std::string                field;
        while ( std::getline( stream_line, field, '\t' ) )
            fields.push_back( field );
        // seqid <tab> seqstart <tab> seqend <tab> binid
        uint32_t binid = std::stoul( fields[3] );
        // add fragmentbin to seq_bins
        seq_bin[fields[0]].push_back( FragmentBin{ std::stoul( fields[1] ), std::stoul( fields[2] ), binid } );
        bin_ids.insert( binid );
    }
}

/**
    read reference sequences from files and store them in reference sequence queue
*/

std::future< void > IBF::parse_ref_seqs(SafeQueue< Seqs > &queue_refs, std::mutex &mtx, interleave::Config &config, Stats &stats)
{
    std::future< void > read_task( std::async( std::launch::async, [=, &queue_refs, &mtx, &stats] {
        // iterate over all reference sequence input files stated in config file
        for ( std::string const& reference_file : config.reference_files )
        {
            seqan::SeqFileIn seqFileIn;
            // open input file in first iteration
            if ( !seqan::open( seqFileIn, seqan::toCString( reference_file ) ) )
            {
                std::cerr << "ERROR: Unable to open the file: " << reference_file << std::endl;
                continue;
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
                    std::scoped_lock lock( mtx );
                    std::cerr << "ERROR: Problems parsing the file: " << reference_file << "[" << e.what() << "]"
                              << std::endl;
                }

                // iterate over all sequences loaded into the sequence bin
                for ( uint64_t i = 0; i < seqan::length( ids ); ++i )
                {
                    stats.totalSeqsFile += 1;
                    // sequence with length < kmer size is invalid and will not be added to the queue
                    if ( seqan::length( seqs[i] ) < config.kmer_size )
                    {
                        if ( config.verbose )
                        {
                            std::scoped_lock lock( mtx );
                            std::cerr << "WARNING: sequence smaller than k-mer size"
                                      << " [" << ids[i] << "]" << std::endl;
                        }
                        stats.invalidSeqs += 1;
                        continue;
                    }
                    std::string cid   = seqan::toCString( ids[i] );
                    // delete everything after first space in seq identifier
                    std::string seqid = cid.substr( 0, cid.find( ' ' ) );
                    // don't add sequences whose id is not in the seq_bin
                    // DEPRECATED: since we build seq_bin based on input reference sequences
                    /*if ( seq_bin.count( seqid ) == 0 )
                    {
                        if ( config.verbose )
                        {
                            std::scoped_lock lock( mtx );
                            std::cerr << "WARNING: sequence not defined on seqid-bin-file"
                                      << " [" << seqid << "]" << std::endl;
                        }
                        stats.invalidSeqs += 1;
                        continue;
                    }
                    */
                    stats.sumSeqLen += seqan::length( seqs[i] );
                    stats.totalBinsBinId += (seqan::length( seqs[i] ) / config.fragment_length) + 1;
                    // add reference sequences to the queue
                    queue_refs.push( Seqs{ seqid, seqs[i] } );
                }
            }
            seqan::close( seqFileIn );
        }
        // notify all threads that the last ref sequences were pushed into the queue
        queue_refs.notify_push_over();
    } ) );

    return read_task;
}


std::vector< std::future< void > > IBF::add_sequences_to_filter(Config &config, uint64_t &binid, SafeQueue< Seqs > &queue_refs)
{
    std::vector< std::future< void > > tasks;
    Tfilter filter = this->filter;
    for ( uint16_t taskNo = 0; taskNo < config.threads_build; ++taskNo )
    {
        tasks.emplace_back( std::async( std::launch::async, [=, &filter, &queue_refs, &binid] {
            while ( true )
            {
                // take ref seq from the queue
                Seqs val = queue_refs.pop();
                if ( val.seqid != "" )
                {
                    // add all fragments of the current reference to the IBF
                    for ( uint64_t i = 0; (i * config.fragment_length - config.kmer_size + 1) < length(val.seq) - 1; i++ )
                    {
                        // consecutive fragments need to overlap by kmer_size-1 nts
                        // otherwise we would miss to include kmers spanning the borders of the fragments
                        uint64_t fragstart = i * config.fragment_length - config.kmer_size + 1;
                        // first fragment starts at sequence index 0
                        if (fragstart < 0) fragstart = 0;
                        uint64_t fragend = (i+1) * config.fragment_length;
                        // make sure that last fragment ends at last position of the reference sequence
                        if (fragend > length(val.seq)) fragend = length(val.seq);
                        //auto [fragstart, fragend, binid] = seq_bin.at( val.seqid )[i];
                        // fragment of the reference is defined as infix
                        // For infixes, we have to provide both the including start and the excluding end position.
                        // fragstart -1 to fix offset
                        // fragend -1+1 to fix offset and not exclude last position
                        seqan::Infix< seqan::Dna5String >::Type fragment = seqan::infix( val.seq, fragstart, fragend );
                        std::cout<<fragment<<::std::endl;
                        seqan::insertKmer( filter, fragment, ++binid );
                    }
                }
                else
                {
                    break;
                }
            }
        } ) );
    }
    return tasks;
}

void IBF::load_filter( Config& config, const std::set< uint64_t >& bin_ids, Stats& stats )
{

    
    // load from disk in case of update
    if ( !config.update_filter_file.empty() )
    {
        // load filter
        seqan::retrieve( this->filter, seqan::toCString( config.update_filter_file ) );
        // totalBinsFile account for all bins, even empty
        stats.totalBinsFile = seqan::getNumberOfBins( this->filter );
        config.kmer_size    = seqan::getKmerSize( this->filter );
        // config.hash_functions = seqan::get...( filter ); // not avail.
        // config.filter_size_bits = seqan::get...( filter ); // not avail.

        // last element (set is ordered) plus one

        uint32_t number_new_bins = *bin_ids.rbegin() + 1;
        if ( number_new_bins > stats.totalBinsFile )
        {
            // just resize if number of bins is bigger than amount on IBF
            // when updating an IBF with empty bins or removing the last bins, this will not be true
            this->filter.resizeBins( number_new_bins );
            stats.newBins = number_new_bins - stats.totalBinsFile;
        } // if new bins are smaller (less bins, sequences removed) IBF still keep all bins but empty

        // Reset bins if complete set of sequences is provided (re-create updated bins)
        if ( config.update_complete )
        {
            std::vector< uint32_t > updated_bins;
            // For all binids in the file provided, only clean bins for the old bins
            // new bins are already cleared
            for ( auto const& binid : bin_ids )
            {
                if ( binid >= stats.totalBinsFile - 1 )
                {
                    break;
                }
                updated_bins.emplace_back( binid );
                // std::cerr << "Cleared: " << binid << std::endl;
            }
            seqan::clear( this->filter, updated_bins, config.threads ); // clear modified bins
        }
    }
    else
    {
        this->filter = Tfilter( stats.totalBinsBinId, config.hash_functions, config.kmer_size, config.filter_size_bits );
        stats.totalBinsFile = seqan::getNumberOfBins( this->filter );
    }

}

void print_time( const interleave::Config& config,
                 const StopClock&          timeIBF,
                 const StopClock&          timeLoadFiles,
                 const StopClock&          timeLoadSeq,
                 const StopClock&          timeBuild,
                 const StopClock&          timeLoadFilter,
                 const StopClock&          timeSaveFilter )
{
    using ::operator<<;

    std::cerr << "IBF-build         start time: " << timeIBF.begin() << std::endl;
    std::cerr << "Loading files     start time: " << timeLoadFiles.begin() << std::endl;
    std::cerr << "Loading files       end time: " << timeLoadFiles.end() << std::endl;
    std::cerr << "Loading sequences start time: " << timeLoadSeq.begin() << std::endl;
    std::cerr << "Loading filter    start time: " << timeLoadFilter.begin() << std::endl;
    std::cerr << "Loading filter      end time: " << timeLoadFilter.end() << std::endl;
    std::cerr << "Building filter   start time: " << timeBuild.begin() << std::endl;
    std::cerr << "Loading sequences   end time: " << timeLoadSeq.end() << std::endl;
    std::cerr << "Building filter     end time: " << timeBuild.end() << std::endl;
    std::cerr << "Saving filter     start time: " << timeSaveFilter.begin() << std::endl;
    std::cerr << "Saving filter       end time: " << timeSaveFilter.end() << std::endl;
    std::cerr << "ganon-build         end time: " << timeIBF.end() << std::endl;
    std::cerr << std::endl;
    std::cerr << " - loading files: " << timeLoadFiles.elapsed() << std::endl;
    std::cerr << " - loading filter: " << timeLoadFilter.elapsed() << std::endl;
    std::cerr << " - loading sequences (1t): " << timeLoadSeq.elapsed() << std::endl;
    std::cerr << " - building filter (" << config.threads_build << "t): " << timeBuild.elapsed() << std::endl;
    std::cerr << " - saving filter: " << timeSaveFilter.elapsed() << std::endl;
    std::cerr << " - total: " << timeIBF.elapsed() << std::endl;
    std::cerr << std::endl;
}

void print_stats( Stats& stats, const StopClock& timeBuild )
{
    double   elapsed_build = timeBuild.elapsed();
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



bool IBF::build( Config config )
{

    if ( !config.validate() )
        return false;

    StopClock timeIBF;
    timeIBF.start();
    StopClock timeLoadFiles;
    StopClock timeLoadFilter;
    StopClock timeLoadSeq;
    StopClock timeBuild;
    StopClock timeSaveFilter;

    if ( config.verbose )
        std::cerr << config;

    //////////////////////////////

    Stats stats;

    //////////////////////////////

    std::mutex                mtx;
    SafeQueue< Seqs > queue_refs(
        config.n_batches ); // config.n_batches*config.n_refs = max. amount of references in memory

    // Start extra thread for reading the input
    timeLoadSeq.start();
    std::future< void > read_task = parse_ref_seqs(queue_refs, mtx, config, stats);

    // divide ref seqs in fragments, that will be placed in different bins of the ibf
    read_task.get();
    // load new or given ibf
    timeLoadFilter.start();
    std::set< uint64_t > bin_ids;
    load_filter( config, bin_ids, stats );
    timeLoadFilter.stop();
    
    // Start execution threads to add kmers
    
    timeBuild.start();
    uint64_t binid = 0;
    std::vector< std::future< void > > tasks = add_sequences_to_filter(config, binid, queue_refs);

    
    timeLoadSeq.stop();

    for ( auto&& task : tasks )
    {
        task.get();
    }
    timeBuild.stop();
    //////////////////////////////

    // Store filter
    timeSaveFilter.start();
    seqan::store( filter, seqan::toCString( config.output_filter_file ) );
    timeSaveFilter.stop();
    //////////////////////////////

    timeIBF.stop();

    if ( !config.quiet )
    {
        std::cerr << std::endl;
        if ( config.verbose )
        {
            print_time(
                config, timeIBF, timeLoadFiles, timeLoadSeq, timeLoadFiles, timeLoadFilter, timeSaveFilter );
        }
        print_stats( stats, timeBuild );
    }
    return true;
}

} // namespace interleave
