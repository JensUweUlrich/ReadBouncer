/*
 * classify.hpp
 *
 *  Created on: 31.03.2021
 *      Author: jens-uwe.ulrich
 */

#include <QMessageBox>

/**
 * Write classified reads to file: reference_name + _classified.fasta
 * @param  config ConfigReader object
 * @param  filter IBF
 * @return std::ofstream of the output file
 */

std::ofstream writeClassifed (ConfigReader config, interleave::IBFMeta filter){

    std::filesystem::path outfile(config.output_dir);
    outfile /= filter.name + "_classfied.fasta";
    std::ofstream outf;
    outf.open(outfile, std::ios::out);

    return outf;

}

/**
 * Write unclassified reads to file unclassified.fasta
 * @param  config ConfigReader object
 * @return path to the output file
 */

std::filesystem::path writeUnclassifed (ConfigReader config){

    std::filesystem::path outfile(config.output_dir);
    outfile /= "unclassified.fasta";

    return outfile;

}

/**
 * Classify reads against deplete and target filters
 * @param  DepletionFilters vector of one or more than one depletion interleaved bloom filter(s)
 * @param  TargetFilters vector of one or more than one target interleaved bloom filter(s)
 * @param  Conf classification config object
 * @param  r read's object (sequence and id)
 * @param  best_filter_index index of best matched filter
 * @param  classified bool classification result
 * @param  targetFastas vector of output files of classified reads
 * @param  id read's id
 * @param  seq read's sequence
 * @param  UnclassifiedOut output file of unclassifed reads
 * @return true if read was classified, false otherwise
 */

bool classify_deplete_target(std::vector<interleave::IBFMeta>& DepletionFilters, std::vector<interleave::IBFMeta>& TargetFilters,
                             interleave::ClassifyConfig& Conf, interleave::Read& r, int best_filter_index,
                            std::vector< std::ofstream>& targetFastas, seqan::CharString& id, seqan::CharString& seq, seqan::SeqFileOut &UnclassifiedOut){

    std::pair<uint64_t, uint64_t> p = r.classify(TargetFilters, DepletionFilters, Conf);
    if (p.first > 0)
    {
        if (p.second > 0)
        {
            Conf.error_rate -= 0.02;
            p = r.classify(TargetFilters, DepletionFilters, Conf);
            Conf.error_rate += 0.02;

            if (p.first > 0 && p.second > 0)
            {
                return false;
            }
            else if (p.first > 0)
            {
                best_filter_index = r.classify(TargetFilters, Conf);
                if (best_filter_index != -1)
                {
                    TargetFilters[best_filter_index].classified += 1;
                    targetFastas[best_filter_index] << ">" << id << std::endl;
                    targetFastas[best_filter_index] << seq << std::endl;
                    return true;
                }
            }
            else
            {
                return false;
            }
        }
        else
        {
           best_filter_index = r.classify(TargetFilters, Conf);
           if (best_filter_index != -1)
           {
               TargetFilters[best_filter_index].classified += 1;

               targetFastas[best_filter_index] << ">" << id << std::endl;
               targetFastas[best_filter_index] << seq << std::endl;
               return true;
            }
        }
    }
    else
    {
        return false;
    }

    return false;

}


// Get fragment start
uint64_t fragment_start(int chunk_length, uint8_t i){

    return (i) * chunk_length;
}

// Get fragment end
uint64_t fragment_end(int chunk_length, uint8_t i){

    return (i+1) * chunk_length;
}

// Get classification results for testing
struct ClassificationResults{

    uint64_t found = 0;
    uint16_t failed = 0;
    uint64_t too_short = 0;
    uint64_t readCounter = 0;

}ClassificationResults_;


/**
*	classify reads from an input file based on given depletion and/or target filters
*	@parser	: toml input parameters
*/

void classify_reads(ConfigReader& config, std::vector<interleave::IBFMeta> DepletionFilters, std::vector<interleave::IBFMeta> TargetFilters)
{
    std::shared_ptr<spdlog::logger> readbouncer_logger = spdlog::get("ReadBouncerLog");
    // create classification config
    interleave::ClassifyConfig Conf{};

    bool deplete = false;
    bool target = false;

    if(DepletionFilters.size() >= 1){

        deplete = true;

    } if(TargetFilters.size() >= 1){

        target = true;

    } if(!deplete && !target) {

        QMessageBox::critical(NULL , "Error", "No depletion or target filters have been provided for classification! ");
                return;
    }

    //std::cout<< "Size of depletion filters: "<< deplete << '\n';
    //std::cout<< "Size of target filters: "<< target << '\n';

    for (std::filesystem::path read_file : config.IBF_Parsed.read_files)
    {

        Conf.strata_filter = -1;
        //Conf.significance = params.kmer_significance;
        Conf.significance = 0.95;// default value
        Conf.error_rate =config.IBF_Parsed.error_rate;

        uint64_t found = 0;
        uint16_t failed = 0;
        uint64_t too_short = 0;
        uint64_t readCounter = 0;
        StopClock::Seconds avgClassifyduration = 0;
        // start stop clock
        StopClock::TimePoint begin = StopClock::Clock::now();


        // initialize classification output files
        seqan::SeqFileOut UnclassifiedOut;

        // only print classified reads to file if option given via command line
        /*std::vector< std::ofstream> targetFastas{};
        for (interleave::IBFMeta f : TargetFilters)
        {
            std::ofstream outf = writeClassifed(config, f);
            targetFastas.emplace_back(std::move(outf));
        }
        */

        std::vector< std::ofstream> targetFastas{};
        for (interleave::IBFMeta f : TargetFilters)
        {
            std::filesystem::path outfile(config.output_dir);
            outfile /= f.name + ".fasta";
            std::ofstream outf;
            outf.open(outfile, std::ios::out);

            targetFastas.emplace_back(std::move(outf));
        }

        std::filesystem::path outfile(config.output_dir);
        outfile /= "unclassified.fasta";
        if (!seqan::open(UnclassifiedOut, seqan::toCString(outfile.string())))
        {
            QMessageBox::critical(NULL , "Unable to open the file: ", QString::fromStdString(outfile.string()));
            return;
        }

        seqan::SeqFileIn seqFileIn;
        if (!seqan::open(seqFileIn, seqan::toCString(read_file.string())))
        {
            QMessageBox::critical(NULL , "Unable to open the file: ", QString::fromStdString(read_file.string()));
            return;
        }

        std::cout << '\n';
        std::cout << "Classification results of: " << read_file.string() << '\n';
        std::cout << '\n';


        while (!seqan::atEnd(seqFileIn))
        {
            interleave::Read r;
            seqan::CharString id;
            seqan::CharString seq;

            try
            {
                seqan::readRecord(id, seq, seqFileIn);
                readCounter++;
            }
            catch (seqan::Exception const& e)
            {
                QString titel = "Error: ";
                titel.push_back(e.what());
                QString msg = "[@ ";
                msg.push_back(QString::fromStdString(toCString(id)));
                msg.push_back("]");
                QMessageBox::critical(NULL , titel, msg);

                continue;
            }

            // read length has to be at least the size of the prefix used for read classification
            if (seqan::length(seq) < config.IBF_Parsed.chunk_length){

                too_short++;
                continue;
            }


            StopClock classifyRead;
            bool classified = false;
            int best_filter_index = -1;
            classifyRead.start();
            try
            {
                // as long as rea
                uint8_t i = 0;
                // try to classify read parser.max_chunks times
                while (i < config.IBF_Parsed.max_chunks)
                {
                    uint64_t fragend = fragment_end(config.IBF_Parsed.chunk_length, i);
                    uint64_t fragstart = fragment_start(config.IBF_Parsed.chunk_length, i);

                    // make sure that last fragment ends at last position of the reference sequence
                    if (fragend > length(seq))
                        fragend = length(seq);

                    seqan::Infix< seqan::CharString >::Type fragment = seqan::infix(seq, fragstart, fragend);
                    seqan::Dna5String fr = (seqan::Dna5String) fragment;
                    r = interleave::Read(id, fr);

                    if (deplete && target)
                    {
                        classified = classify_deplete_target(DepletionFilters, TargetFilters, Conf, r, best_filter_index,targetFastas, id, seq, UnclassifiedOut);

                    }
                    else if (deplete)
                        classified = r.classify(DepletionFilters, Conf) > -1;
                    else
                    {
                        best_filter_index = r.classify(TargetFilters, Conf);
                        if (best_filter_index != -1)
                        {
                            classified = true;
                            TargetFilters[best_filter_index].classified += 1;
                            targetFastas[best_filter_index] << ">" << id << std::endl;
                            targetFastas[best_filter_index] << seq << std::endl;
                        }
                    }
                    if (classified)// || last_frag)
                    {
                        found++;
                        break;
                    }
                    i++;
                }
                if (!classified)
                    seqan::writeRecord(UnclassifiedOut, id, (seqan::Dna5String)seq);

                classifyRead.stop();

            }
            catch (std::exception& e)
            {
                failed++;
                std::stringstream estr;
                estr << "Error classifying Read : " << r.id << "(Len=" << seqan::length(r.sequence) << ")";
                        readbouncer_logger->error(estr.str());
                estr.str("");
                estr << "Error message          : " << e.what();
                        readbouncer_logger->error(estr.str());
                        readbouncer_logger->flush();

                        QMessageBox::critical(NULL , "Error", QString::fromStdString(estr.str()));
            }

            avgClassifyduration += (classifyRead.elapsed() - avgClassifyduration) / readCounter;
            std::chrono::duration< StopClock::Seconds > elapsed = classifyRead.end() - begin;
            if (elapsed.count() > 60.0)
            {
                std::stringstream sstr;
                readbouncer_logger->info("------------------------------- Intermediate Results -------------------------------");
                sstr << "Number of classified reads                         :   " << found;
                readbouncer_logger->info(sstr.str());
                sstr.str("");

                sstr << "Number of of too short reads (len < " << config.IBF_Parsed.chunk_length << ")   :   " << too_short;
                readbouncer_logger->info(sstr.str());
                sstr.str("");
                sstr << "Number of all reads                                :   " << readCounter;
                readbouncer_logger->info(sstr.str());
                sstr.str("");
                sstr << "Average Processing Time Read Classification        :   " << avgClassifyduration;
                readbouncer_logger->info(sstr.str());
                readbouncer_logger->info("-----------------------------------------------------------------------------------");
                readbouncer_logger->flush();
                begin = classifyRead.end();
            }

        }

        for (int i = 0; i < targetFastas.size(); ++i)
            targetFastas[i].close();

        //seqan::close(seqFileIn);
        seqan::close(UnclassifiedOut);

        std::stringstream sstr;
        std::cout << "------------------------------- Final Results -------------------------------" << std::endl;
        std::cout << "Number of classified reads                         :   " << found << std::endl;
        std::cout << "Number of of too short reads (len < " << config.IBF_Parsed.chunk_length << ")           :   " << too_short << std::endl;
        std::cout << "Number of all reads                                :   " << readCounter << std::endl;

        for (interleave::IBFMeta f : TargetFilters)
        {

            std::cout << f.name << "\t : " << f.classified << "\t\t" << ((float) f.classified) / ((float) readCounter) << std::endl;


            //seqan::close(f.outfile);

        }
        std::cout << "Average Processing Time Read Classification        :   " << avgClassifyduration << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;

        // Testing struct
        ClassificationResults_.found = found;
        ClassificationResults_.failed = failed;
        ClassificationResults_.too_short = too_short;
        ClassificationResults_.readCounter = readCounter;

        found = 0;
        failed = 0;
        too_short = 0;
        readCounter = 0;
        avgClassifyduration = 0;

    }
}
