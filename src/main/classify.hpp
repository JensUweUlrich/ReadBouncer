/*
 * classify.hpp
 *
 *  Created on: 31.03.2021
 *      Author: jens-uwe.ulrich
 */

void parse_reads(std::string const& reads_file,
	interleave::TReads& reads,
	uint64_t prefixLength)
{
	seqan::SeqFileIn seqFileIn;
	if (!seqan::open(seqFileIn, seqan::toCString(reads_file)))
	{
		std::cerr << "ERROR: Unable to open the file: " << reads_file << std::endl;
		return;
	}

	seqan::CharString id;
	seqan::CharString seq;

	while (!seqan::atEnd(seqFileIn))
	{
		try
		{
			seqan::readRecord(id, seq, seqFileIn);

			uint64_t fragend = prefixLength;
			// make sure that last fragment ends at last position of the reference sequence
			if (fragend > length(seq)) fragend = length(seq);
			seqan::Infix< seqan::CharString >::Type fragment = seqan::infix(seq, 0, fragend);
			seqan::Dna5String fr = (seqan::Dna5String) fragment;
			reads.emplace_back(interleave::Read(id, fr));

		}
		catch (seqan::Exception const& e)
		{
			std::cerr << "ERROR: " << e.what() << " [@" << id << "]" << std::endl;
			break;
		}
	}
	seqan::close(seqFileIn);
}

std::vector<std::string> split(const string& s, char delim) {
	std::vector<std::string> result;
	std::stringstream ss(s);
	string item;

	while (getline(ss, item, delim)) {
		result.push_back(item);
	}

	return result;
}

/**
*	classify reads from an input file based on given depletion and/or target filters
*	@parser	: toml input parameters
*/
void classify_reads(ConfigReader config, std::vector<interleave::IBFMeta> DepletionFilters, std::vector<interleave::IBFMeta> TargetFilters)
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

            std::cerr<<"[Error] No depletion or target filters have been provided for classification! "<<'\n';
            exit(1);
        }
	//std::cout<< "Size of depletion filters: "<< DepletionFilters.size() << '\n';
	//std::cout<< "Size of target filters: "<< TargetFilters.size() << '\n';

	for (std::filesystem::path read_file : config.IBF_Parsed.read_files)
	{

		// parse input reads
		//interleave::TReads reads;
		//parse_reads(parser.read_file, reads, parser.preLen);

		Conf.strata_filter = -1;
		//Conf.significance = params.kmer_significance;
		Conf.significance = 0.95;
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
		std::vector< std::ofstream> targetFastas{};
		for (interleave::IBFMeta f : TargetFilters)
		{
			std::filesystem::path outfile(config.output_dir);
			outfile /= f.name + ".fasta";
			std::ofstream outf;
			outf.open(outfile, std::ios::out);

			/*seqan::SeqFileOut out;
			if (!seqan::open(out, seqan::toCString(outfile.string())))
			{
				std::cerr << "ERROR: Unable to open the file: " << outfile.string() << std::endl;
				return;
			}*/
			targetFastas.emplace_back(std::move(outf));
		}

		std::filesystem::path outfile(config.output_dir);
		outfile /= "unclassified.fasta";
		if (!seqan::open(UnclassifiedOut, seqan::toCString(outfile.string())))
		{
			std::cerr << "ERROR: Unable to open the file: " << outfile.string() << std::endl;
			return;
		}
		

		seqan::SeqFileIn seqFileIn;
		if (!seqan::open(seqFileIn, seqan::toCString(read_file.string())))
		{
			std::cerr << "ERROR: Unable to open the file: " << read_file.string() << std::endl;
			return;
		}

		std::cout << '\n';
		std::cout << "Classification results of: " << read_file.string() << '\n';
		std::cout << '\n';

		while (!seqan::atEnd(seqFileIn))
		{
			seqan::CharString id;
			seqan::CharString seq;
			interleave::Read r;
			try
			{
				seqan::readRecord(id, seq, seqFileIn);
				readCounter++;
			}
			catch (seqan::Exception const& e)
			{
				std::cerr << "ERROR: " << e.what() << " [@" << id << "]" << std::endl;
				continue;
			}

			// read length has to be at least the size of the prefix used for read classification
			if (seqan::length(seq) < config.IBF_Parsed.chunk_length)			{
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
					uint64_t fragend = (i+1) * config.IBF_Parsed.chunk_length;
					uint64_t fragstart = i * config.IBF_Parsed.chunk_length;
					// make sure that last fragment ends at last position of the reference sequence
					if (fragend > length(seq)) fragend = length(seq);
					seqan::Infix< seqan::CharString >::Type fragment = seqan::infix(seq, fragstart, fragend);
					seqan::Dna5String fr = (seqan::Dna5String) fragment;
					r = interleave::Read(id, fr);
					if (deplete && target)
					{
						//if (r.classify(DepletionFilters, Conf) > -1 )
						//classified = r.classify(DepletionFilters, Conf) > -1 && r.classify(TargetFilters, Conf) == -1;
						std::pair<uint64_t, uint64_t> p = r.classify(TargetFilters, DepletionFilters, Conf);
						if (p.first > 0)
						{
							if (p.second > 0)
							{
								Conf.error_rate -= 0.02;
								p = r.classify(TargetFilters, DepletionFilters, Conf);
								Conf.error_rate += 0.02;
								//std::vector<interleave::TIbf> t;
								//t.emplace_back(DepletionFilters[0].filter);
								//r.classify(t, Conf);

								if (p.first > 0 && p.second > 0)
								{
									// do nothing
								}
								else if (p.first > 0)
								{
									classified = true;
									best_filter_index = r.classify(TargetFilters, Conf);
									if (best_filter_index != -1)
									{
										TargetFilters[best_filter_index].classified += 1;
										
										//seqan::writeRecord(outfiles[best_filter_index], id, (seqan::Dna5String) seq);
										targetFastas[best_filter_index] << ">" << id << std::endl;
										targetFastas[best_filter_index] << seq << std::endl;
										
									}
								}
								else
								{
									classified = false;
									seqan::writeRecord(UnclassifiedOut, id, (seqan::Dna5String)seq);
									break;
								}
							}
							else
							{
								classified = true;
								best_filter_index = r.classify(TargetFilters, Conf);
								if (best_filter_index != -1)
								{
									TargetFilters[best_filter_index].classified += 1;
									//seqan::writeRecord(outfiles[best_filter_index], id, seq);
									targetFastas[best_filter_index] << ">" << id << std::endl;
									targetFastas[best_filter_index] << seq << std::endl;
								}
							}
						}
						else
						{
							classified = false;
							seqan::writeRecord(UnclassifiedOut, id, (seqan::Dna5String)seq);
						}
					}
					else if (deplete)
						classified = r.classify(DepletionFilters, Conf) > -1;
					else
					{
						classified = true;
						best_filter_index = r.classify(TargetFilters, Conf);
						if (best_filter_index != -1)
						{
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

		found = 0;
		failed = 0;
		too_short = 0;
		readCounter = 0;
		avgClassifyduration = 0;

	}
}

