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
			reads.emplace_back(interleave::Read(id, fragment));

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
*	@parser	: command line input parameters
*/
void classify_reads(read_classify_parser& parser)
{
	std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
	// initialize depletion and target filters
	
	interleave::IBF DepleteFilter{};
	std::vector<interleave::IBFMeta> DepletionFilters{};
	std::vector<interleave::IBFMeta> TargetFilters{};

	bool deplete = false;
	bool target = false;

	// parse depletion IBF if given as parameter
	if (parser.ibf_deplete_file.length() > 0)
	{
		std::vector<std::string> vector_files = split(parser.ibf_deplete_file, ',');
		for (std::string n : vector_files)
		{
			interleave::IBFMeta filter{};
			filter.name = std::filesystem::path(n).stem().string();
			interleave::IBF tf{};
			interleave::IBFConfig DepleteIBFconfig{};
			try
			{
				DepleteIBFconfig.input_filter_file = n;
				interleave::FilterStats stats = std::move(tf.load_filter(DepleteIBFconfig));
				filter.filter = std::move(tf.getFilter());
				if (parser.verbose)
					interleave::print_load_stats(stats);
				deplete = true;
			}
			catch (interleave::ParseIBFFileException& e)
			{
				nanolive_logger->error("Error parsing depletion IBF using the following parameters");
				nanolive_logger->error("Depletion IBF file                : " + n);
				nanolive_logger->error("Error message : " + std::string(e.what()));
				nanolive_logger->flush();
				throw;
			}

			DepletionFilters.emplace_back(std::move(filter));
		}
	}


/*	if (parser.ibf_deplete_file.length() > 0)
	{
		try
		{
			DepleteIBFconfig.input_filter_file = parser.ibf_deplete_file;
			interleave::FilterStats stats = DepleteFilter.load_filter(DepleteIBFconfig);
			DepletionFilters.emplace_back(DepleteFilter.getFilter());
			if (parser.verbose)
				interleave::print_load_stats(stats);
			deplete = true;
		}
		catch (interleave::ParseIBFFileException& e)
		{
			nanolive_logger->error("Error parsing depletion IBF using the following parameters");
			nanolive_logger->error("Depletion IBF file                : " + parser.ibf_deplete_file);
			nanolive_logger->error("Error message : " + std::string(e.what()));
			nanolive_logger->flush();
			throw;
		}
	}
*/
	
	// parse target IBF if given as parameter
	if (parser.ibf_target_file.length() > 0)
	{
		std::vector<std::string> vector_files = split(parser.ibf_target_file, ',');
		for (std::string n : vector_files)
		{
			interleave::IBFMeta filter{};
			filter.name = std::filesystem::path(n).stem().string();
			interleave::IBF tf{};
			interleave::IBFConfig TargetIBFconfig{};
			try
			{
				TargetIBFconfig.input_filter_file = n;
				interleave::FilterStats stats = std::move(tf.load_filter(TargetIBFconfig));
				filter.filter = std::move(tf.getFilter());
				if (parser.verbose)
					interleave::print_load_stats(stats);
				target = true;
			}
			catch (interleave::ParseIBFFileException& e)
			{
				nanolive_logger->error("Error parsing target IBF using the following parameters");
				nanolive_logger->error("Target IBF file                : " + n);
				nanolive_logger->error("Error message : " + std::string(e.what()));
				nanolive_logger->flush();
				throw;
			}

			TargetFilters.emplace_back(std::move(filter));
		}
	}

	// parse input reads
	//interleave::TReads reads;
	//parse_reads(parser.read_file, reads, parser.preLen);

	// create classification config
	interleave::ClassifyConfig Conf{};
	Conf.strata_filter = -1;
	Conf.significance = parser.kmer_significance;
	Conf.error_rate = parser.error_rate;

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
	if (parser.out_dir.length() > 0)
	{
		for (interleave::IBFMeta f : TargetFilters)
		{
			std::filesystem::path outfile(parser.out_dir);
			outfile /= f.name + ".fasta";
			if (!seqan::open(f.outfile, seqan::toCString(outfile.string())))
			{
				std::cerr << "ERROR: Unable to open the file: " << outfile.string() << std::endl;
				return;
			}
		}
		std::filesystem::path outfile(parser.out_dir);
		outfile /= "unclassified.fasta";
		if (!seqan::open(UnclassifiedOut, seqan::toCString(outfile.string())))
		{
			std::cerr << "ERROR: Unable to open the file: " << outfile.string() << std::endl;
			return;
		}
	}

	seqan::SeqFileIn seqFileIn;
	if (!seqan::open(seqFileIn, seqan::toCString(parser.read_file)))
	{
		std::cerr << "ERROR: Unable to open the file: " << parser.read_file << std::endl;
		return;
	}

	int both = 0;
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
		if (seqan::length(seq) < parser.preLen)
		{
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
			while (i < parser.max_chunks)
			{ 
				uint64_t fragend = (i+1) * parser.preLen;
				uint64_t fragstart = i * parser.preLen;
				// make sure that last fragment ends at last position of the reference sequence
				bool last_frag;
				if (fragend > length(seq))
				{
					fragend = length(seq);
					last_frag = true;
				}
				seqan::Infix< seqan::CharString >::Type fragment = seqan::infix(seq, fragstart, fragend);
				r = interleave::Read(id, fragment);
				if (deplete && target)
				{
					//if (r.classify(DepletionFilters, Conf) > -1 )
					//classified = r.classify(DepletionFilters, Conf) > -1 && r.classify(TargetFilters, Conf) == -1;
					std::pair<uint64_t, uint64_t> p = r.classify(DepletionFilters, TargetFilters, Conf);
					if (p.first > 0)
					{
						if (p.second > 0)
						{
							Conf.error_rate -= 0.02;
							p = r.classify(DepletionFilters, TargetFilters, Conf);
							Conf.error_rate += 0.02;
							//std::vector<interleave::TIbf> t;
							//t.emplace_back(DepletionFilters[0].filter);
							//r.classify(t, Conf);

							if (p.first > 0 && p.second == 0)
							{
								classified = true;
							}
							else
							{
								classified = false;
								//std::cerr << id << "\t Depletion Filter k-mers: " << p.first << "\t" << "Target Filter k-mers: " << p.second << "\t"
								//	<< r.classify(DepletionFilters, Conf) << "\t" << r.classify(TargetFilters, Conf) << std::endl;
								//both++;
							}
						}
						else
							classified = true;
					}
					else
						classified = false;
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
						if (parser.out_dir.length() > 0)
							seqan::writeRecord(TargetFilters[best_filter_index].outfile, id, seq);
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
//			if (both == 15)
//				break;
		}
		catch (std::exception& e)
		{
			failed++;
			std::stringstream estr;
			estr << "Error classifying Read : " << r.getID() << "(Len=" << r.getSeqLength() << ")";
			nanolive_logger->error(estr.str());
			estr.str("");
			estr << "Error message          : " << e.what();
			nanolive_logger->error(estr.str());
			nanolive_logger->flush();
		}
		
		avgClassifyduration += (classifyRead.elapsed() - avgClassifyduration) / readCounter;
		std::chrono::duration< StopClock::Seconds > elapsed = classifyRead.end() - begin;
		if (elapsed.count() > 60.0)
		{
			std::stringstream sstr;
			nanolive_logger->info("------------------------------- Intermediate Results -------------------------------");
			sstr << "Number of classified reads                         :   " << found;
			nanolive_logger->info(sstr.str());
			sstr.str("");
			sstr << "Number of of too short reads (len < " << parser.preLen << ")   :   " << too_short;
			nanolive_logger->info(sstr.str());
			sstr.str("");
			sstr << "Number of all reads                                :   " << readCounter;
			nanolive_logger->info(sstr.str());
			sstr.str("");
			sstr << "Average Processing Time Read Classification        :   " << avgClassifyduration;
			nanolive_logger->info(sstr.str());
			nanolive_logger->info("-----------------------------------------------------------------------------------");
			nanolive_logger->flush();
			begin = classifyRead.end();
		}

	}
	seqan::close(seqFileIn);
	
	std::stringstream sstr;
	std::cout << "------------------------------- Final Results -------------------------------" << std::endl;
	std::cout << "Number of classified reads                         :   " << found << std::endl;
	std::cout << "Number of of too short reads (len < " << parser.preLen << ")           :   " << too_short << std::endl;
	std::cout << "Number of all reads                                :   " << readCounter << std::endl;
	for (interleave::IBFMeta f : TargetFilters)
	{
		std::cout << f.name << "\t : " << f.classified << "\t\t" << ((float) f.classified) / ((float) readCounter) << std::endl;
		if (parser.out_dir.length() > 0)
			seqan::close(f.outfile);
	}
	std::cout << "Average Processing Time Read Classification        :   " << avgClassifyduration << std::endl;
	std::cout << "-----------------------------------------------------------------------------------" << std::endl;
}