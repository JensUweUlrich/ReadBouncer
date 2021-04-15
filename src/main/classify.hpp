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


/**
*	classify reads from an input file based on given depletion and/or target filters
*	@parser	: command line input parameters
*/
void classify_reads(read_classify_parser& parser)
{
	std::shared_ptr<spdlog::logger> nanolive_logger = spdlog::get("NanoLiveLog");
	// initialize depletion and target filters
	interleave::IBFConfig DepleteIBFconfig{};
	interleave::IBFConfig TargetIBFconfig{};
	interleave::IBF DepleteFilter{};
	interleave::IBF TargetFilter{};
	std::vector<interleave::TIbf> DepletionFilters{};
	std::vector<interleave::TIbf> TargetFilters{};

	bool deplete = false;
	bool target = false;

	// parse depletion IBF if given as parameter
	if (parser.ibf_deplete_file.length() > 0)
	{
		try
		{
			DepleteIBFconfig.input_filter_file = parser.ibf_deplete_file;
			interleave::FilterStats stats = DepleteFilter.load_filter(DepleteIBFconfig);
			DepletionFilters.emplace_back(DepleteFilter.getFilter());
			interleave::print_stats(stats);
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

	// parse target IBF if given as parameter
	if (parser.ibf_target_file.length() > 0)
	{
		try
		{
			TargetIBFconfig.input_filter_file = parser.ibf_target_file;
			interleave::FilterStats stats = TargetFilter.load_filter(TargetIBFconfig);
			TargetFilters.emplace_back(TargetFilter.getFilter());
			interleave::print_stats(stats);
			target = true;
		}
		catch (interleave::ParseIBFFileException& e)
		{
			nanolive_logger->error("Error parsing target IBF using the following parameters");
			nanolive_logger->error("Target IBF file                : " + parser.ibf_target_file);
			nanolive_logger->error("Error message : " + std::string(e.what()));
			nanolive_logger->flush();
			throw;
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
	seqan::SeqFileOut ClassifiedOut;
	seqan::SeqFileOut UnclassifiedOut;

	// only print classified reads to file if option given via command line
	if (parser.classified_file.length() > 0)
	{
		if (!seqan::open(ClassifiedOut, seqan::toCString(parser.classified_file)))
		{
			std::cerr << "ERROR: Unable to open the file: " << parser.classified_file << std::endl;
			return;
		}
	}
	// only print unclassified reads to file if option given via command line
	if (parser.unclassified_file.length() > 0)
	{
		if (!seqan::open(UnclassifiedOut, seqan::toCString(parser.unclassified_file)))
		{
			std::cerr << "ERROR: Unable to open the file: " << parser.unclassified_file << std::endl;
			return;
		}
	}

	seqan::SeqFileIn seqFileIn;
	if (!seqan::open(seqFileIn, seqan::toCString(parser.read_file)))
	{
		std::cerr << "ERROR: Unable to open the file: " << parser.read_file << std::endl;
		return;
	}

	while (!seqan::atEnd(seqFileIn))
	{
		seqan::CharString id;
		seqan::CharString seq;
		interleave::Read r;
		try
		{
			seqan::readRecord(id, seq, seqFileIn);

			uint64_t fragend = parser.preLen;
			// make sure that last fragment ends at last position of the reference sequence
			if (fragend > length(seq)) fragend = length(seq);
			seqan::Infix< seqan::CharString >::Type fragment = seqan::infix(seq, 0, fragend);
			r = interleave::Read(id, fragment);

		}
		catch (seqan::Exception const& e)
		{
			std::cerr << "ERROR: " << e.what() << " [@" << id << "]" << std::endl;
			continue;
		}


	
		readCounter++;
		StopClock classifyRead;
		classifyRead.start();
		try
		{
			// read length has to be at least the size of the prefix used for read classification
			if (r.getSeqLength() < parser.preLen)
			{
				too_short++;
				continue;
			}

			// only classify if read is in depletion filter but NOT in target filter
			if (deplete && target)
			{
				if (r.classify(DepletionFilters, Conf) && !r.classify(TargetFilters, Conf))
				{
					found++;
					if (parser.classified_file.length() > 0)
					{
						seqan::writeRecord(ClassifiedOut, id, seq);
					}
				}
				else if (parser.unclassified_file.length() > 0)
				{
					seqan::writeRecord(UnclassifiedOut, id, seq);
				}
			}
			// only classify if read is in depletion filter
			else if (deplete)
			{
				if (r.classify(DepletionFilters, Conf))
				{
					found++;
					if (parser.classified_file.length() > 0)
					{
						seqan::writeRecord(ClassifiedOut, id, seq);
					}
				}
				else if (parser.unclassified_file.length() > 0)
				{
					seqan::writeRecord(UnclassifiedOut, id, seq);
				}
			}
			// only classify if read is in target filter
			else
			{
				if (r.classify(TargetFilters, Conf))
				{
					found++;
					if (parser.classified_file.length() > 0)
					{
						seqan::writeRecord(ClassifiedOut, id, seq);
					}
				}
				else if (parser.unclassified_file.length() > 0)
				{
					seqan::writeRecord(UnclassifiedOut, id, seq);
				}
			}
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
		classifyRead.stop();
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
	if (parser.classified_file.length() > 0)
		seqan::close(ClassifiedOut);
	if (parser.unclassified_file.length() > 0)
		seqan::close(UnclassifiedOut);
	std::stringstream sstr;
	std::cout << "------------------------------- Final Results -------------------------------" << std::endl;
	std::cout << "Number of classified reads                         :   " << found << std::endl;
	std::cout << "Number of of too short reads (len < " << parser.preLen << ")           :   " << too_short << std::endl;
	std::cout << "Number of all reads                                :   " << readCounter << std::endl;
	std::cout << "Average Processing Time Read Classification        :   " << avgClassifyduration << std::endl;
	std::cout << "-----------------------------------------------------------------------------------" << std::endl;
}