/*
 * GuppyBasecaller.cpp
 *
 *  Created on: 15.11.2021
 *      Author: juu
 */

#include <DeepNanoBasecaller.hpp>
#include <iostream>
#include <sstream>

namespace basecall
{

	DeepNanoBasecaller::DeepNanoBasecaller(std::filesystem::path& weights_file, uint8_t threads)
	{
		this->w_file = weights_file.string();
		this->threads = threads;
	}


	void DeepNanoBasecaller::basecall_m(SafeQueue<RTPair>& basecall_queue,
		SafeQueue<RTPair>& classification_queue,
		SafeMap<uint16_t, uint32_t>& channel_stats,
		Runner& runner)
	{
		m.lock();
		Caller* caller = std::move(create_caller(weights.c_str(), w_file.c_str(), 5, 0.01));
		m.unlock();

		while (true)
		{
			if (!basecall_queue.empty())
			{
				RTPair rp = std::move(basecall_queue.pop());

				// if we see data from this read for the first time
				// basecall signals and push to classification queue
				if (!pending.contains(rp.first.id))
				{
					rp.second.timeBasecallRead.start();
					rp.first.sequence = call_raw_signal(caller, rp.first.raw_signals.data(), rp.first.raw_signals.size());
					rp.second.timeBasecallRead.stop();
					//interleave::Read r = interleave::Read(read.id, sequence, read.channelNr, read.readNr, read.processingTimes);

					// create some active channel stats
					uint32_t chNr = rp.first.channelNr;
					channel_stats.assign(std::pair(chNr, channel_stats[chNr] + 1));

					// only classify reads with length >= 300
					// reads with length < 300 will wait for the next read chunks and combine them
					if (rp.first.sequence.size() < 250)
					{
						//stop clock for overall processing time of the read
						// add readid and read entry to pending map
						pending.insert({ rp.first.id , std::make_pair(0, std::make_pair(std::move(rp.first), rp.second)) });
					}
					else
					{
						classification_queue.push(std::move(rp));
					}

				}
				else
				{
					// if read was not classified after first and second try
					// concatenate sequence data from actual signals with former basecalled sequence data
					// use stop clock from former tries again

					TimeMeasures m = pending[rp.first.id].second.second;
					rp.second.timeBasecallRead.start();
					char* sequence = call_raw_signal(caller, rp.first.raw_signals.data(), rp.first.raw_signals.size());
					rp.second.timeBasecallRead.stop();

					// add elapsed basecall time of first chunks to second chunk
					rp.second.timeBasecallRead.decrementStart(m.timeBasecallRead.runtime());
					// add elapsed overall time of first chunks to second chunk
					rp.second.timeCompleteRead.setBegin(m.timeCompleteRead.begin());

					std::stringstream sstr;
					sstr << pending[rp.first.id].second.first.sequence << sequence;
					rp.first.sequence = sstr.str();
					// push prolonged sequence to classification queue
					// longer read may be classified better
					//interleave::Read r = interleave::Read(read.id, sstr.str(), read.channelNr, read.readNr, read.processingTimes);
					if (rp.first.sequence.size() < 250)
					{
						// add readid and read entry to pending map
						//read.processingTimes.timeCompleteRead.stop();
						//r.setProcessingTimes(read.processingTimes);
						pending.insert({ rp.first.id , std::make_pair(0, std::make_pair(std::move(rp.first), rp.second)) });
					}
					else
					{
						classification_queue.push(std::move(rp));
						pending.erase(rp.first.id);
					}
				}

			}

			if (!runner.isRunning)
				break;

		}
	}


	void DeepNanoBasecaller::basecall_live_reads(SafeQueue<RTPair>& basecall_queue,
												 SafeQueue<RTPair>& classification_queue,
												 SafeMap<uint16_t, uint32_t>& channel_stats,
												 Runner& runner)
	{
		std::vector< std::future< void > > tasks;
		for (uint8_t t = 0; t < threads; ++t)
		{
			tasks.emplace_back(std::async(std::launch::async, &DeepNanoBasecaller::basecall_m, this, std::ref(basecall_queue),
				std::ref(classification_queue), std::ref(channel_stats), std::ref(runner)));
		}

		for (auto& task : tasks)
		{
			task.get();
		}
	}

}

