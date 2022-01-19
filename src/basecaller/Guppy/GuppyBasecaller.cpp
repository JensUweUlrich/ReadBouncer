/*
 * GuppyBasecaller.cpp
 *
 *  Created on: 11.10.2021
 *      Author: juu
 */

#include "GuppyBasecaller.hpp"
#include <iostream>

namespace basecall
{

	

	GuppyBasecaller::GuppyBasecaller(std::string& address, std::string& config) :
		client(new guppy_cpp_client::GuppyCPPClient(std::ref(address), std::ref(config)))
	{
		guppy_cpp_client::Priority p = guppy_cpp_client::Priority::high;
		client->set_priority(p);
		client->set_move_and_trace_enabled(false);
		client->set_state_data_enabled(false);
		client->set_max_reads_queued(1000);
		client->set_reconnect_timeout(30);
		if (client->connect() != guppy_cpp_client::Result::success)
		{
			throw BasecallerException(client->get_error_message());
		}
			
	}


	void GuppyBasecaller::adaptBatchSize(const int queue_size)
	{
		if (queue_size > 0)
		{
			batch_size += queue_size;
		}
		else
		{
			batch_size *= 0.8;
		}
	}

	

	

	

	void GuppyBasecaller::basecall_live_reads(SafeQueue<RTPair>& basecall_queue,
		SafeQueue<RTPair>& classification_queue,
		SafeMap<uint16_t, uint32_t>& channel_stats,
		Runner& runner)
	{
	
		std::unordered_map<uint64_t, RTPair>::iterator pass_it{};
		std::unordered_map < std::string, std::pair<uint8_t, RTPair>>::iterator pend_it{};
		int d = 0;
		while (true)
		{
			if (!basecall_queue.empty())
			{
				StopClock bc_clock;
				bc_clock.start();
				std::queue<RTPair> batch = basecall_queue.pop_elements(batch_size);
				uint16_t success = 0;
				while (!batch.empty())
				{
					
					RTPair rp = batch.front();
					batch.pop();
					// if we see data from this read for the first time
					// basecall signals and push to classification queue

					rp.second.timeBasecallRead.start();

					std::vector<int16_t> sig{};
					for (float f : rp.first.raw_signals)
						sig.push_back(static_cast<int16_t>(f));

					guppy_cpp_client::GuppyRead guppy_read{ rp.first.readTag, rp.first.id, std::move(sig) };
					guppy_cpp_client::PassReadResult result = client->pass_read(guppy_read);
					if (result == guppy_cpp_client::PassReadResult::success)
						++success;
					else if (result == guppy_cpp_client::PassReadResult::queue_full)
					{
						std::cerr << "Error: Basecalling queue full : " << rp.first << std::endl;
					}
					else if (result == guppy_cpp_client::PassReadResult::bad_read)
					{
						std::cerr << "Error: Bad read : " << rp.first << std::endl;
					}
					else
					{
						std::cerr << "Error: Client not accepting read : " << rp.first << std::endl;
					}
					passed_reads.insert(std::make_pair(rp.first.readTag, std::move(rp)));
				}
				bc_clock.stop();
				//std::cout << "Passed " << success << " reads to the basecall server in " << bc_clock.elapsed() << " seconds." << std::endl;
				uint16_t done = 0;
				StopClock get_clock;
				get_clock.start();
			//	std::queue<RTPair> called{};
				StopClock::Seconds getCompletedAll{0.0};
				StopClock::Seconds findEraseAll{ 0.0 };
				StopClock::Seconds newReadAll{ 0.0 };
				StopClock::Seconds oldReadAll{ 0.0 };
				StopClock::Seconds findPendAll{ 0.0 }; 
				int newreads = 0;
				int oldreads = 0;
				int numWaiting = 0;
				while (done < success)
				{
					StopClock getCompleted;
					getCompleted.start();
					std::vector<guppy_cpp_client::CalledRead> completed_reads{};
					client->get_completed_reads(completed_reads);
					getCompleted.stop();
					getCompletedAll += getCompleted.elapsed();

					if (completed_reads.empty())
					{
						numWaiting++;
					//	std::this_thread::sleep_for(std::chrono::milliseconds(5));
						continue;
					}

					
//					StopClock forLoop;
//					forLoop.start();
					for (guppy_cpp_client::CalledRead called_read : completed_reads)
					{
						StopClock findErase;
						findErase.start();
						pass_it = passed_reads.find(called_read.read_tag);
						RTPair rp = std::move((*pass_it).second);
						rp.second.timeBasecallRead.stop();
						passed_reads.erase(pass_it);
						done++;
						findErase.stop();
						findEraseAll += findErase.elapsed();
						if (called_read.sequence.size() < 1)
						{
							pending.insert({ rp.first.id , std::make_pair(0, std::move(rp)) });
							continue;
						}
						

						// if we see data from this read for the first time
						// basecall signals and push to classification queue
						StopClock findPend;
						findPend.start();
						pend_it = pending.find(called_read.read_id);
						findPend.stop();
						findPendAll += findPend.elapsed();
						if (pend_it == pending.end())
						{
							newreads++;
							StopClock newRead;
							newRead.start();
							rp.first.sequence = called_read.sequence;

							// create some active channel stats
							uint32_t chNr = rp.first.channelNr;
							channel_stats.assign(std::pair(chNr, channel_stats[chNr] + 1));

							// only classify reads with length >= 300
							// reads with length < 300 will wait for the next read chunks and combine them
							if (called_read.seqlen < 200)
							{
								//stop clock for overall processing time of the read
								// add readid and read entry to pending map
								pending.insert({ called_read.read_id , std::make_pair(0, std::move(rp)) });
							}
							else
							{
								classification_queue.push(std::move(rp));
							}
							newRead.stop();
							newReadAll += newRead.elapsed();
						}
						else
						{
							oldreads++;
							StopClock oldRead;
							oldRead.start();
							// if read chunk was too short
							// concatenate sequence data from actual signals with former basecalled sequence data
							// use stop clock from former tries again
							TimeMeasures m = (*pend_it).second.second.second;

							// add elapsed basecall time of first chunks to second chunk
							rp.second.timeBasecallRead.decrementStart(m.timeBasecallRead.runtime());
							// add elapsed overall time of first chunks to second chunk
							rp.second.timeCompleteRead.setBegin(m.timeCompleteRead.begin());

							std::stringstream sstr;
							sstr << (*pend_it).second.second.first.sequence << called_read.sequence;
							rp.first.sequence = sstr.str();
							// push prolonged sequence to classification queue
							// longer read may be classified better
							if (rp.first.sequence.size() < 200)
							{
								// add readid and read entry to pending map
								pending.insert({ called_read.read_id , std::make_pair(0, std::move(rp)) });
							}
							else
							{
								classification_queue.push(std::move(rp));
								pending.erase(pend_it);
							}
							oldRead.stop();
							oldReadAll += oldRead.elapsed();
						}
					}
	//				forLoop.stop();
	//				std::cout << "For Loop " << done << " reads in " << forLoop.elapsed() << " seconds." << std::endl;
				}
		//		classification_queue.push_elements(called);
		//		get_clock.stop();
		//		std::cout << "Get completed " << done << " reads in " << getCompletedAll << " seconds." << std::endl;
		//		std::cout << "Num Waiting " << numWaiting << " call." << std::endl;
		//		std::cout << "Find & Erase Passed " << done << " reads in " << findEraseAll << " seconds." << std::endl;
		//		std::cout << "Find Pending " << done << " reads in " << findPendAll << " seconds." << std::endl;
		//		std::cout << "New Read " << newreads << " reads in " << newReadAll << " seconds." << std::endl;
		//		std::cout << "Old Reads " << oldreads << " reads in " << oldReadAll << " seconds." << std::endl;
		//		std::cout << "Called " << done << " reads in " << get_clock.elapsed() << " seconds." << std::endl;
				//adaptBatchSize(basecall_queue.size());

			}

			

			if (!runner.isRunning)
				break;
		}
		
		client->disconnect();

		
	}

}

