/*
 * GuppyBasecaller.hpp
 *
 *  Created on: 11.10.2021
 *      Author: juu
 */

#pragma once

#include <string>
#include <sstream>
#include <GuppyCPPClient.hpp>
#include <SafeQueue.hpp>
#include <SafeMap.hpp>
#include <filesystem>
#include <future>
#include <ont_read.hpp>
#include <runner.hpp>
#include <Basecaller.hpp>

using namespace interfaces;

namespace basecall
{

	class GuppyBasecaller : public Basecaller
	{
		private:
			std::unique_ptr<guppy_cpp_client::GuppyCPPClient> client;
			std::unordered_map<uint64_t, RTPair> passed_reads{};
			std::unordered_map<std::string, std::pair<uint8_t, RTPair> > pending{};
			std::mutex callerMutex;
			std::mutex passMutex;
			uint16_t batch_size{ 200 };

			void adaptBatchSize(const int queue_size);

		public:
			GuppyBasecaller(std::string& address, std::string& config);
			~GuppyBasecaller() = default;

			/// Copying of object is not allowed.
			GuppyBasecaller(const GuppyBasecaller&) = delete;

			/// Copying of object is not allowed.
			GuppyBasecaller& operator=(const GuppyBasecaller&) = delete;

			void basecall_live_reads(SafeQueue<RTPair>& basecall_queue,
				SafeQueue<RTPair>& classification_queue,
				SafeMap<uint16_t, uint32_t>& channel_stats,
				Runner& runner);

	}
	;

} //namespace

