/*
 * GuppyBasecaller.hpp
 *
 *  Created on: 11.10.2021
 *      Author: juu
 */

#pragma once

#include <string>
#include <future>
#include <SafeQueue.hpp>
#include <SafeMap.hpp>
#include <filesystem>
#include <ont_read.hpp>
#include <runner.hpp>
#include <DeepNano2.h>
#include <Basecaller.hpp>

using namespace interfaces;

namespace basecall
{

	class DeepNanoBasecaller : public Basecaller
	{
	private:
		std::string weights = "48";
		std::string w_file{};
		std::mutex m;
		uint8_t threads = 1;
		SafeMap<std::string, std::pair<uint8_t, RTPair> > pending{};

		void basecall_m(SafeQueue<RTPair>& basecall_queue,
			SafeQueue<RTPair>& classification_queue,
			SafeMap<uint16_t, uint32_t>& channel_stats,
			Runner& runner);

	public:
		DeepNanoBasecaller(std::filesystem::path& weights_file, uint8_t threads);
		~DeepNanoBasecaller() = default;

		/// Copying of object is not allowed.
		DeepNanoBasecaller(const DeepNanoBasecaller&) = delete;

		/// Copying of object is not allowed.
		DeepNanoBasecaller& operator=(const DeepNanoBasecaller&) = delete;

		void basecall_live_reads(SafeQueue<RTPair>& basecall_queue,
			SafeQueue<RTPair>& classification_queue,
			SafeMap<uint16_t, uint32_t>& channel_stats,
			Runner& runner);
	};

} //namespace

