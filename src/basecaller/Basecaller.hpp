/*
 * GuppyBasecaller.hpp
 *
 *  Created on: 15.11.2021
 *      Author: juu
 */

#pragma once

#include <string>
#include <SafeQueue.hpp>
#include <SafeMap.hpp>
#include <ont_read.hpp>
#include <runner.hpp>

using namespace interfaces;

namespace basecall
{

	class Basecaller
	{
		public:
			virtual void basecall_live_reads(SafeQueue<RTPair>& basecall_queue,
											 SafeQueue<RTPair>& classification_queue,
											 SafeMap<uint16_t, uint32_t>& channel_stats,
											 Runner& runner) = 0;
	};

} //namespace

