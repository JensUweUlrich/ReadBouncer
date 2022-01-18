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

	class BasecallerException: public std::exception
	{
		private:
			std::string error_message
			{ };

		public:

			explicit BasecallerException();

			explicit BasecallerException(const std::string &msg) :
							error_message(msg)
			{
			}
			virtual ~BasecallerException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}

	};


} //namespace

