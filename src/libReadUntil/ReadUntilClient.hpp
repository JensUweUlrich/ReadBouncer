/*
 * ReadUntilClient.hpp
 *
 * Singleton client class for
 *
 *  Created on: 28.10.2019
 *      Author: jens
 */

#include <string>
#include <sstream>
#include <thread>
#include <chrono>
#include <grpcpp/grpcpp.h>
#include "Acquisition.hpp"
#include "Instance.hpp"
#include "Data.hpp"
#include "Device.hpp"
#include "Manager.hpp"

#include "spdlog/spdlog.h"
#include "spdlog/sinks/rotating_file_sink.h"

#ifndef LIBREADUNTIL_READUNTILCLIENT_HPP_
#define LIBREADUNTIL_READUNTILCLIENT_HPP_

namespace readuntil
{
	class ReadUntilClient
	{

		private:

			std::shared_ptr<::grpc::Channel> channel;
			std::string mk_host{ "127.0.0.1" };
			uint16_t mk_port{ 8000 };
			//std::shared_ptr<spdlog::sinks::daily_file_sink_mt> daily_sink = std::make_shared<spdlog::sinks::daily_file_sink_mt>("RUClientLog", 23, 59);
			std::shared_ptr<spdlog::logger> connection_logger;
			ReadUntilClient() = default;
			ReadUntilClient(const ReadUntilClient&) = delete;

		public:

			inline std::shared_ptr<::grpc::Channel> getChannel()
			{
				return channel;
			}

			static ReadUntilClient& getClient()
			{
				static ReadUntilClient instance;
				return instance;
			}

			inline void setHost(const std::string &newHost)
			{
				mk_host = newHost;
			}

			inline void setPort(const uint16_t &newPort)
			{
				mk_port = newPort;
			}

			bool connect();

			MinKnowService* getMinKnowService(const MinKnowServiceType type);

			~ReadUntilClient()
			{
				channel.reset();
			}

			ReadUntilClient& operator=(const ReadUntilClient&) = delete;

	};
}
#endif /* LIBREADUNTIL_READUNTILCLIENT_HPP_ */

