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
#include "Instance.hpp"

#ifndef LIBREADUNTIL_READUNTILCLIENT_HPP_
#define LIBREADUNTIL_READUNTILCLIENT_HPP_

namespace readuntil
{
	class ReadUntilClient
	{

		private:

			std::shared_ptr<::grpc::Channel> channel;
			std::string mk_host
			{ "127.0.0.1" };
			uint16_t mk_port
			{ 8002 };

			ReadUntilClient() = default;
			~ReadUntilClient() = default;
			ReadUntilClient(const ReadUntilClient&) = delete;

		public:

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

			void connect();

			MinKnowService* getMinKnowService(const MinKnowServiceType type);

			ReadUntilClient& operator=(const ReadUntilClient&) = delete;

	};
}
#endif /* LIBREADUNTIL_READUNTILCLIENT_HPP_ */

