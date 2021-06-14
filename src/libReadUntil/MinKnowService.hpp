/*
 * MinKnowServiceFactory.hpp
 *
 * abstract factory class for minknow api services
 *
 *  Created on: 08.11.2019
 *      Author: jens-uwe ulrich
 */

#include <grpcpp/grpcpp.h>

#include "MinKnowServiceType.hpp"

#ifndef LIBREADUNTIL_MINKNOWSERVICE_HPP_
#define LIBREADUNTIL_MINKNOWSERVICE_HPP_

namespace readuntil
{
	class MinKnowService
	{
		protected:
			::grpc::ClientContext context;

		public:
			MinKnowService() = default;
			MinKnowService& operator=(MinKnowService &other)
			{
				return *this;
			}
			~MinKnowService()
			{
				// ensures signal to server, that the stream is done
                context.TryCancel();
			}
	};
}

#endif // define LIBREADUNTIL_MINKNOWSERVICE_HPP_
