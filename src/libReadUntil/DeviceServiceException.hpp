/*
 * DeviceServiceException.hpp
 *
 *  Created on: 14.11.2019
 *      Author: jens
 */

#include "ReadUntilClientException.hpp"
#include <string>

#ifndef LIBREADUNTIL_DEVICESERVICEEXCEPTION_HPP_
#define LIBREADUNTIL_DEVICESERVICEEXCEPTION_HPP_

namespace readuntil
{

	class DeviceServiceException: public ReadUntilClientException
	{

		private:
			std::string error_message
			{ };

		public:

			DeviceServiceException() : ReadUntilClientException()
			{

			}
			DeviceServiceException(const std::string &msg) :
							ReadUntilClientException(msg)
			{
			}
			~DeviceServiceException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}

	};

}
#endif /* LIBREADUNTIL_DEVICESERVICEEXCEPTION_HPP_ */
