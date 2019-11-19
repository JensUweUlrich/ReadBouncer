/*
 * ReadUntilClientException.hpp
 *
 *  Created on: 11.11.2019
 *      Author: jens
 */

#include <exception>
#include <string>

#ifndef LIBREADUNTIL_READUNTILCLIENTEXCEPTION_HPP_
#define LIBREADUNTIL_READUNTILCLIENTEXCEPTION_HPP_

namespace readuntil
{

	class ReadUntilClientException: public std::exception
	{
		private:
			std::string error_message
			{ };

		public:

			explicit ReadUntilClientException();

			explicit ReadUntilClientException(const std::string &msg) :
							error_message(msg)
			{
			}
			virtual ~ReadUntilClientException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}

	};

} // end namespace readuntil

#endif /* LIBREADUNTIL_READUNTILCLIENTEXCEPTION_HPP_ */
