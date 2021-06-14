/*
 * DataServiceException.hpp
 *
 *  Created on: 19.11.2019
 *      Author: Jens-Uwe Ulrich
 */

#include "ReadUntilClientException.hpp"
#include <string>

#ifndef LIBREADUNTIL_DATASERVICEEXCEPTION_HPP_
#define LIBREADUNTIL_DATASERVICEEXCEPTION_HPP_

namespace readuntil
{

	class DataServiceException: public ReadUntilClientException
	{

		private:
			std::string error_message
			{ };

		public:

			DataServiceException() : ReadUntilClientException()
			{

			}
			DataServiceException(const std::string &msg) :
							ReadUntilClientException(msg)
			{
			}
			~DataServiceException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}

	};

	class FailedSetupMessageException : public DataServiceException
	{

	private:
		std::string error_message
		{ };

	public:

		FailedSetupMessageException() : DataServiceException()
		{

		}
		FailedSetupMessageException(const std::string& msg) :
			DataServiceException(msg)
		{
		}
		~FailedSetupMessageException() throw ()
		{
		}

		virtual const char* what() const throw ()
		{
			return error_message.c_str();
		}

	};

	class FailedActionRequestException : public DataServiceException
	{

	private:
		std::string error_message
		{ };

	public:

		FailedActionRequestException() : DataServiceException()
		{

		}
		FailedActionRequestException(const std::string& msg) :
			DataServiceException(msg)
		{
		}
		~FailedActionRequestException() throw ()
		{
		}

		virtual const char* what() const throw ()
		{
			return error_message.c_str();
		}

	};
}



#endif /* LIBREADUNTIL_DATASERVICEEXCEPTION_HPP_ */
