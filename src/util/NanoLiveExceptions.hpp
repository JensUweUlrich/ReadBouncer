#pragma once

/*
 * Exceptions.hpp
 *
 *  Created on: 23.02.2021
 *      Author: Jens-Uwe Ulrich
 */

#include <exception>
#include <string>


class FileNotFoundException: public std::exception
{
	private:
		std::string error_message{ };
	public:

	explicit FileNotFoundException();

	explicit FileNotFoundException(const std::string &msg) :
						error_message(msg)
	{
	}
	virtual ~FileNotFoundException() throw ()
	{
	}

	virtual const char* what() const throw ()
	{
		return error_message.c_str();
	}
};



