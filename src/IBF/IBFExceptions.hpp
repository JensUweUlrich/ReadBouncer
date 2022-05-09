#pragma once

/*
 * InterleavedBloomFilterException.hpp
 *
 *  Created on: 23.02.2021
 *      Author: Jens-Uwe Ulrich
 */

#include <exception>
#include <string>

// Qt exceptions
#include <QMessageBox>


namespace interleave
{

    class IBFException: public std::exception
	{
		private:
			std::string error_message
			{ };

		public:

			explicit IBFException();

			explicit IBFException(const std::string &msg) :
							error_message(msg)
			{
			}
			virtual ~IBFException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}

	};

	class IBFBuildException: public IBFException
	{

		private:
			std::string error_message
			{ };

		public:

			IBFBuildException() : IBFException()
			{

			}
			IBFBuildException(const std::string &msg) :
							IBFException(msg)
			{
			}
			~IBFBuildException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}

	};

    class IBFClassifyException: public IBFException
	{

		private:
			std::string error_message
			{ };

		public:

			IBFClassifyException() : IBFException()
			{

			}
			IBFClassifyException(const std::string &msg) :
							IBFException(msg)
			{
			}
			~IBFClassifyException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}
	}; 

    class ShortReadException: public IBFClassifyException
	{

		private:
			std::string error_message
			{ };

		public:

			ShortReadException() : IBFClassifyException()
			{

			}
			ShortReadException(const std::string &msg) :
							IBFClassifyException(msg)
			{
			}
			~ShortReadException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}
	}; 

    class CountKmerException: public IBFClassifyException
	{

		private:
			std::string error_message
			{ };

		public:

			CountKmerException() : IBFClassifyException()
			{

			}
			CountKmerException(const std::string &msg) :
							IBFClassifyException(msg)
			{
			}
			~CountKmerException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}
	}; 

	class InvalidConfigException: public IBFBuildException
	{

		private:
			std::string error_message
			{ };

		public:

			InvalidConfigException() : IBFBuildException()
			{

			}
			InvalidConfigException(const std::string &msg) :
							IBFBuildException(msg)
			{
			}
			~InvalidConfigException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}

	};

    class NullFilterException: public IBFBuildException
	{

		private:
			std::string error_message
			{ };

		public:

			NullFilterException() : IBFBuildException()
			{

			}
			NullFilterException(const std::string &msg) :
							IBFBuildException(msg)
			{
			}
			~NullFilterException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}

	};

    class InsertSequenceException: public IBFBuildException
	{

		private:
			std::string error_message
			{ };

		public:

			InsertSequenceException() : IBFBuildException()
			{

			}
			InsertSequenceException(const std::string &msg) :
							IBFBuildException(msg)
			{
			}
			~InsertSequenceException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}

	};

    class StoreFilterException: public IBFBuildException
	{

		private:
			std::string error_message
			{ };

		public:

			StoreFilterException() : IBFBuildException()
			{

			}
			StoreFilterException(const std::string &msg) :
							IBFBuildException(msg)
			{
			}
			~StoreFilterException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}

	};

    class FileParserException: public IBFBuildException
	{

		private:
			std::string error_message
			{ };

		public:

			FileParserException() : IBFBuildException()
			{

			}
			FileParserException(const std::string &msg) :
							IBFBuildException(msg)
			{
			}
			~FileParserException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}

	};

   class MissingReferenceFilesException: public FileParserException
	{

		private:
			std::string error_message
			{ };

		public:

			MissingReferenceFilesException() : FileParserException()
			{

			}
			MissingReferenceFilesException(const std::string &msg) :
							FileParserException(msg)
			{
			}
			~MissingReferenceFilesException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}
	};

    class MissingIBFFileException: public FileParserException
	{

		private:
			std::string error_message
			{ };

		public:

			MissingIBFFileException() : FileParserException()
			{

			}
			MissingIBFFileException(const std::string &msg) :
							FileParserException(msg)
			{
			}
			~MissingIBFFileException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}
	};

    class ParseIBFFileException: public FileParserException
	{

		private:
			std::string error_message
			{ };

		public:

			ParseIBFFileException() : FileParserException()
			{

			}
			ParseIBFFileException(const std::string &msg) :
							FileParserException(msg)
			{
			}
			~ParseIBFFileException() throw ()
			{
			}

			virtual const char* what() const throw ()
			{
				return error_message.c_str();
			}
	};  

}

