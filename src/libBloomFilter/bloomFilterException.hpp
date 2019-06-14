#include <exception>
#include <string>

class BloomFilterException : public std::exception
{
	protected:
		std::string error_message{};

	public:
		BloomFilterException(const std::string &msg) : error_message(msg) {}
		~BloomFilterException(){}


		virtual const char* what() const throw()
		{
			return error_message.c_str();
		}



};
