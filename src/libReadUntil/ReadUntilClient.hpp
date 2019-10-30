/*
 * ReadUntilClient.hpp
 *
 *  Created on: 28.10.2019
 *      Author: jens
 */

//#include <minknow>

#ifndef LIBREADUNTIL_READUNTILCLIENT_HPP_
#define LIBREADUNTIL_READUNTILCLIENT_HPP_

class ReadUntilClient
{

	private:
		std::string mk_host{"127.0.0.1"};
		uint8_t mk_port{8002};

	public:
		ReadUntilClient();
		virtual ~ReadUntilClient();
};

#endif /* LIBREADUNTIL_READUNTILCLIENT_HPP_ */
