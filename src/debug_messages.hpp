/*
 * debug_messages.hpp
 *
 *  Created on: 12.11.2019
 *      Author: jens
 */

#ifndef DEBUG_MESSAGES_HPP_
#define DEBUG_MESSAGES_HPP_

#ifdef DEBUG
	#define DEBUGVAR( os, var ) \
	(os) << "DEBUG: " << __FILE__ << "(" << __LINE__ << ") "\
	<< #var << " = [" << (var) << "]" << std::endl

	#define DEBUGMESSAGE( os, msg ) \
	(os) << "DEBUG: " << __FILE__ << "(" << __LINE__ << ") " \
    << msg << std::endl
#endif // DEBUG
#endif /* DEBUG_MESSAGES_HPP_ */
