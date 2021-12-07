/*
 * MinKnowServiceType.hpp
 *
 * Enum class for MinKnow Services
 *
 *  Created on: 08.11.2019
 *      Author: jens-uwe ulrich
 */

#ifndef LIBREADUNTIL_MINKNOWSERVICETYPE_HPP_
#define LIBREADUNTIL_MINKNOWSERVICETYPE_HPP_

namespace readuntil
{
	enum MinKnowServiceType
	{
		ACQUISITION,
		ANALYSIS_CONFIGURATION,
		DATA,
		DEVICE,
		INSTANCE,
		KEYSTORE,
		LOG,
		MANAGER,
		MINION_DEVICE,
		PROMETHION_DEVICE,
		PROTOCOL,
		RPC_OPTIONS,
		STATISTICS
	};
}

#endif
