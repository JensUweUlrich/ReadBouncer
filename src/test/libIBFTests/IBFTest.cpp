#include "gtest/gtest.h"

// IBF headers
#include "IBF.hpp"
#include "IBFConfig.hpp"

// Testing methods
#include "mockIBF.hpp"
#include "read.hpp"
#include "createfilter.hpp"
#include "ibfconfigtest.hpp"







int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);
	
	return RUN_ALL_TESTS();
}
