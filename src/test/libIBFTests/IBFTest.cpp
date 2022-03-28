#include "IBF.hpp"
#include "gtest/gtest.h"

#include "mockIBF.hpp"
#include "read.hpp"
#include "createfilter.hpp"
#include "IBFConfig.hpp"




int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);
	
	return RUN_ALL_TESTS();
}
