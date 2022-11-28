#include "gtest/gtest.h"

// IBF headers
#include "IBF.hpp"
#include "IBFConfig.hpp"

// Testing methods
#include "mockIBF.hpp"
#include "read.hpp"
#include "createfilter.hpp"
#include "ibfconfigtest.hpp"



/*
IBF.size() = 10
add_seq_to_filter(IBF, "ACGT")
FILTER_stat;
std::string toClassify{ "ACGT" };
toClassify.find(IBF);
if (!notFound) {

	gTest::Fail();
}
*/

/*
std::future<void> get(), asycn nachschauen 
*/



int main(int argc, char** argv)
{
	std::cout << "test1" << std::endl;
	testing::InitGoogleTest(&argc, argv);
	
	return RUN_ALL_TESTS();
}
