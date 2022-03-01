#include "IBF.hpp"
#include "gtest/gtest.h"
#include "gmock/gmock.h"


class IBFTest: public ::testing::Test
{
	protected:

		interleave::IBF *ibf;
        
		// prepare the objects for each test.
		void SetUp() override
		{
			ibf = new interleave::IBF();
		}
        // release any resources we allocated in SetUp()
		void TearDown() override
		{
			delete ibf;
		}

    public:



};

TEST_F (IBFTest, ParseRefSeqsTest){

    
}

int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
