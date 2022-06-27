
// ReadUntil library
#include "ReadUntilClient.hpp"
#include "Data.hpp"
// IBF library
#include "IBF.hpp"

// tomel parser
#include "configReader.hpp"

//include google tests library 
#include "gtest/gtest.h"

//classify library to be tested
#include "../main/classify.hpp"

// command line parser
//#include "../main/parser.hpp"


// subcommand related function
#include "../main/ibfbuild.hpp"
#include "../main/parser.hpp"

#include "minknowtest.hpp"



int main(int argc, char** argv)
{
	// parse configuration file 
	// tomlFile = parse_config(argc, argv);
	std::string const tomlFile = parse_config(argc, argv);
	testing::InitGoogleTest(&argc, argv);
	
	return RUN_ALL_TESTS();
}
