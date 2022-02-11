//#include <gtest/gtest.h>
#include <queue>

/*
The above uses both ASSERT_* and EXPECT_* assertions.
 The rule of thumb is to use EXPECT_* when you want the test to continue to reveal more errors after the assertion failure,
and use ASSERT_* when continuing after failure doesn’t make sense.
*/

#include <vector>
#include <regex>
#include <iostream>
 
std::string canonicalpath(const std::string &path)
{
    if (path.length() <= 1)
        return path;
 
    std::string sep = path[0] == '/' ? "/" : "";
 
    std::vector<std::string> entries;
    std::smatch match;
    std::regex re("[^/]+");
    for (auto p = path; std::regex_search(p, match, re); p = match.suffix()) {
        if (match.str() == ".." && !entries.empty()
                && !(sep == "" && entries.back() == ".."))
            entries.pop_back();
        else
            entries.push_back(match.str());
    }
 
    std::string cpath;
    for (auto s: entries) {
        cpath += sep + s;
        sep = "/";
    }
    return cpath;
}

/*
TEST(canonicalTests, relativePath) {
    EXPECT_STREQ(canonicalpath("abc/de/").data(), "abc/de");
    EXPECT_STREQ(canonicalpath("abc/../de").data(), "de");
    EXPECT_STREQ(canonicalpath("../../abc").data(), "../../abc");
    EXPECT_STREQ(canonicalpath("abc/../../../de").data(), "../../de");
    EXPECT_STREQ(canonicalpath("abc/../de/../fgh").data(), "fgh");
}
 
TEST(canonicalTests, absolutePath) {
    EXPECT_STREQ(canonicalpath("/abc/de/").data(), "/abc/de");
    EXPECT_STREQ(canonicalpath("/abc/../de").data(), "/de");
    EXPECT_STREQ(canonicalpath("/../../abc").data(), "/abc");
    EXPECT_STREQ(canonicalpath("/abc/../../../de").data(), "/de");
    EXPECT_STREQ(canonicalpath("/abc/../de/../fgh").data(), "/fgh");
}
 
TEST(canonicalTests, boundaryCase) {
    EXPECT_STREQ(canonicalpath("").data(), "");
    EXPECT_STREQ(canonicalpath("/").data(), "/");
}*/