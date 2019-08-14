cmake_minimum_required(VERSION 3.4)

project(googletest-download NONE)

include(ExternalProject)
ExternalProject_Add(googletest
  GIT_REPOSITORY    https://github.com/google/googletest.git
  GIT_TAG           master
  SOURCE_DIR        "/home/jens/Development/MinHashContainment/build/test/googletest-src"
  BINARY_DIR        "/home/jens/Development/MinHashContainment/build/test/googletest-build"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)


