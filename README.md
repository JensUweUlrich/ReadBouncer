# NanoLiveTk

## Overview
c++ Toolkit making use of Oxford Nanopore's Read Until functionality. Instead of using the Python based Read Until API, NanoLiveTk directly uses Google's remote procedure calls to pull live read data from the MinION device by communicating with the MinKNOW API. 

## Installation

### Compilation From Source
NanoLiveTk has the following dependencies

* [CMake](https://cmake.org/) for building
* [C++ Compiler](https://nuwen.net/mingw.html) like MinGW or GCC
* [Perl](http://strawberryperl.com/releases.html)
* [Go](https://golang.org/dl/)
* [NASM Compiler](https://www.nasm.us/)
* [OpenSSL](https://bintray.com/vszakats/generic/openssl)
