Test GUI with Qt6 for NanoLiveTk 

c++ Toolkit making use of Oxford Nanopore's Read Until functionality. Instead of using the Read Until API NanoLiveTk directly uses grpc calls to pull live read data by communicating with the MinKNOW API. 

Installation

From Source

Requirements

* CMake
* C++ Compiler (e.g. MinGW https://nuwen.net/mingw.html)
* Perl (e.g. http://strawberryperl.com/releases.html)
* Go (https://golang.org/dl/)
* NASM Compiler (https://www.nasm.us/)
* OpenSSL (https://bintray.com/vszakats/generic/openssl)
* Qt >= 5 (https://www.qt.io/download)

-------------------------------------------------------------------------------------------------------------------------
Requirements for Qt: 

* CDB Debugger 
* MSVC 2019  x86_64
* QMake 
* CMake tool for QT 6.0.1
* UI Designer 
-------------------------------------------------------------------------------------------------------------------------

Just for test (build only with Release) -- Build with PowerShell and Qt_designer: 

```
cmake.exe  -DCMAKE_PREFIX_PATH="c:/QT/QT5.12.10/6.1.0/msvc2019_64/lib/cmake/" ../[$PREFIX] CMAKE_BINARY_DIR
cmake.exe --build . --target ALL_BUILD

```

Note: avoid using `cmake.exe --build . --target all` 

