# NanoLIVE

## Overview
C++ based tool for live classification of Nanopore reads (aka adaptive sequencing) on Windows without the need for GPUs. The Toolkit uses Oxford Nanopore's Read Until functionality to unblock reads that match to a given reference sequence database. The database is indexed as Interleaved Bloom Filter for fast classification.
* In a first step the reference sequences were used to build an interleaved bloom filter using the [SeqAn library](https://github.com/seqan/seqan3)
* After Starting a sequencing run NanoLIVE pulls signales from the sequencer via ONT's [MinKNOW API](https://github.com/nanoporetech/minknow_api) by using Google's remote procedure calls
* Signals are basecalled in real-time with [DeepNano-blitz](https://github.com/fmfi-compbio/deepnano-blitz)
* Basecalled reads are matched to the reference database and matched reads are unblocked 

## <a name="uguide"></a>Users' Guide

### <a name="install"></a>Installation

The easiest way is to download the provided installer EXE and simply click through the installation process. 

### <a name="compile"></a>Compilation From Source
NanoLiveTk has the following dependencies that should be installed before compiling the source code

* [CMake](https://cmake.org/) for building
* [C++ Compiler](https://visualstudio.microsoft.com/) Currently we only support MS Visual Studio Compiler
* [Perl](http://strawberryperl.com/releases.html)
* [Go](https://golang.org/dl/)
* [NASM Compiler](https://www.nasm.us/)

Then just need to clone the repository to your computer, create a "build" directory within cloned directory and let cmake do the work for you

```
git clone https://github.com/JensUweUlrich/NanoLive.git
cd .\NanoLive\
mkdir build
cd .\build\
cmake.exe  ..\src
cmake.exe --build . --config Release
cmake.exe --build . --config Release --target package
```
The last step creates the "NanoLive-0.3.0-win64.exe" within the build directory, which is a simple installer for Windows that leads you through the installation process.



