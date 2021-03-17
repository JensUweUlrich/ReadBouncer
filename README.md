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


### <a name="general"></a>General usage

#### <a name="ibfbuild"></a>Building the database
Before you can use NanoLIVE for adaptive sequencing, the reference database has to be build by using the subcommand <b>ibfbuild</b>. In this step you have to provide a reference sequence file in FASTA format, an output file in which the Interleaved Bloom Filter (IBF) shall be stored ad the size of the kmers used to build the IBF. 
For better classification sensitivity the reference sequences are fragmented and each fragment represents exactly one bin in the IBF. By default we recommend a fragment size of 100000 basepairs. 
The optimal filter size is calculated automatically (for human genome it is approx. 4 GigaBytes). But if you are short of disk space or have a large reference sequence set, you can state the filter size here in MegaBytes. <b>Note that this can negatively influence the read classification accuracy.</b>

```
Build Interleaved Bloom Filter with given references sequences

OPTIONS, ARGUMENTS:
  ibfbuild                Build Interleaved Bloom Filter with given references sequences

  -?, -h, --help
  -v, --verbose           Show additional output as to what we are doing.
  -o, --output-file <output-file>
                          Output file of Interleaved Bloom Filter
  -i, --input-reference <input-reference>
                          Reference sequence file (fasta format) used to build the IBF; reads matching this reference will be filtered out
  -k, --kmer-size <kmer-size>
                          Kmer size used for building the Interleaved Bloom Filter
  -t, --threads <threads> Number of building threads
  -f, --fragment-size <fragment-size>
                          Length of fragments from the reference that are put in one bin of the IBF
  -s, --filter-size <filter-size>
                          IBF size in MB
```

#### <a name="classify"></a>Classify Query Reads

If you like to test NanoLIVE's read classification with a set of Nanopore reads, you can use the <b>classify</b> subcommand. 

```
classify                classify nanopore reads based on a given IBF file

  -?, -h, --help
  -v, --verbose           Show additional output as to what we are doing.
  -r, --read-file <read-file>
                          File with reads to classify (FASTA or FASTQ format)
  -i, --ibf-file <ibf-file>
                          Interleaved Bloom Filter file
  -s, --significance <probability>
                          significance level for confidence interval of number of errorneous kmers (default is 0.95
  -e, --error-rate <err>  exepected per read sequencing error rate (default is 0.1
  -t, --threads <threads> Number of classification threads
```

