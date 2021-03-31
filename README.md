# NanoLIVE

## Table of Contents

- [Overview](#started)
- [Users' Guide](#uguide)
  - [Installation](#install)
    - [Compilation From Source](#compile)
  - [General usage](#general)
    - [Building the database](#ibfbuild)
    - [Classify Query Reads](#classify)
    - [Live Depletion of Nanopore Reads](#deplete)
  - [Use cases](#ucase)
    - [Deplete Host Background Reads](#host-depletion)
  - [Advanced features](#advanced)
  - [Algorithm overview](#algo)
  - [Getting help](#help)
  - [Citing minimap2](#cite)
- [Developers' Guide](#dguide)
- [Limitations](#limit)

## <a name="overview"></a>Overview
C++ based tool for live classification of Nanopore reads (aka adaptive sequencing) on Windows without the need for GPUs. The Toolkit uses Oxford Nanopore's Read Until functionality to unblock reads that match to a given reference sequence database. The database is indexed as Interleaved Bloom Filter for fast classification.
* In a first step the reference sequences were used to build an interleaved bloom filter using the [SeqAn library](https://github.com/seqan/seqan3)
* After Starting a sequencing run NanoLIVE pulls signals from the sequencer via ONT's [MinKNOW API](https://github.com/nanoporetech/minknow_api) by using Google's remote procedure calls ([gRPC](https://grpc.io/))
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

Then just need to clone the repository to your computer, create a <b>build</b> directory within cloned directory and let cmake do the work for you

```
git clone https://github.com/JensUweUlrich/NanoLive.git
cd .\NanoLive\
mkdir build
cd .\build\
cmake.exe  ..\src
cmake.exe --build . --config Release
cmake.exe --build . --config Release --target package
```
The last step creates the <b>NanoLive-0.3.0-win64.exe</b> within the build directory, which is a simple installer for Windows that leads you through the installation process.


### <a name="general"></a>General usage

#### <a name="ibfbuild"></a>Building the database
Before you can use NanoLIVE for adaptive sequencing, the reference database has to be build by using the subcommand <b>ibfbuild</b>. In this step you have to provide a reference sequence file in FASTA format, an output file in which the Interleaved Bloom Filter (IBF) shall be stored and the size of the kmers used to build the IBF.  

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

<b>--fragment-size</b>
The reference sequence is fragmented in subsequences of this length, where every fragment is stored in a separate bin of the interleaved bloom filter. Fragments overlap by 500 nucleotides since the expected length of pulled read information from MinKNOW is between 200 and 400 nucleotides. The fragmentation leads to better classification specificity when using smaller ```kmer-size``` values. By default we recommend a fragment size of 100000 basepairs. But for smaller expected sequencing errors, larger ```kmer-size``` values can be used and thus the size of the fragments can be increased as well.

<b>--kmer-size</b>
For every fragment we compute all kmers of this size, calculate hash values for the kmers and add those to the interleaved bloom filter like described by [Dadi et. al., 2018](https://academic.oup.com/bioinformatics/article/34/17/i766/5093228). For sequencing error rates of around 10% we recommend using k=13. But with ONT's continuous improvements in single read accuracy, higher values are feasible as well.

<b>--filter-size</b>
The optimal filter size is calculated automatically (for human genome it is approx. 4 GigaBytes). But if you are short of disk space or have a large reference sequence set, you can state the maximum filter size here in MegaBytes. <b>Note that this can negatively impact the read classification accuracy.</b>

#### <a name="classify"></a>Classify Query Reads

If you like to test NanoLIVE's read classification with a set of Nanopore reads, you can use the <b>classify</b> subcommand. You only have to provide an interleaved bloom filter (IBF) file and some reads as FASTA or FASTQ file. We than take the first 500 nucleotides from each read and map those against each bin of the IBF. If a certain number of kmers are shared between the 500 nucleotides and the bin, the read is classified as match.

```
classify                classify nanopore reads based on a given IBF file

  -?, -h, --help
  -v, --verbose           Show additional output as to what we are doing.
  -r, --read-file <read-file>
                          File with reads to classify (FASTA or FASTQ format)
  -i, --ibf-file <ibf-file>
                          Interleaved Bloom Filter file
  -s, --significance <probability>
                          significance level for confidence interval of number of errorneous kmers (default is 0.95)
  -e, --error-rate <err>  exepected per read sequencing error rate (default is 0.1)
  -t, --threads <threads> Number of classification threads
```
<b>--significance and --error-rate</b>
For more accurate classification of reads we are calculating the expected number of mutated kmers for each 500 nucleotide start sequence of the read based on the expected sequencing ```error rate```. Than a confidence interval for the mutated kmers is calculated as described by [Blanca et. al., 2021](https://www.biorxiv.org/content/10.1101/2021.01.15.426881v2) and the minimum number of matching kmers is calculated based on the upper bound of the confidence interval. The significance level of the confidence interval can be specified, but is set to 95% by default.

#### <a name="deplete"></a>Live Depletion of Nanopore Reads

When nanopore reads of a certain organism (or even more) are not of interest, these reads can be unblocked with <b>live-deplete</b>. Therefore signals are requested from ONT's MinKNOW software, basecalled and matched against the previously build IBF of the reference sequence set. If a read is classified as match, an unblock request is sent back to the MinKNOW software. This unblocks the corresponding pore of the read and another DNA molecule can enter this pore for sequencing.

```
  live-deplete            Live classification and rejection of nanopore reads

  -?, -h, --help
  -v, --verbose           Show additional output as to what we are doing.
  -d, --device <device>   Device or FlowCell name for live analysis
  -c, --host <host>       IP address on which MinKNOW software runs
  -p, --port <port>       MinKNOW communication port
  -i, --ibf-file <ibf-file>
                          Interleaved Bloom Filter file
  -s, --significance <probability>
                          significance level for confidence interval of number of errorneous kmers (default is 0.95)
  -e, --error-rate <err>  exepected per read sequencing error rate (default is 0.1)
  -w, --weights <weights> Deep Nano Weights (default is 48)
```
<b>--host and --port</b>
This is the IP adress and the TCP/IP port on which the MinKNOW software is hosted. NanoLIVE will exchange data with the MinKNOW software via this communication channel. It is recommended to test the communication before starting the sequencing run.

<b>--device</b>
This is the name of the FlowCell for which we want to do the live depletion. 

<b>--ibf-file, --significance and --error-rate</b>
Before NanoLIVE starts to communicate with MinKNOW, the IBF of the reference sequence set will be loaded. <b>Note that you should always start NanoLIVE before starting the seqencing run in MinKNOW because the IBF has to be loaded first</b>. NanoLIVE will tell you when you can start the sequencing run.
For more accurate classification of reads we are calculating the expected number of mutated kmers for each 500 nucleotide start sequence of the read based on the expected sequencing ```error rate```. Than a confidence interval for the mutated kmers is calculated as described by [Blanca et. al., 2021](https://www.biorxiv.org/content/10.1101/2021.01.15.426881v2) and the minimum number of matching kmers is calculated based on the upper bound of the confidence interval. The significance level of the confidence interval can be specified, but is set to 95% by default. 

<b>--weights</b>
For CPU based real-time basecalling of Nanopore reads, NanoLIVE integrates [DeepNano-blitz](https://github.com/fmfi-compbio/deepnano-blitz). This basecaller uses recurrent neural networks (RNNs) for signal-to-nucleotide translation. There are different sizes of RNNs available, e.g. 48, 56, 64, 80, 96 and 256. In general, the smaller the RNN the faster basecalling is performed. But on the other hand higher RNN weight values provide higher base call accuracy. <b>Note that necessary basecalling speed could only be supported for maximum weights value of 80 on a an Intel Core i7 2,8 GHz processor. Therefore we recommend values smaller than 80 to keep up with sequencing speed.</b>

#### <a name="test"></a>Testing Connection NanoLIVE

Before starting a sequencing run with live depletion, it's recommended to test the connection between NanoLIVE and MinKNOW by using subcommand `connection-test`.

```
connection-test         Test connection to a working MinKNOW instance

  -?, -h, --help
  -v, --verbose           Show additional output as to what we are doing.
  -d, --device <device>   Device or FlowCell name for live analysis (required)
  -c, --host <host>       IP address on which MinKNOW software runs (default: localhost)
  -p, --port <port>       MinKNOW communication port (default: 9501)
  -u, --unblock-all       Unblock all reads
```
<b>--host and --port</b>
This is the IP adress and the TCP/IP port on which the MinKNOW software is hosted. NanoLIVE will exchange data with the MinKNOW software via this communication channel. It is recommended to test the communication before starting the sequencing run.

<b>--device</b>
This is the name of the FlowCell for which we want to do the live depletion. 

<b>--unblock-all</b>
If you have a bulk FAST5 file at hand, you can simulate a Nanopore sequencing run and try out if MinKNOW accepts messages for unblocking pores from NanoLIVE. You can simply replay the sequencing run from the bulk fast5 file and state NanoLIVE to send unblock messages for all reads. If you observe that lots of reads have lengths below 1kb unblocking works as expected.

### <a name="ucase"></a>Use Cases

### <a name="unblockall"></a>Test NanoLIVE-to-MinKNOW interaction
Before using NanoLIVE in a real experiment, we recommend first running a playback experiment to test unblock speed.

1. Download an open access bulk FAST5 file from 
[here](http://s3.amazonaws.com/nanopore-human-wgs/bulkfile/PLSP57501_20170308_FNFAF14035_MN16458_sequencing_run_NOTT_Hum_wh1rs2_60428.fast5). 
This file is 21Gb so make sure you have plenty of space.
2. To configure a run for playback, you need to find and edit a sequencing TOML file. These are typically located in `C:\Program Files\OxfordNanopore\MinKNOW\conf\package\sequencing`. Edit a file such as sequencing_MIN106_DNA.toml and under the entry `[custom_settings]` 
add a field: 
    ```text
    simulation = "/full/path/to/your_bulk.FAST5"
    ```
3. Open a Windows Power Shell with administrator privileges and go to the MinKNOW binary directory, located in `C:\Program Files\OxfordNanopore\MinKNOW\bin`, and call the config editor with the following two commands:
```
.\config_editor.exe --conf application --filename ../conf/app_conf --set data_generation.simulated_device=1
.\config_editor.exe --conf application --filename ../conf/app_conf --set device.simulator_device_count=1
```
4. In the MinKNOW GUI, right click on a sequencing position and select `Reload Scripts` (In some cases you need to reboot Windows). Your version of MinKNOW will now show a simulated device named `MS00000` that will playback the bulkfile rather than live sequencing.
5. Open a Windows Power Shell and go to your working directory where NanoLIVE result files shall be stored. Than call NanoLIVE with subcommand `connection-test` and correct parameters for host, port and device name
```
full\path\to\NanoLIVE\root\directory\bin\NanoLive.exe connection-test --host 127.0.0.1 --port 9501 --device MS00000
```
The output should state that the connection could be successfully established and that you can continue with live-depletion.
5. For testing the unblock functionality you should start a simulation with `--unblock-all` option 
```
full\path\to\NanoLIVE\root\directory\bin\NanoLive.exe connection-test --host 127.0.0.1 --port 9501 --device MS00000 --unblock-all
```
When NanoLIVE says that it successfully established a connection, you can start a sequencing run on the the device, which will playback the run from the bulkfile.
6. Open the read length histogram after 5 minutes and have a look at the read counts plot.
![alt text](images/unblock_all.png "Unblock All Image")
7. Now zoom in to the histogram so that only read counts for read lengths up to 5kb are shown. You should see a peak for read counts between 500b and 1 kb like the one in the figure below.
![alt text](images/unblock_all_5kb.png "Unblock All Image (5kb)")
If that's the case you can go on with testing basecalling and classification

### <a name="host-depletion"></a>Live-Basecalling and read classification

In order to test if read depletion works on your machine, you can start a `live-depletion` playback run with the bulk FAST5 file from the test above. If you already set up the playback functionality, you only need to download the reference sequence of one or more human chromosomes from e.g. the [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39) as FASTA file. In the example below, chromosome 3 reads shall be depleted.

1. Before live-depletion, we need to build an Interleaved Bloom Filter (IBF) for the reference sequence(s) we aim to deplete. 
```
full\path\to\NanoLIVE\root\directory\bin\NanoLive.exe ibfbuild -o path\to\output\directory\hg38p13_chr3.ibf -i path\to\reference\file\hg38p13_chr3.fasta -k 13 -f 100000
```
2. Now you can start live-depletion of chromosome 3 with the following subcommand from you working directory
 ```
full\path\to\NanoLIVE\root\directory\bin\NanoLive.exe live-deplete -i path\to\output\directory\hg38p13_chr3.ibf -d MS00000
```
3. Start a sequencing run on the simulated device as you did above. Open the read length histogram after 15 minutes and have a look at the read counts plot. When you zoom into the region for reads up to 5kb length, you should see a plot like this:
![alt text](images/unblock_all.png "Unblock All Image")

4. After stopping the run, NanoLIVE will provide you with some statistics about the number of classified (unblocked) and unclassified reads, which will be sequenced until the end. You will also see average overall processing times as well as for basecalling and classification. You should aim for overall processing times for classified reads below one second. The average processing time for basecalling and classification should be below 0.02 seconds. Otherwise you will experience 