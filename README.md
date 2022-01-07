# ReadBouncer (Toml feature) 

## Table of Contents

- [Overview](#started)
- [Users' Guide](#uguide)
  - [Installation](#install)
    - [Compilation From Source](#compile)
      - [Compilation on Windows](#wincompile)
      - [Compilation on Linux](#linuxcompile)
  - [General usage](#general)
    - [Building the database](#ibfbuild)
    - [Classify Query Reads](#classify)
    - [Live Depletion of Nanopore Reads](#deplete)
  - [Use cases](#ucase)
    - [Classify already sequenced reads](#classifyreads)
    - [Test unblocking all reads](#unblockall)
    - [Deplete Host Background Reads](#host-depletion)

## <a name="overview"></a>Overview
C++ based tool for live classification of Nanopore reads (aka adaptive sampling) on Windows and Linux without the need for GPUs. The Toolkit uses Oxford Nanopore's Read Until functionality to unblock reads that match to a given reference sequence database. The database is indexed as Interleaved Bloom Filter for fast classification.
* In a first step the reference sequences are used to build an Interleaved Bloom Filter using the [SeqAn library](https://github.com/seqan/seqan3)
* Interleaved Bloom Filters can be used as depletion or enrichment target
* After Starting a sequencing run, ReadBouncer receives raw current signals from the sequencer via ONT's [MinKNOW API](https://github.com/nanoporetech/minknow_api) using Google's remote procedure calls ([gRPC](https://grpc.io/))
* Signals are basecalled in real-time with [DeepNano-blitz](https://github.com/fmfi-compbio/deepnano-blitz)
* Basecalled reads are matched to the reference sequences in the Interleaved Bloom Filters and unblock/discard messages are sent back to the sequencer according to the read classification 


## <a name="uguide"></a>Users' Guide

### <a name="install"></a>Installation

The easiest way is to download the provided installer files for [Windows](https://owncloud.hpi.de/s/v7D7OoFVnmVymzM) or [Linux](https://owncloud.hpi.de/s/RBpVMqZdGEG0j7r) and simply click through the installation process. 

### <a name="compile"></a>Compilation From Source

ReadBouncer has the following dependencies that should be installed before compiling the source code

* [CMake](https://cmake.org/) for building
* [C++ Compiler](https://visualstudio.microsoft.com/) Currently we support MSVC (on Windows) and GCC (on Linux)
* [Perl](http://strawberryperl.com/releases.html)
* [Go](https://golang.org/dl/)
* [NASM Compiler](https://www.nasm.us/)
* [NSIS](https://nsis.sourceforge.io/Main_Page) Needed for creation of installer executables
* [UUID](https://linux.die.net/man/3/uuid) On Linux machines, `uuid` and `uuid-dev` need to be installed (e.g. `sudo apt install uuid uuid-dev`)
* [TAR](http://gnuwin32.sourceforge.net/packages/gtar.htm) On Windows, make sure tar is installed on the system

Then just need to clone the repository to your computer, create a <b>build</b> directory within cloned directory and let cmake do the work for you

#### <a name="wincompile"></a>Compilation on Windows

Compiling the source code yourself on Windows just requires the following command line calls

```
git clone https://github.com/JensUweUlrich/ReadBouncer.git
cd .\ReadBouncer\
mkdir build
cd .\build\
cmake.exe  ..\src
cmake.exe --build . --config Release
cmake.exe --build . --config Release --target package
```
The last step creates the <b>ReadBouncer-1.0.0-win64.exe</b> within the build directory, which is a simple installer for Windows that leads you through the installation process.

#### <a name="linuxcompile"></a>Compilation on Linux

Compiling the source code Linux works similar to Windows. Just open a terminal, change to your working directory of choice and use the following commands

```
git clone https://github.com/JensUweUlrich/ReadBouncer.git
cd ReadBouncer
mkdir build
cd build
cmake  ../src
cmake --build . --config Release
cmake --build . --config Release --target package
```
The last step creates the <b>ReadBouncer-1.0.0-Linux.sh</b> within the build directory, which is a simple command line installer for Linux that leads you through the installation process. You can also skip the last `cmake` step and just call `sudo make install`, which installs ReadBouncer in your `/usr/local/` directory. 


### <a name="general"></a>General usage

```

usage         = "build", "target", "classify"           #atm only one of those
output_dir    = "path/to/write/output/files/to"         #all generated output files will be stored here
log_directory = "path/to/write/log/files/to"            #all generated log files will be stored here

```

#### <a name="ibfbuild"></a>Building the database
Before you can use ReadBouncer for selective sequencing, the reference database has to be build by using the usage <b>build</b> using the config.toml file. In this step you have to provide the reference sequence(s) as a comma-separated list of FASTA files (target/depletion files) and the size of the kmers used to build the IBF.  

```
Build Interleaved Bloom Filter with given references sequences

[IBF]

kmer_size     = X                       #(unsigned integer with default 13)
fragment_size = X                       #(unsigned integer with default 100000)
threads       = X                       #(unsigned integer with default 3) classification threads
target_files  = "xxx.fasta,xxx.fasta"
deplete_files = "xxx.fasta,xxx.fasta" 

```
<b> kmer_size</b><br>
For every fragment we compute all k-mers of this size, calculate hash values for the k-mers and add those to the Interleaved Bloom Filter like described by [Dadi et. al., 2018](https://academic.oup.com/bioinformatics/article/34/17/i766/5093228). For sequencing error rates of around 10% we recommend using k=13. But with ONT's continuous improvements in single read accuracy, higher values are feasible as well.

<b> fragment_size</b><br>
The reference sequence is fragmented in subsequences of this length, where every fragment is stored in a separate bin of the Interleaved Bloom Filter. Fragments overlap by 500 nucleotides since the expected length of received read information from MinKNOW before classification is between 200 and 400 nucleotides. The fragmentation leads to better classification specificity when using smaller ```kmer-size``` values. By default we recommend a fragment size of 100,000 basepairs. But for smaller expected sequencing errors, larger ```kmer-size``` values can be used and thus the size of the fragments can be increased as well.

#### <a name="classify"></a>Classify Query Reads

If you like to test ReadBouncer's read classification with a set of Nanopore reads, you can use the <b>classify</b> usage. You only have to provide one or more than one Interleaved Bloom Filter (IBF) file or as FASTA files (automatically create ibf file for every given fasta file) and some reads as FASTA or FASTQ file. When providing an IBF file as `depletion-file` ReadBouncer takes the prefix from each read and maps it against each bin of the depletion-IBF. If the number of k-mers that are shared between the prefix and at least one bin of the IBF exceeds a certain threshold, the read is classified as match. The threshold is calculated as described in the publication. If you provide a target IBF file as well, ReadBouncer will not classify reads if the prefix matches the depletion IBF but not the target IBF. 

```
classify                classify nanopore reads based on a given IBF file

[IBF]

kmer_size           = X                        #(unsigned integer with default 13)
fragment_size       = X                        #(unsigned integer with default 100000)
threads             = X                        #(unsigned integer with default 3)
target_files        = "xxx.fasta,xxx.fasta"
deplete_files       = "xxx.fasta,xxx.fasta" 
read_files          = "xxx.fasta/xxx.fastq"    #can be a comma-separated list of fasta or fastq files
exp_seq_error_rate  = 0.1                      #(unsigned float between 0 and 1 default 0.1)
chunk_length        = X                        #(unsigned integer with default 250)
max_chunks          = X                        #(unsigned integer with default 5)


```
<b> chunk_length and max_chunks</b><br>
These parameters decide about the number of nucleotides from the beginning of each read that is used for classification. ReadBouncer tries to classify reads in in chunk-wise manner, when using e.g. 2 chunks with 360 bp length, ReadBouncer takes the first 360 bp of the read for classification and will add the subsequent 360 bp of that read for another classification try if the firts try did not succeed. By default, ReadBouncer only takes the first 360 nucleotides since this represents approximately 0.8 seconds of sequencing.

<b> depletion_files and target_files</b><br>
The depletion-file contains the IBF for references we would want to deplete in a real experiment. If we also have a priori knowledge about potential organisms in our sample we definitely want to sequence, it is also possible to build an IBF of the reference sequences of those organisms and provide ReadBouncer the IBF as a target filter file. This can improve the specificity of the classification.

<b> exp_seq_error_rate </b><br>
For more accurate classification of reads we are calculating the expected number of mutated k-mers for each read prefix based on the expected sequencing ```error rate```. Than a confidence interval for the mutated k-mers is calculated as described by [Blanca et. al., 2021](https://www.biorxiv.org/content/10.1101/2021.01.15.426881v2) and the minimum number of matching k-mers is calculated based on the upper bound of the confidence interval. The significance level of the confidence interval is set to 95% by default.



#### <a name="deplete"></a>Live Depletion of Nanopore Reads

When nanopore reads of a certain organism (or even more) are not of interest, these reads can be unblocked with <b>deplete</b>. Therefore raw current signals are requested from ONT's MinKNOW software, basecalled and matched against the previously build IBF of the reference sequence set. If a read is classified as match with the depletion IBF (and not with the target IBF), an unblock request is sent back to the MinKNOW software. This results in a rejection of the read and another DNA molecule can enter the corresponding pore for sequencing.

```
  target            Live classification and rejection of nanopore reads that match the depletion filter, and not the taget filter (if provided)
  
  [IBF]
  
  target_files        = "xxx.fasta,xxx.fasta"
  deplete_files       = "xxx.fasta,xxx.fasta" 
  threads             = X                      #(unsigned integer with default 3) classification threads
  
  [MinKNOW]
  
  host                = xxx.xxx.xxx            #(ip address or name of the computer hosting MinKNOW) 
  port                = X                      #(port number used fo grpc communication by by MinKNOW instance)
  flowcell            = X                      #(name of the flowcell used)
  
  [Basecaller]
  
  caller              = DeepNano/Guppy         #(default is DeepNano)
  host                = xxx.xxx.xxx            #(ip address or name of the computer hosting Guppy Basecall Server)
  port                = X                      #(port number on which the basecall server is running on the host)
  threads             = X                      #(unsigned integer with default 3) basecall threads
  
```

<b> depletion_files and target_files</b><br>
The depletion-file contains the IBF for references we would want to deplete in a real experiment or you can provide the file as a FASTA file and assign the parameters to build the IBF. ReadBouncer builds automatically an IBF for each provided FASTA file. If we also have a priori knowledge about potential organisms or genomic regions of interest in our sample that should not be rejected, it is also possible to build an IBF of the reference sequences of those organisms or genomic regions and provide ReadBouncer this information as target IBF file. This can improve the specificity of the classification.

<b> exp_seq_error_rate </b><br>
For more accurate classification of reads we are calculating the expected number of mutated k-mers for each read prefix based on the expected sequencing ```error rate```. Than a confidence interval for the mutated k-mers is calculated as described by [Blanca et. al., 2021](https://www.biorxiv.org/content/10.1101/2021.01.15.426881v2) and the minimum number of matching k-mers is calculated based on the upper bound of the confidence interval. The significance level of the confidence interval is set to 95% by default.

<b> [MinKNOW] host and port </b><br>
This is the IP adress and the TCP/IP port on which the MinKNOW software is hosted. ReadBouncer will exchange data with the MinKNOW software via this communication channel. It is recommended to test the communication before starting the sequencing run.

<b> flowcell</b><br>
This is the name of the FlowCell for which we want to do the live depletion. 

<b> [Basecaller] caller </b><br>
Basecaller used during adaptive sampling (default: DeepNano).

<b> [Basecaller] host </b><br>
IP address of guppy basecall server (default: 127.0.0.1).

<b> [Basecaller] port </b><br>
TCP/IP port of guppy basecall server (default: 5555).

<b> [Basecaller] threads </b><br>
Number of threads used for base calling.

### <a name="ucase"></a>Use Cases (TODO from here --> to be checked)

### <a name="classifyreads"></a>Classify already sequenced reads
Sometimes it can be useful to find all reads of an organism in a set of reads that were already sequenced without aligning the sequences. ReadBouncer offers this functionality by using the `classify` subcommand. The following steps describe how to classify all bacterial reads from a Zymo Mock Community that was sequenced on a MinION device.

1. Download the bacterial reference sequences of the Zymo Mock Community from [here](https://owncloud.hpi.de/s/di1lwRsvkXAr4XN) and store it in your working directory.
2. Build an Interleaved Bloom Filter (IBF) file from those reference sequence set by providing the necessary parameters in the config.toml file.

```
usage         = "build"
output_dir    = "full\path\to\ReadBouncer\output_dir\"
log_directory = "full\path\to\ReadBouncer\"

[IBF]

kmer_size     = 13                       #(unsigned integer with default 13)
fragment_size = 100000                   #(unsigned integer with default 100000)
threads       = 3                        #(unsigned integer with default 3) classification threads
target_files  = "path\to\reference\file\ZmcBacterialReferences.fasta"
deplete_files = "" 
```
Using command line: 
```
full\path\to\ReadBouncer\root\directory\bin\ReadBouncer.exe  full\path\to\ReadBouncer\config.toml 
```
3. Use the `classify` subcommand get all reads that origin from one of the 7 bacteria in the Zymo Mock Community mix. (You can use the following [sample read set]() to simply test the feature)

### <a name="unblockall"></a>Test ReadBouncer-to-MinKNOW interaction
Before using ReadBouncer in a real experiment, we recommend running a playback experiment to test unblock speed first.

1. Download an open access bulk FAST5 file from 
[here](http://s3.amazonaws.com/nanopore-human-wgs/bulkfile/PLSP57501_20170308_FNFAF14035_MN16458_sequencing_run_NOTT_Hum_wh1rs2_60428.fast5). 
This file is 21Gb so make sure you have plenty of space.
2. To configure a run for playback, you need to find and edit a sequencing TOML file. These are typically located in `C:\Program Files\OxfordNanopore\MinKNOW\conf\package\sequencing` on Windows or `/opt/ont/minknow/conf/package/sequencing` on Linux. Edit a file such as sequencing_MIN106_DNA.toml and under the entry `[custom_settings]` 
add a field: 
    ```text
    simulation = "/full/path/to/your_bulk.FAST5"
    ```
3. On Windows, open a Windows Power Shell with administrator privileges and go to the MinKNOW binary directory, located in `C:\Program Files\OxfordNanopore\MinKNOW\bin`, and call the config editor with the following two commands:
```
.\config_editor.exe --conf application --filename ..\conf\app_conf --set data_generation.simulated_device=1
.\config_editor.exe --conf application --filename ..\conf\app_conf --set device.simulator_device_count=1
```
On Linux, open a terminal and change to the MinKNOW binary directory, located in `/opt/ont/minknow/bin`, and call the config editor with the following two commands:
```
sudo ./config_editor --conf application --filename ../conf/app_conf --set data_generation.simulated_device=1
sudo ./config_editor --conf application --filename ../conf/app_conf --set device.simulator_device_count=1
```

4. In the MinKNOW GUI, right click on a sequencing position and select `Reload Scripts` (In some cases you need to reboot your operating system). Your MinKNOW instance will now show a simulated device named `MS00000` that will playback the bulkfile rather than live sequencing.
5. Open a Windows Power Shell (or terminal) and go to your working directory where ReadBouncer result files shall be stored. Then provide the necessary parameters in the config.toml file for the live-depletion test. ReadBouncer will test the conncetion to MinKNOW automatically.  

```
[MinKNOW]

host               = "localhost"          #(ip address or name of the computer hosting MinKNOW)
port               = 9501                 #(port number used fo grpc communication by by MinKNOW instance)
flowcell           = "MS00000"            #(name of the flowcell used)

```
The output should state that the connection could be successfully established and that you can continue with live-depletion.
6. When ReadBouncer says that it successfully established a connection, you can start a sequencing run on the the device, which will playback the run from the bulkfile. If you want to use Guppy basecaller, provide the parameters in the configuration file: 

```
[Basecaller]
caller             = "Guppy"       #DeepNano/Guppy (default is DeepNano)
host               = "127.0.0.1"   #(ip address or name of the computer hosting Guppy Basecall Server)
port               = "5555"        #(port number on which the basecall server is running on the host)
threads            = 3             #(unsigned integer with default 3)
```
7. Open the read length histogram after 5 minutes and have a look at the read counts plot.
<p align="center">
  <img src="images/unblock_all.PNG" width="750" title="Unblock All Image">
</p>
8. Now zoom in to the histogram so that only read counts for read lengths up to 5kb are shown. You should see a peak for read counts between 500b and 1 kb like the one in the figure below.
<p align="center">
  <img src="images/unblock_all_5kb.PNG" width="750" title="Unblock All Image (5kb)">
</p>
If that's the case you can go on with testing basecalling and classification

### <a name="host-depletion"></a>Live-Basecalling and read classification

In order to test if read depletion works on your machine, you can start a `depletion` playback run with the bulk FAST5 file from the test above. If you already set up the playback functionality, you only need to download the reference sequence of one or more human chromosomes from e.g. the [Telomere-to-telomere consortium](https://github.com/marbl/CHM13) as FASTA file. In the example below, we aim to deplete all human reads.

1. Before depletion, we need to build an Interleaved Bloom Filter (IBF) for the reference sequence(s) we aim to deplete. 

```
usage         = "build"
output_dir    = "full\path\to\ReadBouncer\output_dir\"
log_directory = "full\path\to\ReadBouncer\"

[IBF]

kmer_size     = 13                      #(unsigned integer with default 13)
fragment_size = 100000                  #(unsigned integer with default 100000)
threads       = 3                       #(unsigned integer with default 3) classification threads
target_files  = "path\to\reference\file\chm13.fasta"
deplete_files = "" 
```
Using command line: 

```
full\path\to\ReadBouncer\root\directory\bin\ReadBouncer.exe  full\path\to\ReadBouncer\config.toml 
```

2. Now you can start depletion of human reads with the following subcommand from you working directory
 ```
full\path\to\ReadBouncer\root\directory\bin\ReadBouncer.exe deplete -d path\to\output\directory\chm13.ibf -f MS00000 -c 3 -b 3 -v
```
3. Start a sequencing run on the simulated device as you did above. Open the read length histogram after 15 minutes and have a look at the read counts plot. When you zoom into the region for reads up to 5kb length, you should see a plot like this:
<p align="center">
  <img src="images/live_unblock.png" width="750" title="Live Unblock Image">
</p>

4. After stopping the run, ReadBouncer will provide you with some statistics about the number of classified (unblocked) and unclassified reads, which will be sequenced until the end. You will also see average overall processing times as well as for basecalling and classification. You should aim for overall processing times for classified reads below one second. The average processing time for basecalling and classification should be below 0.01 seconds. Otherwise you will experience bigger lengths of unblocked reads.
