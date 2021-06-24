# INTERCHANGE

NTERCHANGE for horIzoNtal TransfER CHAracterization in Non-assembled GEnome is pipeline for the detection of horizontal tranfer in non-assembled genomes by characterization of reads from conserved regions between 2 species.


## Table of Contents

- [Requirements](#req)
- [User's Guide](#uguide)
  - [Installation](#install)
  - [Usage](#usage)
  - [Example](#example)
- [Limitations](#limit)


## <a name="req"></a>Requirements

- Python3.8 with the following packages
  - [biopython version 1.78](https://biopython.org/)
  - [pandas version 1.2.0](https://pandas.pydata.org/)

  ```bash
  # Command if you work with multiple version of Python installed in parallel
  python3 -m pip install biopython # default Python 3
  python3.8 -m pip install biopython # specifically Python 3.8
  python3 -m pip install pandas # default Python 3
  python3.8 -m pip install pandas # specifically Python 3.8
  ```

- [Pigz version 2.4](https://zlib.net/pigz/)
  - A tool for gzip that exploits multiple processors and multiple cores to the hilt when compressing data

  ```bash
  # Installation with Ubuntu/Debian
  sudo apt-get install pigz
  ```

- [GenomeTools versio 1.6.1](http://genometools.org/)
    - The GenomeTools genome analysis system is a free collection of bioinformatics tools

    ```bash
    # Installation with Ubuntu/Debian
    sudo apt-get install genometools
    ```

    ```bash
    # download GenomeTools versio 1.6.1
    wget http://genometools.org/pub/genometools-1.6.1.tar.gz
    # decompress
    tar zxvf genometools-1.6.1.tar.gz
    cd genometools-1.6.1
    # Install in current directory
    sudo make install

    ```
- [PRINSEQ LITE version 0.20.4](https://github.com/uwb-linux/prinseq)
  - PRINSEQ is a tool to preprocess genomic or metagenomic sequence data in FASTA or FASTQ format
  - The lite version is a standalone perl script (prinseq-lite.pl) that does not
require any non-core perl modules for processing

  ```bash
  # download PRINSEQ LITE version 0.20.4
  wget https://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz
  # decompress
  tar zxvf prinseq-lite-0.20.4.tar.gz
  # make it executable
  cd prinseq-lite-0.20.4
  sudo chmod +x prinseq-lite.pl
  ```

- [SPAdes version 3.15.2](https://github.com/ablab/spades)
  - SPAdes is an assembly toolkit containing various assembly pipelines

  ```bash
  # download SPAdes version 3.15.2
  wget http://cab.spbu.ru/files/release3.15.2/SPAdes-3.15.2-Linux.tar.gz
  # decompress
  tar zxvf SPAdes-3.15.2-Linux.tar.gz

  ```

- [DIAMOND version 0.9.29](https://github.com/bbuchfink/diamond)
    - DIAMOND is a sequence aligner for protein and translated DNA searches with high performance analysis of big sequence data

    ```bash
    # Installation with Ubuntu/Debian
    sudo apt-get install diamond-aligner
    ```

    ```bash
    # download DIAMOND version 0.9.29
    wget http://github.com/bbuchfink/diamond/releases/download/v2.0.9/diamond-linux64.tar.gz
    # decompress
    tar zxvf diamond-linux64.tar.gz
    ```

- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  -

  ```bash
  # Installation with Ubuntu/Debian
  sudo apt-get install ncbi-blast+
  ```

  ```bash
  # download BLAST version 2.9.0
  wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz
  # decompress
  tar -zxvf ncbi-blast-2.8.1+-x64-linux.tar.gz
  # add BLAST location to system PATH
  export PATH=$HOME/tools/BLAST/ncbi-blast-2.8.1+/bin:$PATH
  ```

The pipeline has not been tested with other versions of the above programs, but newer versions probably work by checking that the options used are still working


Hardware requirements: this pipeline is developed for Linux operating system
With the test dataset, we used:
    - x86-64 CPUs
    - 32 Go of system memory
    - of free hard drive space


## <a name="uguide"></a>Users' Guide

### <a name="install"></a>Installation

In bash compatible terminal:

```bash
# download INTERCHANGE version 1.0
wget https://github.com/emaubin/INTERCHANGE/archive/main.zip
# decompress
tar -zxvf INTERCHANGE_linux_v1.0.tar.gz
```
### <a name="usage"></a>Usage
