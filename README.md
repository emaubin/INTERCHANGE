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
  ```
  # Command if you work with multiple version of Python installed in parallel
  python3 -m pip install biopython # default Python 3
  python3.8 -m pip install biopython # specifically Python 3.8
  python3 -m pip install pandas # default Python 3
  python3.8 -m pip install pandas # specifically Python 3.8
  ```

- [Pigz version 2.4](https://zlib.net/pigz/)
  - A tool for gzip that exploits multiple processors and multiple cores to the hilt when compressing data
  ```bash
    sudo apt-get install pigz
  ```

- [GenomeTools versio 1.6.1](http://genometools.org/)
    - The GenomeTools genome analysis system is a free collection of bioinformatics tools

    ```bash
    # Installation using Ubuntu package management with following command, requires **root permissions**
    sudo apt-get install genometools
    ```

    ```
    # Installation with the source distribution
    wget
    tar zxvf genometools-1.6.1.tar.gz
    cd genometools-1.6.1
    sudo make install
    ```
- [PRINSEQ LITE version 0.20.4](https://github.com/uwb-linux/prinseq)
  - PRINSEQ is a tool to preprocess genomic or metagenomic sequence data in FASTA or FASTQ format
  - The lite version is a standalone perl script (prinseq-lite.pl) that does not
require any non-core perl modules for processing
  ```
  wget https://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz
  tar zxvf prinseq-lite-0.20.4.tar.gz
  cd prinseq-lite-0.20.4
  sudo chmod +x prinseq-lite.pl
  ```
- [SPAdes version 3.14.1](https://github.com/ablab/spades)
  - SPAdes is an assembly toolkit containing various assembly pipelines
  ```
  git clone https://github.com/ablab/spades
  ```
- [DIAMOND version 0.9.29](https://github.com/bbuchfink/diamond)
    - DIAMOND is a sequence aligner for protein and translated DNA searches with high performance analysis of big sequence data
    ```bash
    sudo apt-get install diamond-aligner
    ```






## <a name="uguide"></a>Users' Guide

### <a name="install"></a>Installation
