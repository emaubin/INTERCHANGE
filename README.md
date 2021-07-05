# INTERCHANGE

NTERCHANGE for horIzoNtal TransfER CHAracterization in Non-assembled GEnome is a pipeline for the detection of horizontal tranfers in non-assembled genomes by characterization of reads from conserved regions between 2 species.


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

- [ncbi blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  - A suite of command-line tools to run BLAST

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

- [minimap2](https://github.com/lh3/minimap2)
  - Minimap2 is a sequence alignment programm that aligns DNA or mRNA sequences against a large reference database.

  ```bash
  # Installation with Ubuntu/Debian
  sudo apt-get install minimap2
  ```

  ```bash
  # download Minimap2
  curl -L https://github.com/lh3/minimap2/releases/download/v2.20/minimap2-2.20_x64-linux.tar.bz2
  # decompress
  tar -jxvf minimap2-2.20_x64-linux/minimap2
  ```
- [Samtools Version: 1.10](http://www.htslib.org/)
  - Samtools is a suite of programs for interacting with high-throughput sequencing data

  ```bash
  # Installation with Ubuntu/Debian
  sudo apt-get install samtools
  ```

  Or download Samtools [here](http://www.htslib.org/download/)
  ```bash
  cd samtools-1.x
  ./configure --prefix=/where/to/install
  make
  make install
  ```

The pipeline has not been tested with other versions of the above programs, but newer versions probably work by checking that the options used still exist


Hardware requirements: this pipeline is developed for Linux/Unix operating system
With the test dataset, we used:
    - x86-64 CPUs
    - 32 Go of system memory
    - of free hard drive space


## <a name="uguide"></a>Users' Guide

### <a name="install"></a>Installation

In bash compatible terminal:

```bash
# download INTERCHANGE version 1.0
wget https://github.com/emaubin/INTERCHANGE/archive/V.1.0.zip
# decompress
unzip  V.1.0.zip
```
### <a name="usage"></a>Usage

You must start by filling in the dependencies_paths.txt with the paths to each tool and databases used as indicated in the file. Then, run [step]_param.py scripts whose name start with numbers in the corresponding order. Adapting this pipeline to other datasets, hardware configuration, and automating all procedures require modifications to the code.
To know the arguments needed for each step you can use the help for each script as follows:

```bash
python3 /INTERCHANGE-V.1.0/scripts/1.Genome_format/format_param.py -h
"""
usage: python3 format_param.py -i -p

Script to prepare and format input data for INTERCHANGE

Positional arguments:
  -i TABLE             Input file containing Table of species.
  -p PATHS             File of tools paths.

Settings:
  -t THREAD            Number of CPU for gzip/gunzip. Default [2]

Output options:
  -o OUTPUT_DIRECTORY  Output directory for INTERCHANGE results. Default: /INTERCHANGE_results in current directory

Other:
  -h, --help           Show this help message and exit.
  -v, --version        Show program's version number and exit.
"""
```
