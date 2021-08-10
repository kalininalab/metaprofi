# MetaProFi benchmarks

This directory contains all the commands, analyses and results obtained by MetaProFi.

## Computing setup

- Dell server
  - CPU: AMD EPYC 7702 2.0 GHz
  - RAM: 1.5 TB
  - Storage: Intel SSD DC P4610 3.2 TB (2.9 TiB)
  - OS: CentOS 7 (kernel 3.10.0-1160.24.1.el7.x86_64)
- I/O
  - All input files were downloaded onto a non-RAID NVME NFS file system
  - All outputs were written to the Intel SSD DC P4610 3.2 TB disk

## Installation

MetaProFi was installed as a conda environment

#### Create conda environment and install MetaProFi

  ``` bash
  # Step 1: Install miniconda
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh

  # Step 2: Download the repo
  git clone https://github.com/kalininalab/metaprofi.git

  # Step 3: Create a new environment and install MetaProFi
  conda create --name metaprofi python==3.7.7 pigz
  conda activate metaprofi
  pip install /path/to/metaprofi/git/repo/directory/
  ```

## Results

Each sub-folder in this directory further contains specific details on how the data was downloaded and what commands were run to obtain results discussed in the paper

- [UniProtKB-organism level index](https://github.com/kalininalab/metaprofi/tree/master/benchmarks/uniprotkb_organism_level_index)
- [UniProtKB-sequence level index](https://github.com/kalininalab/metaprofi/tree/master/benchmarks/uniprotkb_sequence_level_index)
- [Tara Oceans index](https://github.com/kalininalab/metaprofi/tree/master/benchmarks/tara_oceans_index)
- [RNA-seq index](https://github.com/kalininalab/metaprofi/tree/master/benchmarks/rna_seq_index)