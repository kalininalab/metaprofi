# UniProtKB organism-level bacterial dataset index construction

This directory contains all the information on how we constructed **organism-level** MetaProFi index for the UniProtKB bacterial dataset.

### Data download

``` bash
# In July 2021, data was downloaded from the current release section like so from UniProt's FTP site (Now running this will download the latest dataset, with additional sequences)

# SwissProt dataset
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# TrEMBL dataset
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
```

* List of bacterial sequence accession ids were downloaded using the UniProt's search interface with the keyword `taxonomy:bacteria`

```
https://www.uniprot.org/uniprot/?query=taxonomy%3Abacteria&sort=score
```

### Data preparation

- Both datasets were parsed and every bacterial sequence was extracted based on their organism name in their FASTA header `OS=OrganismName` to a `OrganismName` specific FASTA file using this [custom script](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/uniprotkb_organism_level_index/uniprotkb_parser.py) with a minimum sequence length criteria (sequence length >= _k_-mer size)

- Size of the compressed files
    - SwissProt: 87 MiB
    - TrEMBL:  50 GiB
    - uniprot-taxonomy_bacteria.list.gz: 611 MiB

_Note: Don't forget to change the file and ouput path in the `uniprotkb_parser.py`_

``` bash
# Activate metaprofi environment
conda activate metaprofi

# Creates input files for organism-level and sequence-level indexing
python uniprotkb_parser.py
```

### MetaProFi commands

- [config](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/uniprotkb_organism_level_index/config.yml) and [input](https://github.com/kalininalab/metaprofi/blob/master/benchmarks/uniprotkb_organism_level_index/input_data.txt) files were created

``` bash
# Activate conda environment
conda activate metaprofi

# Run MetaProFi
metaprofi build input_data.txt config.yml
```

### Results

#### Bloom filter matrix construction

| Time (min) | RAM (GiB) | CPU cores | Disk (GiB) |
| --- | --- | --- | --- |
| 38.05 | <60 | 64 | 139 |

#### Index construction

| Time (min) | RAM (GiB) | CPU cores | Disk (GiB) |
| --- | --- | --- | --- |
| 971 | <60 | 64 | 135 |