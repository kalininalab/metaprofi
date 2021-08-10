"""UniProtKB parser"""

import re
from pyfastx import Fasta
import gzip
import os
from pathlib import Path
import subprocess

# File and output paths (Change me if required)
bacterial_accessions = "/data/uniprotkb_datasets/uniprot-taxonomy_bacteria.list.gz"
swissprot = "/data/uniprotkb_datasets/uniprot_sprot.fasta.gz"
trembl = "/data/uniprotkb_datasets/uniprot_trembl.fasta.gz"
org_level = "/data/uniprotkb_datasets/organism_level"
seq_level = "/data/uniprotkb_datasets/seq_level"
metaprofi_input_file = "/data/uniprotkb_datasets/input_data.txt"

# Create directories if it does not exist
Path(org_level).mkdir(parents=True, exist_ok=True)
Path(seq_level).mkdir(parents=True, exist_ok=True)

# Minimum sequence length criteria (change me if required)
# Set to the length of the k-mer
min_seq_len = 11

# Regex
# 1: To get the UniProt accession id (UniqueIdentifier) from header
uniprot_id = re.compile(r".*\|(.*)\|")

# 2: To get the OS name from the FASTA header
organism_name = re.compile(r".* OS=(.*?)OX=.*")

# 3: To remove '(', ')' and fill gaps with '_'
sub_pattern = re.compile(r"[^A-Za-z0-9 ]+")

# 1: Read the bacterial accession ids into the memory (might consume approx. 14 GiB of memory)
accession_ids = set()
with gzip.open(bacterial_accessions, "rt") as accessions:
    for line in accessions:
        accession_ids.add(line.strip())

# Sequence-level index input file
# 2: Extract bacterial sequences from both swissprot and trembl and combine
# them into one file
for header, seq in Fasta(trembl, build_index=False, full_name=True):
    if len(seq) >= min_seq_len:
        if uniprot_id.match(header).group(1) in accession_ids:
            with open(f"{seq_level}/uniprot_bacteria.fasta", "w") as out_file:
                _ = out_file.write(f">{header}\n{seq}\n")

# Compress the uniprot_bacteria.fasta file
_ = subprocess.run(["pigz", "-p4", f"{seq_level}/uniprot_bacteria.fasta"])

# Organism-level index input files
# 3: Extract sequences and group them by their organism name
for header, seq in Fasta(
    f"{seq_level}/uniprot_bacteria.fasta.gz",
    build_index=False,
    full_name=True,
):
    oname = organism_name.match(header).group(1).strip()
    oname = sub_pattern.sub("", oname).replace(" ", "_")
    with open(f"{org_level}/{oname}.fasta", "a") as out_file:
        _ = out_file.write(f">{header}\n{seq}\n")

# MetaProFi input file construction
files_list = os.listdir(org_level)
with open(metaprofi_input_file, "w") as out_file:
    for f in files_list:
        # File name as the sample identifier
        sample_id = f.split(".")[0]
        _ = out_file.write(f"{sample_id}: {org_level}/{f}\n")
