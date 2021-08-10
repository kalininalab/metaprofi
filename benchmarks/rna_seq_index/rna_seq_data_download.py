"""Downloads RNA-seq dataset from SRA"""

import subprocess

download_dir = "/data/rna_seq/"
accessions_list = []

# Get the accessions list
with open("/data/rna_seq/sra_accessions.txt", "r") as inf:
    for line in inf:
        accessions_list.append(line.strip)

# Download from SRA using 8 threads
for acc_id in accessions_list:
    _ = subprocess.run(
        [
            "parallel-fastq-dump",
            "-s",
            f"{acc_id}",
            "-t",
            "8",
            "-O",
            f"{download_dir}",
            "--gzip",
            "--fasta",
            "0",
        ]
    )
