#!/bin/bash

# Create FASTA folder
mkdir -p FASTA

# Find all fastq files and convert them to FASTA using seqkit
for fastq_file in *.fastq; do
    fasta_file="FASTA/${fastq_file%.fastq}.fasta"
    seqkit fq2fa "$fastq_file" -o "$fasta_file"
done
