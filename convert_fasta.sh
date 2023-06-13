#!/bin/bash

# Create FASTA folder
mkdir -p FASTA

# Log file for error messages
error_log="error.log"

# Find all fastq files and convert them to FASTA using seqkit
for fastq_file in *.fastq; do
    fasta_file="FASTA/${fastq_file%.fastq}.fasta"
    
    # Run seqkit fq2fa and redirect stderr to the error log file
    seqkit fq2fa "$fastq_file" -o "$fasta_file" 2>> "$error_log"
    
    # Check if the FASTA file was generated
    if [ -f "$fasta_file" ]; then
        echo "Conversion successful: $fastq_file -> $fasta_file"
    else
        echo "Conversion failed: $fastq_file" >> "$error_log"
    fi
done
