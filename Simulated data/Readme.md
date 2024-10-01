# Simulated adata generation scripts

This repository contains a collection of scripts for generating simulated dataset for (ITS) sequences in fungal species. These tools cover various aspects of ITS analysis, from sequence trimming to simulation and filtering.

## Scripts Overview

1. [ITSX.sh](#itsx.sh)
2. [cutadapt.sh](#cutadapt.sh)
3. [fasta_header.py](#fasta_header.py)
4. [find_fullITS_fasta.sh](#find_fullits_fasta.sh)
5. [fq2fa.sh](#fq2fa.sh)
6. [simulation.py](#simulation.py)
7. [table_filter_species.py](#table_filter_species.py)

## Detailed Script Descriptions

### ITSX.sh

This Bash script runs ITSx on all FASTA files in the current directory. It creates separate output directories for each input file.

### cutadapt.sh

This script uses Cutadapt to trim ITS and SSU sequences from FASTA files.

**Features:**
- Trims ITS sequences with specific adapters
- Trims SSU sequences with specific adapters
- Outputs trimmed sequences to separate FASTA files

### fasta_header.py

This Python script extracts headers with common words from a FASTA file and writes them to a TSV file.

**Usage:**
```bash
python fasta_header.py <input_fasta> <output_tsv>
```

**Features:**
- Finds common words in species names
- Extracts headers containing these common words
- Outputs results to a TSV file

### find_fullITS_fasta.sh

This Bash script finds and copies full ITS FASTA files from a source directory to a destination directory.

**Features:**
- Searches for full and partial ITS FASTA files
- Excludes ITS1 and ITS2 specific files
- Copies matched files to a specified destination directory

### fq2fa.sh

simple bash script for changing the fastq files into the FASTA format

### simulation.py

This Python script simulates sequencing reads from ITS sequences using SimLoRD.

**Features:**
- Parses input FASTA file
- Generates simulated reads for each sequence
- Creates a report of sequence counts in output FASTQ files

### table_filter_species.py

This Python script filters a CSV table and generates a FASTA file with species names in the headers.

**Features:**
- Filters rows in a CSV file
- Generates a new FASTA file with updated headers including species names