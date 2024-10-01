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


### fq2fa.sh

simple bash script for changing the fastq files into the FASTA format

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

# Chimera Generator Script

This Python script generates chimeric reads by combining sequences from multiple FASTA files. It's designed to simulate chimeric sequences that can occur in sequencing experiments, particularly useful for testing chimera detection algorithms.

## Features

1. Generates chimeric sequences from input FASTA files
2. Creates a detailed log of chimera information
3. Maintains original sequence abundance ratios
4. Includes both intra-file and inter-file chimeras
5. Optionally reverses and complements some chimeras
6. Avoids generating too many short chimeras

## Dependencies

- Python 3.x
- Biopython


## Input

- One or more FASTA files in the same directory as the script

## Output

1. A new directory named `chimeric_reads` containing:
   - FASTA files with original and chimeric sequences
   - A TSV file `chimera_info.tsv` with detailed information about each chimera

## Detailed Workflow

1. **Setup:**
   - Identifies all FASTA files in the current directory
   - Creates an output directory `chimeric_reads`

2. **Chimera Generation:**
   - For each input FASTA file:
     - Determines the number of chimeras to generate (1-3% of total reads)
     - Generates chimeras by:
       - Selecting two parent sequences
       - Choosing a random breakpoint
       - Combining the sequences at the breakpoint
     - Occasionally reverse complements the chimeric sequence (every 25th chimera)

3. **Abundance Ratio Preservation:**
   - Calculates original abundance ratios of sequences
   - Assigns new chimeras an abundance ratio between 0.1 and 2/3 of their main parent's original ratio

4. **Quality Control:**
   - Avoids generating chimeras where either component is less than 50 base pairs
   - Limits the number of short chimeras to 10% of total chimeras
   - Ensures chimera length is not too different from parent sequences

5. **Output Generation:**
   - Inserts chimeric sequences randomly into the original sequence list
   - Writes new FASTA files with original and chimeric sequences
   - Creates a detailed TSV log with information about each chimera

## Chimera Information Log

The `chimera_info.tsv` file contains the following information for each chimera:

1. Chimera ID
2. First parent sequence ID
3. First parent sequence file
4. Second parent sequence ID
5. Second parent sequence file
6. Breakpoint position
7. Whether the chimera was reverse complemented
8. Assigned abundance ratio
9. Chimera length

## Note

This script is designed for generating test data and should not be used to create actual biological sequences for analysis. The chimeras generated are artificial and may not represent real biological chimeras accurately.