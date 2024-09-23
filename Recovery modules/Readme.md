## RECOVERY MODULES FOR THE CHIMERIC READS

### BLAST
#### Dependencies: Python packages: Biopython, pandas

To apply the BLAST module to the chimeric queries, the *BLASTn.sh* module must run first on them; later, the *BLAST_recovery.py* module can rescue the nonchimeric reads. The module can be run by the command below specifying the directories related to each output:

```bash
python BLAST_recovery.py --chimeras_dir ./chimeras --blast_output_dir ./blast_output --nonchimeric_dir .
```

### ReChime
Please read the ReChime_RUNME.txt file for the instructions in the **ReChime_v1.zip** file before running the ReChime module. 


# 🧬 BlasCh recovery module for recovering False positive chimeras 

[![Python Version](https://img.shields.io/badge/python-3.6%2B-blue)](https://www.python.org/downloads/)

Efficient chimera detection for long-read sequencing data using multiprocessing.

## 🚀 Features

- **Multiprocessing**: Blazing fast processing of multiple BLAST XML files
- **Resource Monitoring**: Tracks CPU and memory usage for optimization
- **Detailed Reporting**: Generates comprehensive overall and per-file statistics
- **Adaptive Performance**: Automatically utilizes all available CPU cores
- **Robust Error Handling**: Implements logging for seamless debugging


## 🏃‍♂️ Usage

1. Place your input FASTA files in the `./input` directory
2. Ensure BLAST XML result files are in the current directory
3. Run the script:
   ```
   python BlasCh.py
   ```
4. Find results in the `./rescued_reads` directory

## ⚙ Configuration

Modify these variables at the script's beginning to customize behavior:

```python
HIGH_IDENTITY_THRESHOLD = 99.0
HIGH_COVERAGE_THRESHOLD = 99.0
SIGNIFICANT_COVERAGE_THRESHOLD = 80.0
SIGNIFICANT_IDENTITY_THRESHOLD = 80.0
```
# BLAST Alignment Thresholds Explanation

The script uses two sets of thresholds to categorize BLAST alignments:

1. High Identity Thresholds:
   - Identity >= 99% and Coverage >= 99%

2. Significant Alignments Thresholds:
   - Coverage >= 80% and Identity >= 80%

## Why Two Sets of Thresholds?

These two sets of thresholds serve different purposes in the chimera detection process:

1. High Identity Thresholds (99%/99%):
   - Purpose: To identify nearly perfect matches.
   - Interpretation: Alignments meeting these criteria suggest that the query sequence is almost identical to a known sequence in the database.
   - Use in classification: Used to detect potential false positive chimeras or uncertain chimeras.

2. Significant Alignments Thresholds (80%/80%):
   - Purpose: To identify meaningful, but not necessarily perfect, matches.
   - Interpretation: Alignments meeting these criteria suggest that a substantial portion of the query sequence is similar to a known sequence, but allows for some differences.
   - Use in classification: Used to detect potential absolute chimeras or uncertain chimeras.

## How the Thresholds Work

For each alignment in the BLAST results:

1. Calculate query coverage: (alignment length/query sequence length) * 100
2. Calculate identity percentage: (number of identical matches/alignment length) * 100
3. Compare these values to the thresholds:
   - If both values meet or exceed the high identity thresholds (99%/99%), categorize them as a "high identity" alignment.
   - If both values meet or exceed the significant alignment thresholds (80%/80%) but don't meet the high identity thresholds, categorize them as a "significant" alignment.
   - If neither set of thresholds is met, the alignment is not considered further in the chimera classification process.

## Impact on Classification

- Sequences with one or more high identity alignments are classified as either false positive chimeras or uncertain chimeras, depending on whether the alignments are to different species.
- Sequences with multiple significant alignments (but no high identity alignments) are classified as absolute chimeras.
- Sequences with only one significant alignment are classified as uncertain chimeras.
- Sequences with no high identity or significant alignments are classified as non-chimeric.

This two-tiered approach allows the script to distinguish between nearly identical matches and merely significant matches, providing a more nuanced classification of potential chimeric sequences.

## Initial Alignment Check (in Blasch_modified version):

The script first parses the BLAST results, where each sequence (query) is aligned against a database of reference sequences (hits).
For each query, the alignments are evaluated based on two key metrics:
Query Coverage: The percentage of the query sequence that aligns with the reference (hit) sequence.
Identity Percentage: The percentage of identical matches in the alignment.
High-Identity Alignment:

Sequences with high identity (≥ 99%) and high query coverage (≥ 99%) are considered for stricter classification as chimeras.
If a sequence has multiple high-identity alignments or if the alignment is to a different taxonomic group, it's flagged as an absolute chimera or a false positive chimera.
Borderline Sequence Criteria: If no high-identity alignment is found (meaning the sequence does not meet the 99% identity and coverage threshold), the script evaluates the remaining alignments using less stringent criteria:

Max Identity and Coverage:
The script calculates the maximum identity percentage and maximum query coverage across all alignments for that sequence.
If any alignment has a minimum identity percentage of 80% and query coverage of 80%, the sequence is flagged as borderline.

## 📂 Directory Structure

```
.
├── input/                  # Input FASTA files
├── rescued_reads/          # Output directory
├── BlasCh.py
└── README.md
```

## 📊 Output

- Classified sequences in separate FASTA files
- Detailed report: `chimera_detection_report.txt`
- Log file with process info and resource usage
