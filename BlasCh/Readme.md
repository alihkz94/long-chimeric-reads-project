# 🧬 BlasCh recovery module for recovering False positive chimeras 

[![Python Version](https://img.shields.io/badge/python-3.6%2B-blue)](https://www.python.org/downloads/)

Efficient chimera detection and recovery for long-read sequencing data using multiprocessing.

## 🚀 Features

- **Multiprocessing**: Blazing fast processing of multiple BLAST XML files
- **Resource Monitoring**: Tracks CPU and memory usage for optimization
- **Detailed Reporting**: Generates comprehensive overall and per-file statistics
- **Adaptive Performance**: Automatically utilizes all available CPU cores
- **Robust Error Handling**: Implements logging for seamless debugging
- **Sequence Classification**: Categorizes sequences as non-chimeric, chimeric, or borderline

## 🏃‍♂️ Usage

1. Place your input FASTA files in the `./input` directory
2. Ensure BLAST XML result files are in the current directory
3. Run the script:
   ```
   python BlasCh.py
   ```
4. Find results in the `./rescued_reads` directory

## ⚙ Configuration

The script uses the following thresholds for classification:

```python
HIGH_IDENTITY_THRESHOLD = 99.0
HIGH_COVERAGE_THRESHOLD = 99.0
```

These thresholds are used to determine high-quality matches against the database.

## 📂 Directory Structure

```
.
├── input/                  # Input FASTA files
├── rescued_reads/          # Output directory for classified sequences
├── temp/                   # Temporary directory for intermediate results
├── temp_2/                 # Temporary directory for sequence details CSV files
├── BlasCh.py
└── README.md
```

## 📊 Output

- Classified sequences in separate FASTA files:
  - `*_non_chimeric.fasta`: Non-chimeric sequences
  - `*_borderline.fasta`: Borderline sequences
  - `*_chimeric.fasta`: Chimeric sequences (not directly output, but classified)
- Detailed report: `chimera_detection_report.txt` in the `rescued_reads` directory
- Sequence details CSV files in the `temp_2` directory

## 🧠 Classification Logic

The script classifies sequences based on the following criteria:

1. **Non-chimeric**: 
   - No significant non-self hits, or
   - High-quality match against the database (≥99% identity and ≥99% coverage)
2. **Chimeric**: Multiple non-self alignments
3. **Borderline**: Single alignment, requires further analysis

## 🛠 Dependencies

- BioPython
- psutil

## 📝 Logging

The script uses Python's logging module to provide detailed information about the process, including CPU and memory usage.

## ⚠️ Error Handling

The script implements robust error handling and will log any errors that occur during processing.
