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
[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)

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
