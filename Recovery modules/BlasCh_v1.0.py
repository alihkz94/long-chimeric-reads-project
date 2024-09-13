"""
Chimera recovery module for Long-Read Sequencing

UPDATES:
- Implemented multiprocessing for faster processing of multiple BLAST XML files.
- Added CPU and memory usage tracking using psutil.
- Improved error handling and logging with detailed process information.
- Optimized file I/O operations and memory management.
- Enhanced report generation with overall and per-file statistics.

Description:
    This script processes BLAST XML results to identify and classify chimeric sequences in long-read sequencing data.
    It categorizes sequences as false positive chimeras, absolute chimeras, uncertain chimeras, or non-chimeric sequences
    based on specified alignment criteria. The script utilizes multiprocessing for efficient handling of large datasets
    and includes system resource usage monitoring.

Usage:
    python BlasCh.py

Requirements:
    - Python 3.6+
    - BioPython
    - psutil

Dependencies:
    - os
    - shutil
    - multiprocessing
    - time
    - psutil
    - Bio.Blast.NCBIXML (from BioPython)
    - Bio.SeqIO (from BioPython)
    - collections (defaultdict)
    - logging

Input:
    - FASTA files in the specified input directory
    - BLAST XML result files in the current directory

Output:
    - Classified sequences in separate FASTA files within the output directory
    - Detailed report file summarizing overall and per-file results
    - Log file with process information and system resource usage

Configuration:
    Modify the following variables at the beginning of the script to customize paths and thresholds:
    - input_dir: Directory containing input FASTA files
    - directory: Directory containing BLAST XML files (current directory by default)
    - temp_dir: Directory for storing temporary files
    - output_dir: Directory for storing output files
    - NUM_PROCESSES: Number of processes to use for multiprocessing (automatically set to CPU count)
    - HIGH_IDENTITY_THRESHOLD: Threshold for high identity percentage (99.0)
    - HIGH_COVERAGE_THRESHOLD: Threshold for high query coverage (99.0)
    - SIGNIFICANT_COVERAGE_THRESHOLD: Threshold for significant query coverage (80.0)
    - SIGNIFICANT_IDENTITY_THRESHOLD: Threshold for significant identity percentage (80.0)

Author: ALI HAKIMZADEH
Version: 1.0
Date: 2024-09-11
"""
# Load Libraries 
import os
import shutil
import multiprocessing
import time
import psutil  # For CPU and memory usage tracking
from Bio.Blast import NCBIXML
from Bio import SeqIO
from collections import defaultdict
import logging

# Configure logging to output process information
logging.basicConfig(level=logging.DEBUG, format='%(processName)s: %(message)s')

# Constants for classification criteria
HIGH_IDENTITY_THRESHOLD = 99.0
HIGH_COVERAGE_THRESHOLD = 99.0
SIGNIFICANT_COVERAGE_THRESHOLD = 80.0
SIGNIFICANT_IDENTITY_THRESHOLD = 80.0

# Dynamically determine the number of CPUs to use for multiprocessing
NUM_PROCESSES = multiprocessing.cpu_count()

def log_system_usage():
    """Logs the current CPU and RAM usage."""
    process = psutil.Process(os.getpid())
    cpu_percent = psutil.cpu_percent(interval=1)
    memory_info = process.memory_info()
    logging.info(f"CPU usage: {cpu_percent}%")
    logging.info(f"Memory usage: {memory_info.rss / (1024 ** 2):.2f} MB")

def load_fasta_sequences(fasta_file):
    """Load sequences from a FASTA file into a dictionary."""
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
    
    fasta_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_sequences[record.id] = str(record.seq)
    return fasta_sequences

def extract_query_id(blast_query_def):
    """Extract the query ID from the BLAST query definition."""
    return blast_query_def

def parse_blast_results(args):
    """Parse BLAST XML results and classify sequences into different categories."""
    xml_file, temp_dir, fasta_file = args
    logging.debug(f"Starting processing for {xml_file}")
    
    fasta_sequences = load_fasta_sequences(fasta_file)

    false_positive_chimeras = set()
    absolute_chimeras = set()
    uncertain_chimeras = set()
    non_chimeric_sequences = set()

    with open(xml_file) as result_handle:
        blast_records = list(NCBIXML.parse(result_handle))
        
        for blast_record in blast_records:
            query_id = extract_query_id(blast_record.query)
            if query_id not in fasta_sequences:
                continue

            significant_alignments = []
            high_identity_alignments = []
            
            for alignment in blast_record.alignments:
                hit_id = extract_query_id(alignment.hit_def)
                if hit_id != query_id:  # Exclude self-hits
                    for hsp in alignment.hsps:
                        query_coverage = (hsp.align_length / blast_record.query_length) * 100
                        identity_percentage = (hsp.identities / hsp.align_length) * 100
                        
                        if identity_percentage >= HIGH_IDENTITY_THRESHOLD and query_coverage >= HIGH_COVERAGE_THRESHOLD:
                            high_identity_alignments.append((alignment, hsp))
                        elif query_coverage >= SIGNIFICANT_COVERAGE_THRESHOLD and identity_percentage >= SIGNIFICANT_IDENTITY_THRESHOLD:
                            significant_alignments.append((alignment, hsp))

            # Classification logic
            if high_identity_alignments:
                if all(extract_query_id(align[0].hit_def).split('_')[0] != query_id.split('_')[0] for align in high_identity_alignments):
                    false_positive_chimeras.add(query_id)
                else:
                    uncertain_chimeras.add(query_id)
            elif len(significant_alignments) > 1:
                absolute_chimeras.add(query_id)
            elif significant_alignments:
                uncertain_chimeras.add(query_id)
            else:
                non_chimeric_sequences.add(query_id)

        # Write intermediate results to the temporary folder
        base_filename = os.path.basename(fasta_file).replace(".chimeras.fasta", "")
        write_sequences_to_file(false_positive_chimeras, fasta_sequences, os.path.join(temp_dir, f"{base_filename}_false_positive_chimeras.fasta"))
        write_sequences_to_file(non_chimeric_sequences, fasta_sequences, os.path.join(temp_dir, f"{base_filename}_non_chimeric.fasta"))
    
    logging.debug(f"Completed processing for {xml_file}")
    return false_positive_chimeras, absolute_chimeras, uncertain_chimeras, non_chimeric_sequences

def write_sequences_to_file(seq_ids, fasta_sequences, output_file):
    """Write sequences to a FASTA file."""
    with open(output_file, 'w') as out_file:
        for seq_id in sorted(seq_ids):
            out_file.write(f">{seq_id}\n{fasta_sequences[seq_id]}\n")

def combine_sequences(temp_dir, base_filename, output_dir):
    """Combine non-chimeric and false-positive chimera sequences into a final output."""
    false_positive_file = os.path.join(temp_dir, f"{base_filename}_false_positive_chimeras.fasta")
    non_chimeric_file = os.path.join(temp_dir, f"{base_filename}_non_chimeric.fasta")
    recovered_file = os.path.join(output_dir, f"{base_filename}_recovered.fasta")

    with open(recovered_file, 'w') as outfile:
        for fname in [non_chimeric_file, false_positive_file]:
            if os.path.exists(fname):
                with open(fname) as infile:
                    shutil.copyfileobj(infile, outfile)

def generate_report(all_false_positive_chimeras, all_absolute_chimeras, all_uncertain_chimeras, all_non_chimeric_sequences, file_results, output_dir):
    """
    Generate a detailed report summarizing the results, sorted by file name.
    
    Args:
        all_false_positive_chimeras (set): Set of all false positive chimeras.
        all_absolute_chimeras (set): Set of all absolute chimeras.
        all_uncertain_chimeras (set): Set of all uncertain chimeras.
        all_non_chimeric_sequences (set): Set of all non-chimeric sequences.
        file_results (dict): Dictionary containing results for each processed file.
        output_dir (str): Directory for output files.
    """
    report = {
        "False Positive Chimeras": len(all_false_positive_chimeras),
        "Absolute Chimeras": len(all_absolute_chimeras),
        "Uncertain Chimeras": len(all_uncertain_chimeras),
        "Non-Chimeric Sequences": len(all_non_chimeric_sequences)
    }

    total_sequences = sum(report.values())

    with open(os.path.join(output_dir, "chimera_detection_report.txt"), 'w') as report_file:
        report_file.write("Chimera Detection Report\n")
        report_file.write("========================\n\n")
        report_file.write("Overall Summary:\n")
        report_file.write("----------------\n")
        for category, count in report.items():
            percentage = (count / total_sequences) * 100 if total_sequences > 0 else 0
            report_file.write(f"{category}: {count} ({percentage:.2f}%)\n")
        report_file.write(f"\nTotal Sequences Processed: {total_sequences}\n\n")
        
        report_file.write("Detailed Results by File:\n")
        report_file.write("-------------------------\n")
        
        # Sort the file_results dictionary by keys (file names)
        sorted_file_results = dict(sorted(file_results.items()))
        
        for file, results in sorted_file_results.items():
            report_file.write(f"\nFile: {file}\n")
            file_total = sum(results.values())
            for category, count in results.items():
                percentage = (count / file_total) * 100 if file_total > 0 else 0
                report_file.write(f"  {category}: {count} ({percentage:.2f}%)\n")
            report_file.write(f"  Total: {file_total}\n")


def process_all_xml_files(directory, temp_dir, output_dir, input_dir):
    """Process all BLAST XML files and generate final outputs."""
    # Clean up existing directories if they exist
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    # Recreate directories after cleanup
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    args_list = []
    file_results = defaultdict(lambda: defaultdict(int))
    
    all_false_positive_chimeras = set()
    all_absolute_chimeras = set()
    all_uncertain_chimeras = set()
    all_non_chimeric_sequences = set()

    start_time = time.time()  # Track the start time
    log_system_usage()  # Log initial system usage

    for filename in os.listdir(directory):
        if filename.endswith("_blast_results.xml"):
            xml_file = os.path.join(directory, filename)
            fasta_file = os.path.join(input_dir, filename.replace("_blast_results.xml", ".chimeras.fasta"))
            args_list.append((xml_file, temp_dir, fasta_file))

    logging.debug(f"Using {NUM_PROCESSES} processes for multiprocessing")
    
    with multiprocessing.Pool(processes=NUM_PROCESSES) as pool:
        results = pool.map(parse_blast_results, args_list)

    for i, result in enumerate(results):
        false_positive_chimeras, absolute_chimeras, uncertain_chimeras, non_chimeric_sequences = result
        all_false_positive_chimeras.update(false_positive_chimeras)
        all_absolute_chimeras.update(absolute_chimeras)
        all_uncertain_chimeras.update(uncertain_chimeras)
        all_non_chimeric_sequences.update(non_chimeric_sequences)
        
        file_results[args_list[i][0]] = {
            "False Positive Chimeras": len(false_positive_chimeras),
            "Absolute Chimeras": len(absolute_chimeras),
            "Uncertain Chimeras": len(uncertain_chimeras),
            "Non-Chimeric Sequences": len(non_chimeric_sequences)
        }

    # Combine sequences and generate final report
    for filename in os.listdir(input_dir):
        if filename.endswith(".chimeras.fasta"):
            base_filename = os.path.basename(filename).replace(".chimeras.fasta", "")
            combine_sequences(temp_dir, base_filename, output_dir)
    
    generate_report(all_false_positive_chimeras, all_absolute_chimeras, all_uncertain_chimeras, all_non_chimeric_sequences, file_results, output_dir)

    log_system_usage()  # Log final system usage

    # Clean up temporary folder after processing
#    shutil.rmtree(temp_dir)

    end_time = time.time()
    elapsed_time = end_time - start_time
    logging.info(f"Total time taken: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    input_dir = "./xml_uchime/input"
    directory = "./xml_uchime"
    temp_dir = "./xml_uchime/temp"
    output_dir = "./xml_uchime/rescued_reads"

    try:
        process_all_xml_files(directory, temp_dir, output_dir, input_dir)
        print("Processing complete. Check the rescued_reads folder for final results.")
    except Exception as e:
        logging.error(f"An error occurred during processing: {e}")
