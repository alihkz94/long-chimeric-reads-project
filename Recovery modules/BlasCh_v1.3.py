"""
Chimera recovery module for Long-Read Sequencing

UPDATES:
- Implemented multiprocessing for faster processing of multiple BLAST XML files.
- Added CPU and memory usage tracking using psutil.
- Improved error handling and logging with detailed process information.
- Optimized file I/O operations and memory management.
- Enhanced report generation with overall and per-file statistics.
- Implemented a check for multiple alignments (HSPs) in the first hit to classify absolute chimeras.
- Refined classification criteria based on high-identity and high-coverage alignments.
- Classification for borderline sequences was added based on adjusted identity and coverage thresholds.

Description:
    This script processes BLAST XML results to identify and classify chimeric sequences in long-read sequencing data.
    It categorizes sequences as non-chimeric sequences, absolute chimeras, and borderline sequences
    based on specified alignment criteria. The script utilizes multiprocessing to efficiently handle large datasets
    and includes monitoring system resource usage.

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
    - csv

Input:
    - FASTA files in the specified input directory
    - BLAST XML result files in the specified directory

Output:
    - Classified sequences in separate FASTA files within the output directory
    - Detailed report file summarizing overall and per-file results
    - Log file with process information and system resource usage
    - CSV file with sequence details including query coverage, identity, and classification

Configuration:
    Modify the following variables at the beginning of the script to customize paths and thresholds:
    - input_dir: Directory containing input FASTA files
    - directory: Directory containing BLAST XML files (current directory by default)
    - temp_dir: Directory for storing temporary files
    - output_dir: Directory for storing output files
    - NUM_PROCESSES: Number of processes to use for multiprocessing (automatically set to CPU count)
    - HIGH_IDENTITY_THRESHOLD: Threshold for high identity percentage (99.0)
    - HIGH_COVERAGE_THRESHOLD: Threshold for high query coverage (99.0)
    - BORDERLINE_COVERAGE_THRESHOLD: Threshold for borderline coverage (80.0)
    - BORDERLINE_IDENTITY_THRESHOLD: Threshold for borderline identity (91.0)

Author: ALI HAKIMZADEH
Version: 1.3
Date: 2024-09-25
"""

# Load Libraries 
import os
import shutil
import multiprocessing
import time
import psutil
from Bio.Blast import NCBIXML
from Bio import SeqIO
from collections import defaultdict
import logging
import csv

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(processName)s: %(message)s')

# Constants for classification criteria
HIGH_IDENTITY_THRESHOLD = 99.0
HIGH_COVERAGE_THRESHOLD = 99.0
BORDERLINE_COVERAGE_THRESHOLD = 80.0
BORDERLINE_IDENTITY_THRESHOLD = 91.0

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
    xml_file, temp_dir, temp_2_dir, fasta_file = args
    logging.debug(f"Starting processing for {xml_file}")
    
    fasta_sequences = load_fasta_sequences(fasta_file)

    non_chimeric_sequences = set()
    absolute_chimeras = set()
    borderline_chimeras = set()

    sequence_details = []

    with open(xml_file) as result_handle:
        blast_records = list(NCBIXML.parse(result_handle))
        
        for blast_record in blast_records:
            query_id = extract_query_id(blast_record.query)
            if query_id not in fasta_sequences:
                continue

            if blast_record.alignments:
                first_hit = blast_record.alignments[0]
                first_hsp = first_hit.hsps[0]
                
                query_coverage = min((first_hsp.align_length / blast_record.query_length) * 100, 100)
                identity_percentage = (first_hsp.identities / first_hsp.align_length) * 100
                
                sequence_details.append([query_id, query_coverage, identity_percentage])

                # Check if the first hit contains multiple high-scoring segment pairs (HSPs)
                if len(first_hit.hsps) > 1:
                    # Multiple HSPs usually indicate a chimera
                    absolute_chimeras.add(query_id)
                    logging.debug(f"{query_id} classified as absolute chimera due to multiple HSPs in the first hit")
                elif identity_percentage >= HIGH_IDENTITY_THRESHOLD and query_coverage >= HIGH_COVERAGE_THRESHOLD:
                    hit_id = extract_query_id(first_hit.hit_def)
                    if "size=" not in hit_id:
                        # Hit is from the database, indicating a potential chimera
                        absolute_chimeras.add(query_id)
                        logging.debug(f"{query_id} classified as absolute chimera (high-identity hit is from database)")
                    else:
                        # Likely a self-hit (non-chimeric)
                        non_chimeric_sequences.add(query_id)
                        logging.debug(f"{query_id} classified as non-chimeric (high-identity hit from same sample)")
                elif query_coverage >= BORDERLINE_COVERAGE_THRESHOLD and identity_percentage >= BORDERLINE_IDENTITY_THRESHOLD:
                    borderline_chimeras.add(query_id)
                    logging.debug(f"{query_id} classified as borderline chimera")
                else:
                    absolute_chimeras.add(query_id)
                    logging.debug(f"{query_id} classified as absolute chimera (low identity/coverage)")
            else:
                non_chimeric_sequences.add(query_id)
                logging.debug(f"{query_id} classified as non-chimeric (no significant hits)")

        # Write intermediate results to the temporary folder
        base_filename = os.path.basename(fasta_file).replace(".chimeras.fasta", "")
        write_sequences_to_file(non_chimeric_sequences, fasta_sequences, os.path.join(temp_dir, f"{base_filename}_non_chimeric.fasta"))
    
        # Write the sequence details (ID, coverage, identity) to CSV in temp_2
        csv_file_path = os.path.join(temp_2_dir, f"{base_filename}_sequence_details.csv")
        with open(csv_file_path, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(["Sequence ID", "Query Coverage (%)", "Identity Percentage (%)"])
            csvwriter.writerows(sequence_details)

    logging.debug(f"Completed processing for {xml_file}")
    return non_chimeric_sequences, absolute_chimeras, borderline_chimeras

def write_sequences_to_file(seq_ids, fasta_sequences, output_file):
    """Write sequences to a FASTA file."""
    with open(output_file, 'w') as out_file:
        for seq_id in sorted(seq_ids):
            out_file.write(f">{seq_id}\n{fasta_sequences[seq_id]}\n")

def generate_report(all_non_chimeric_sequences, all_absolute_chimeras, all_borderline_chimeras, file_results, output_dir):
    """Generate a detailed report summarizing the results, sorted by file name."""
    report = {
        "Non-Chimeric Sequences": len(all_non_chimeric_sequences),
        "Absolute Chimeras": len(all_absolute_chimeras),
        "Borderline Chimeras": len(all_borderline_chimeras)
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
        
        sorted_file_results = dict(sorted(file_results.items()))
        
        for file, results in sorted_file_results.items():
            report_file.write(f"\nFile: {file}\n")
            file_total = sum(results.values())
            
            for category, count in results.items():
                percentage = (count / file_total) * 100 if file_total > 0 else 0
                report_file.write(f"  {category}: {count} ({percentage:.2f}%)\n")
                
            report_file.write(f"  Total: {file_total}\n")

def process_all_xml_files(directory, temp_dir, temp_2_dir, output_dir, input_dir):
    """Process all BLAST XML files and generate final outputs."""
    # Clean up existing directories if they exist
    for dir in [temp_dir, temp_2_dir, output_dir]:
        if os.path.exists(dir):
            shutil.rmtree(dir)
        os.makedirs(dir, exist_ok=True)

    args_list = []
    file_results = defaultdict(lambda: defaultdict(int))
    
    all_non_chimeric_sequences = set()
    all_absolute_chimeras = set()
    all_borderline_chimeras = set()

    start_time = time.time()
    log_system_usage()

    for filename in os.listdir(directory):
        if filename.endswith("_blast_results.xml"):
            xml_file = os.path.join(directory, filename)
            fasta_file = os.path.join(input_dir, filename.replace("_blast_results.xml", ".chimeras.fasta"))
            args_list.append((xml_file, temp_dir, temp_2_dir, fasta_file))

    logging.debug(f"Using {NUM_PROCESSES} processes for multiprocessing")
    
    with multiprocessing.Pool(processes=NUM_PROCESSES) as pool:
        results = pool.map(parse_blast_results, args_list)

    for i, result in enumerate(results):
        non_chimeric_sequences, absolute_chimeras, borderline_chimeras = result
        all_non_chimeric_sequences.update(non_chimeric_sequences)
        all_absolute_chimeras.update(absolute_chimeras)
        all_borderline_chimeras.update(borderline_chimeras)
        
        file_results[args_list[i][0]] = {
            "Non-Chimeric Sequences": len(non_chimeric_sequences),
            "Absolute Chimeras": len(absolute_chimeras),
            "Borderline Chimeras": len(borderline_chimeras)
        }

    # Combine non-chimeric sequences into final output
    for filename in os.listdir(input_dir):
        if filename.endswith(".chimeras.fasta"):
            base_filename = os.path.basename(filename).replace(".chimeras.fasta", "")
            non_chimeric_file = os.path.join(temp_dir, f"{base_filename}_non_chimeric.fasta")
            if os.path.exists(non_chimeric_file):
                shutil.copy(non_chimeric_file, os.path.join(output_dir, f"{base_filename}_recovered.fasta"))
    
    generate_report(all_non_chimeric_sequences, all_absolute_chimeras, all_borderline_chimeras, file_results, output_dir)

    log_system_usage()

    end_time = time.time()
    elapsed_time = end_time - start_time
    logging.info(f"Total time taken: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    input_dir = "./input"
    directory = "."
    temp_dir = "./temp"
    temp_2_dir = "./temp_2"
    output_dir = "./rescued_reads"

    try:
        process_all_xml_files(directory, temp_dir, temp_2_dir, output_dir, input_dir)
        print("Processing complete. Check the rescued_reads folder for final results.")
    except Exception as e:
        logging.error(f"An error occurred during processing: {e}")
