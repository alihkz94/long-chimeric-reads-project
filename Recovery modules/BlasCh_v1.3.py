"""
Chimera recovery module for Long-Read Sequencing

UPDATES:
- Implemented multiprocessing for faster processing of multiple BLAST XML files.
- Added CPU and memory usage tracking using psutil.
- Improved error handling and logging with detailed process information.
- Optimized file I/O operations and memory management.
- Enhanced report generation with overall and per-file statistics.
- Implemented check for multiple alignments (HSPs) in the first hit to classify absolute chimeras.
- Refined classification criteria based on high-identity and high-coverage alignments.
- Added classification for borderline sequences based on adjusted identity and coverage thresholds.

Description:
    This script processes BLAST XML results to identify and classify chimeric sequences in long-read sequencing data.
    It categorizes sequences as false positive chimeras, absolute chimeras and borderline sequences
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

Author: ALI HAKIMZADEH
Version: 1.3
Date: 2024-09-24
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

# Configure logging to output process information
logging.basicConfig(level=logging.DEBUG, format='%(processName)s: %(message)s')

# Constants for classification criteria
HIGH_IDENTITY_THRESHOLD = 99.0  # Identity threshold to classify a hit as high-quality
HIGH_COVERAGE_THRESHOLD = 99.0  # Coverage threshold to classify a hit as high-quality

# Dynamically determine the number of CPUs to use for multiprocessing
NUM_PROCESSES = multiprocessing.cpu_count()

def log_system_usage():
    """Logs the current CPU and RAM usage for monitoring system resource consumption."""
    process = psutil.Process(os.getpid())
    cpu_percent = psutil.cpu_percent(interval=1)
    memory_info = process.memory_info()
    logging.info(f"CPU usage: {cpu_percent}%")
    logging.info(f"Memory usage: {memory_info.rss / (1024 ** 2):.2f} MB")

def load_fasta_sequences(fasta_file):
    """Load sequences from a FASTA file into a dictionary where the keys are sequence IDs."""
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
    
    fasta_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_sequences[record.id] = str(record.seq)
    return fasta_sequences

def extract_query_id(blast_query_def):
    """Extract the query ID from the BLAST query definition."""
    return blast_query_def  # Simply return the query definition as is (can be customized)

def parse_blast_results(args):
    """Parse BLAST XML results and classify sequences into different categories based on alignment details."""
    xml_file, temp_dir, temp_2_dir, fasta_file = args
    logging.debug(f"Starting processing for {xml_file}")
    
    # Load the FASTA sequences for later classification
    fasta_sequences = load_fasta_sequences(fasta_file)

    # Sets for storing classified sequences
    false_positive_chimeras = set()
    absolute_chimeras = set()
    borderline_sequences = set()

    # Store details for sequences (ID, coverage, identity) to write to CSV later
    sequence_details = []
    multiple_alignment_hits = []

    with open(xml_file) as result_handle:
        blast_records = list(NCBIXML.parse(result_handle))
        
        for blast_record in blast_records:
            query_id = extract_query_id(blast_record.query)
            if query_id not in fasta_sequences:
                continue  # Skip if the query is not in the FASTA file (shouldn't happen)

            high_identity_alignments = []
            all_alignments = []
            
            # Process BLAST alignments for the current query
            if blast_record.alignments:
                first_hit = blast_record.alignments[0]

                # Check if the first hit contains multiple high-scoring segment pairs (HSPs)
                if len(first_hit.hsps) > 1:
                    # Multiple HSPs usually indicate a chimera
                    absolute_chimeras.add(query_id)
                    logging.debug(f"{query_id} classified as absolute chimera due to multiple HSPs in the first hit")
                else:
                    # Examine the hit further to determine its origin
                    hit_id = extract_query_id(first_hit.hit_def)
                    if "size=" not in hit_id:
                        # Hit is from the database, indicating a potential chimera
                        false_positive_chimeras.add(query_id)
                        logging.debug(f"{query_id} classified as false positive chimera (hit is from database)")
                    else:
                        # Likely a self-hit (non-chimeric)
                        logging.debug(f"{query_id} has a high-identity hit from the same sample (non-chimeric)")

            # Check all alignments for multiple hits and high-identity alignments
            for alignment in blast_record.alignments:
                hit_id = extract_query_id(alignment.hit_def)
                if hit_id != query_id:  # Exclude self-hits
                    for hsp in alignment.hsps:
                        # Calculate query coverage and identity percentage
                        query_coverage = min((hsp.align_length / blast_record.query_length) * 100, 100)
                        identity_percentage = (hsp.identities / hsp.align_length) * 100
                        
                        # Record sequence details for CSV output
                        sequence_details.append([query_id, query_coverage, identity_percentage])
                        all_alignments.append((alignment, hsp, query_coverage, identity_percentage))

                        if identity_percentage >= HIGH_IDENTITY_THRESHOLD and query_coverage >= HIGH_COVERAGE_THRESHOLD:
                            high_identity_alignments.append((alignment, hsp))

            # Further classification based on alignment quality
            if high_identity_alignments:
                if len(high_identity_alignments) > 1:
                    # Multiple high-identity alignments indicate a chimera
                    absolute_chimeras.add(query_id)
                    logging.debug(f"{query_id} classified as absolute chimera due to multiple high-identity alignments")
                else:
                    hit_id = extract_query_id(high_identity_alignments[0][0].hit_def)
                    if "size=" not in hit_id:
                        false_positive_chimeras.add(query_id)
                        logging.debug(f"{query_id} classified as false positive chimera (hit is from database)")
                    else:
                        logging.debug(f"{query_id} has a high-identity hit from the same sample (non-chimeric)")
            elif len(all_alignments) > 1:
                # If there are multiple alignments but not high identity, classify as borderline or absolute chimera
                max_identity = max(align[3] for align in all_alignments)
                max_coverage = max(align[2] for align in all_alignments)
                if max_identity >= 80.0 and max_coverage >= 91.0:
                    borderline_sequences.add(query_id)
                    logging.debug(f"{query_id} classified as borderline sequence")
                else:
                    absolute_chimeras.add(query_id)
                    logging.debug(f"{query_id} classified as absolute chimera")

        # Write intermediate results to the temporary folder
        base_filename = os.path.basename(fasta_file).replace(".chimeras.fasta", "")
        write_sequences_to_file(false_positive_chimeras, fasta_sequences, os.path.join(temp_dir, f"{base_filename}_false_positive_chimeras.fasta"))
        write_sequences_to_file(borderline_sequences, fasta_sequences, os.path.join(temp_dir, f"{base_filename}_borderline_sequences.fasta"))
    
        # Write sequence details (ID, coverage, identity) to CSV
        csv_file_path = os.path.join(temp_2_dir, f"{base_filename}_sequence_details.csv")
        with open(csv_file_path, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(["Sequence ID", "Query Coverage (%)", "Identity Percentage (%)", "Classification"])
            for detail in sequence_details:
                if detail[0] in false_positive_chimeras:
                    classification = "False Positive Chimera"
                elif detail[0] in absolute_chimeras:
                    classification = "Absolute Chimera"
                elif detail[0] in borderline_sequences:
                    classification = "Borderline"
                else:
                    classification = "Unknown"
                csvwriter.writerow(detail + [classification])

    logging.debug(f"Completed processing for {xml_file}")
    return false_positive_chimeras, absolute_chimeras, borderline_sequences

def write_sequences_to_file(seq_ids, fasta_sequences, output_file):
    """Write sequences with specified IDs to a FASTA file."""
    with open(output_file, 'w') as out_file:
        for seq_id in sorted(seq_ids):
            out_file.write(f">{seq_id}\n{fasta_sequences[seq_id]}\n")

def combine_sequences(temp_dir, base_filename, output_dir):
    """Combine false-positive chimera and borderline sequences into a single final output file."""
    false_positive_file = os.path.join(temp_dir, f"{base_filename}_false_positive_chimeras.fasta")
    borderline_file = os.path.join(temp_dir, f"{base_filename}_borderline_sequences.fasta")
    recovered_file = os.path.join(output_dir, f"{base_filename}_recovered.fasta")

    with open(recovered_file, 'w') as outfile:
        for fname in [false_positive_file, borderline_file]:
            if os.path.exists(fname):
                with open(fname) as infile:
                    shutil.copyfileobj(infile, outfile)

def generate_report(all_false_positive_chimeras, all_absolute_chimeras, all_borderline_sequences, file_results, output_dir):
    """Generate a detailed report summarizing the results by file."""
    report = {
        "False Positive Chimeras": len(all_false_positive_chimeras),
        "Absolute Chimeras": len(all_absolute_chimeras),
        "Borderline Sequences": len(all_borderline_sequences)
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
            
            recovered_sequences = results["False Positive Chimeras"] + results["Borderline Sequences"]
            
            for category, count in results.items():
                percentage = (count / file_total) * 100 if file_total > 0 else 0
                report_file.write(f"  {category}: {count} ({percentage:.2f}%)\n")
                
            report_file.write(f"  Total: {file_total}\n")
            report_file.write(f"  Total Recovered Sequences: {recovered_sequences}\n")

def process_all_xml_files(directory, temp_dir, temp_2_dir, output_dir, input_dir):
    """Process all BLAST XML files in the directory and generate final outputs."""
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    if os.path.exists(temp_2_dir):
        shutil.rmtree(temp_2_dir)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(temp_2_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    args_list = []
    file_results = defaultdict(lambda: defaultdict(int))
    
    all_false_positive_chimeras = set()
    all_absolute_chimeras = set()
    all_borderline_sequences = set()

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
        false_positive_chimeras, absolute_chimeras, borderline_sequences = result
        all_false_positive_chimeras.update(false_positive_chimeras)
        all_absolute_chimeras.update(absolute_chimeras)
        all_borderline_sequences.update(borderline_sequences)
        
        file_results[args_list[i][0]] = {
            "False Positive Chimeras": len(false_positive_chimeras),
            "Absolute Chimeras": len(absolute_chimeras),
            "Borderline Sequences": len(borderline_sequences)
        }

    # Combine sequences and generate final report
    for filename in os.listdir(input_dir):
        if filename.endswith(".chimeras.fasta"):
            base_filename = os.path.basename(filename).replace(".chimeras.fasta", "")
            combine_sequences(temp_dir, base_filename, output_dir)
    
    generate_report(all_false_positive_chimeras, all_absolute_chimeras, all_borderline_sequences, file_results, output_dir)

    log_system_usage()  # Log final system usage

    end_time = time.time()
    elapsed_time = end_time - start_time
    logging.info(f"Total time taken: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    # Set input and output directories for processing
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
