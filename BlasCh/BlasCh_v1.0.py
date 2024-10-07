"""
Chimera recovery module for Long-Read Sequencing

Description:
    This script processes BLAST XML results to identify and classify chimeric sequences in long-read sequencing data.
    It categorizes sequences as false positive chimeras, absolute chimeras, and borderline sequences
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
    - CSV file with sequence details, including query coverage, identity, and classification

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
Version: 1.0
Date: 2024-09-29
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
logging.basicConfig(level=logging.INFO, format='%(processName)s: %(message)s')

# Constants for classification criteria
HIGH_IDENTITY_THRESHOLD = 99.0
HIGH_COVERAGE_THRESHOLD = 99.0

# Dynamically determine the number of CPUs to use for multiprocessing
NUM_PROCESSES = max(1, multiprocessing.cpu_count()-1)

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

def analyze_blast_hits(blast_record, query_id):
    """Analyze BLAST hits for a given blast record."""
    hits_info = []
    self_hit = None
    for alignment in blast_record.alignments:
        if not alignment.hsps:
            continue  # Skip if there are no HSPs
        hsp = alignment.hsps[0]  # Consider only the best HSP for each alignment
        hit_id = extract_query_id(alignment.hit_def)
        
        # Skip self-hits (criterion 1)
        if hit_id == query_id:
            self_hit = {
                "hit_id": hit_id,
                "identity": (hsp.identities / hsp.align_length) * 100,
                "coverage": min((hsp.align_length / blast_record.query_length) * 100, 100),
                "is_same_sample": "size=" in hit_id
            }
            continue
        
        hits_info.append({
            "hit_id": hit_id,
            "identity": (hsp.identities / hsp.align_length) * 100,
            "coverage": min((hsp.align_length / blast_record.query_length) * 100, 100),
            "is_same_sample": "size=" in hit_id
        })
    
    return hits_info, self_hit

def classify_sequence(hits_info, self_hit):
    """Classify sequence based on hit information and criteria."""
    if not hits_info:
        return "non_chimeric", "No significant non-self hits"
    
    # Criterion 2: Certain false positive (non-chimeric)
    high_quality_db_hit = any(
        hit["identity"] >= HIGH_IDENTITY_THRESHOLD and 
        hit["coverage"] >= HIGH_COVERAGE_THRESHOLD and 
        not hit["is_same_sample"]
        for hit in hits_info
    )
    if high_quality_db_hit:
        return "non_chimeric", "High-quality match against database"
    
    # Criterion 3: Certain true positive chimeras (multiple alignments)
    if len(hits_info) > 1:
        return "chimeric", "Multiple non-self alignments"
    
    # Criterion 4: Middle-ground sequences (borderline)
    return "borderline", "Single alignment, requires further analysis"

def parse_blast_results(args):
    """Parse BLAST XML results and classify sequences into different categories."""
    xml_file, temp_dir, temp_2_dir, fasta_file = args
    logging.info(f"Starting processing for {xml_file}")
    
    try:
        fasta_sequences = load_fasta_sequences(fasta_file)
    except FileNotFoundError as e:
        logging.error(e)
        return set(), set(), set()
    except Exception as e:
        logging.error(f"Error loading FASTA file {fasta_file}: {e}")
        return set(), set(), set()

    non_chimeric_sequences = set()
    chimeric_sequences = set()
    borderline_sequences = set()

    sequence_details = []

    try:
        with open(xml_file) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for blast_record in blast_records:
                query_id = extract_query_id(blast_record.query)
                if query_id not in fasta_sequences:
                    continue

                hits_info, self_hit = analyze_blast_hits(blast_record, query_id)
                classification, reason = classify_sequence(hits_info, self_hit)
                
                if classification == "non_chimeric":
                    non_chimeric_sequences.add(query_id)
                elif classification == "chimeric":
                    chimeric_sequences.add(query_id)
                else:
                    borderline_sequences.add(query_id)
                
                logging.debug(f"{query_id} classified as {classification} ({reason})")
                
                # Add details for all hits
                for i, hit in enumerate(hits_info, 1):
                    sequence_details.append([
                        query_id, 
                        hit["coverage"], 
                        hit["identity"], 
                        classification,
                        f"Hit {i}",
                        "Same sample" if hit["is_same_sample"] else "Database"
                    ])
                
                # Add self-hit information if available
                if self_hit:
                    sequence_details.append([
                        query_id,
                        self_hit["coverage"],
                        self_hit["identity"],
                        "Self-hit",
                        "Self-hit",
                        "Same sample"
                    ])

                # Periodically write sequence details to prevent memory overflow
                if len(sequence_details) >= 1000:
                    write_sequence_details(sequence_details, temp_2_dir, fasta_file)
                    sequence_details.clear()

        # Write any remaining sequence details
        if sequence_details:
            write_sequence_details(sequence_details, temp_2_dir, fasta_file)
            sequence_details.clear()

        # Write sequences to output files only if they are non-empty
        base_filename = os.path.basename(fasta_file).replace(".chimeras.fasta", "")
        write_sequences_to_file(non_chimeric_sequences, fasta_sequences, os.path.join(temp_dir, f"{base_filename}_non_chimeric.fasta"))
        write_sequences_to_file(borderline_sequences, fasta_sequences, os.path.join(temp_dir, f"{base_filename}_borderline.fasta"))
        write_sequences_to_file(chimeric_sequences, fasta_sequences, os.path.join(temp_dir, f"{base_filename}_chimeric.fasta"))
    
        logging.info(f"Completed processing for {xml_file}")
        return non_chimeric_sequences, chimeric_sequences, borderline_sequences
    except Exception as e:
        logging.error(f"Error processing {xml_file}: {e}")
        return set(), set(), set()

def write_sequences_to_file(seq_ids, fasta_sequences, output_file):
    """Write sequences to a FASTA file only if seq_ids is not empty."""
    if not seq_ids:
        logging.info(f"No sequences to write for {output_file}. Skipping file creation.")
        return
    
    try:
        with open(output_file, 'w') as out_file:
            for seq_id in sorted(seq_ids):
                out_file.write(f">{seq_id}\n{fasta_sequences[seq_id]}\n")
        logging.info(f"Wrote {len(seq_ids)} sequences to {output_file}")
    except Exception as e:
        logging.error(f"Error writing to {output_file}: {e}")

def write_sequence_details(details, temp_2_dir, fasta_file):
    """Write sequence details to a CSV file in batches."""
    if not details:
        logging.debug("No sequence details to write. Skipping.")
        return
    
    try:
        base_filename = os.path.basename(fasta_file).replace(".chimeras.fasta", "")
        csv_file_path = os.path.join(temp_2_dir, f"{base_filename}_sequence_details.csv")
        file_exists = os.path.isfile(csv_file_path)
        with open(csv_file_path, 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            if not file_exists:
                csvwriter.writerow(["Sequence ID", "Query Coverage (%)", "Identity Percentage (%)", "Classification", "Hit Type", "Hit Origin"])
            csvwriter.writerows(details)
        logging.debug(f"Wrote {len(details)} sequence details to {csv_file_path}")
    except Exception as e:
        logging.error(f"Error writing sequence details to {csv_file_path}: {e}")

def generate_report(all_non_chimeric_sequences, all_absolute_chimeras, all_borderline_chimeras, file_results, output_dir):
    """Generate a detailed report summarizing the results, sorted by file name."""
    try:
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
                report_file.write(f"\nFile: {os.path.basename(file)}\n")
                file_total = sum(results.values())
                
                for category, count in results.items():
                    percentage = (count / file_total) * 100 if file_total > 0 else 0
                    report_file.write(f"  {category}: {count} ({percentage:.2f}%)\n")
                    
                report_file.write(f"  Total: {file_total}\n")
        logging.info("Report generated successfully.")
    except Exception as e:
        logging.error(f"Error generating report: {e}")

def process_all_xml_files(directory, temp_dir, temp_2_dir, output_dir, input_dir):
    """Process all BLAST XML files and generate final outputs."""
    # Clean up existing directories if they exist
    for dir_path in [temp_dir, temp_2_dir, output_dir]:
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)
        os.makedirs(dir_path, exist_ok=True)

    args_list = []
    file_results = defaultdict(lambda: defaultdict(int))
    
    all_non_chimeric_sequences = set()
    all_absolute_chimeras = set()
    all_borderline_chimeras = set()

    start_time = time.time()
    log_system_usage()

    # Prepare arguments for multiprocessing
    for filename in os.listdir(directory):
        if filename.endswith("_blast_results.xml"):
            xml_file = os.path.join(directory, filename)
            fasta_file = os.path.join(input_dir, filename.replace("_blast_results.xml", ".chimeras.fasta"))
            if os.path.exists(fasta_file):
                args_list.append((xml_file, temp_dir, temp_2_dir, fasta_file))
            else:
                logging.warning(f"FASTA file for {xml_file} not found. Skipping.")

    logging.info(f"Using {NUM_PROCESSES} processes for multiprocessing.")
    
    with multiprocessing.Pool(processes=NUM_PROCESSES) as pool:
        results = pool.map(parse_blast_results, args_list)

    # Merge all results
    for i, result in enumerate(results):
        non_chimeric_sequences, absolute_chimeras, borderline_chimeras = result
        all_non_chimeric_sequences.update(non_chimeric_sequences)
        all_absolute_chimeras.update(absolute_chimeras)
        all_borderline_chimeras.update(borderline_chimeras)
        
        xml_file = args_list[i][0]
        file_results[xml_file] = {
            "Non-Chimeric Sequences": len(non_chimeric_sequences),
            "Absolute Chimeras": len(absolute_chimeras),
            "Borderline Chimeras": len(borderline_chimeras)
        }

    # Combine non-chimeric and borderline sequences into final output
    for filename in os.listdir(input_dir):
        if filename.endswith(".chimeras.fasta"):
            base_filename = os.path.basename(filename).replace(".chimeras.fasta", "")
            non_chimeric_file = os.path.join(temp_dir, f"{base_filename}_non_chimeric.fasta")
            borderline_file = os.path.join(temp_dir, f"{base_filename}_borderline.fasta")
            chimeric_file = os.path.join(temp_dir, f"{base_filename}_chimeric.fasta")
            
            # Copy non-chimeric sequences if the file exists
            if os.path.exists(non_chimeric_file):
                shutil.copy(non_chimeric_file, os.path.join(output_dir, f"{base_filename}_non_chimeric.fasta"))
                logging.info(f"Copied {non_chimeric_file} to output directory.")
            else:
                logging.info(f"No non-chimeric sequences for {base_filename}. Skipping file copy.")
            
            # Copy borderline sequences if the file exists
            if os.path.exists(borderline_file):
                shutil.copy(borderline_file, os.path.join(output_dir, f"{base_filename}_borderline.fasta"))
                logging.info(f"Copied {borderline_file} to output directory.")
            else:
                logging.info(f"No borderline sequences for {base_filename}. Skipping file copy.")
            
            # Copy chimeric sequences if the file exists
            if os.path.exists(chimeric_file):
                shutil.copy(chimeric_file, os.path.join(output_dir, f"{base_filename}_chimeric.fasta"))
                logging.info(f"Copied {chimeric_file} to output directory.")
            else:
                logging.info(f"No chimeric sequences for {base_filename}. Skipping file copy.")

    generate_report(all_non_chimeric_sequences, all_absolute_chimeras, all_borderline_chimeras, file_results, output_dir)

    log_system_usage()

    end_time = time.time()
    elapsed_time = end_time - start_time
    logging.info(f"Total time taken: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    input_dir = "./input"      # Directory containing FASTA files
    directory = "."            # Directory containing XML files
    temp_dir = "./temp"
    temp_2_dir = "./temp_2"
    output_dir = "./rescued_reads"

    try:
        process_all_xml_files(directory, temp_dir, temp_2_dir, output_dir, input_dir)
        print("Processing complete. Check the rescued_reads folder for final results.")
    except Exception as e:
        logging.error(f"An error occurred during processing: {e}")
