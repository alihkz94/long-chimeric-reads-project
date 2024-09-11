"""
False positive chimeras recovery for long-read sequencing

UPDATES: 
-This version includes reporting the number of rescued sequences and the time the module takes to finish. 
-It handles CPU management better by running one CPU for each chunk for BLAST.
-BLASTING files now in alphabetical order 

Description:
    This script is designed to recover false positive chimeric sequences from a set of FASTA files based on BLAST results.
    It processes input FASTA files, performs BLAST searches, and identifies non-chimeric sequences based on
    specified criteria. The script then saves the non-chimeric sequences and performs additional analysis
    on trimmed sequences.

Usage:
    python BlasCh_v1.py

Requirements:
    - Python 3.6+
    - BioPython
    - BLAST+ (Blastn command-line tool)

Dependencies:
    - os
    - logging
    - shutil
    - subprocess
    - multiprocessing
    - Bio (from BioPython)
    - time

Input:
    - FASTA files in the specified input directory
    - BLAST result files in the specified blast output directory
    - BLAST database specified by the 'db' variable

Output:
    - Filtered non-chimeric sequences in the 'rescued_reads' directory
    - Trimmed sequences and their BLAST results in 'begin' and 'end' directories

Configuration:
    Modify the following variables at the beginning of the script to customize paths and thresholds:
    - input_dir: Directory containing input FASTA files
    - blast_output_dir: Directory containing BLAST result files
    - output_dir_begin: Directory for storing trimmed sequences (first 100 bases)
    - output_dir_end: Directory for storing trimmed sequences (last 100 bases)
    - rescued_dir: Directory for storing non-chimeric sequences
    - db: Path to the BLAST database
    - header: BLAST output format specification

Author: Ali Hakimzadeh
Version: 0.1.0
Date: 2024-08-28
"""

# Import necessary libraries
import os
import logging
import shutil
import subprocess
import multiprocessing as mp
from Bio import SeqIO
import argparse
import time

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Argument parser for user-specified options
def parse_args():
    parser = argparse.ArgumentParser(description='Chimeric Sequence Rescue Script')
    parser.add_argument('--cpus', type=int, default=8, help='Number of CPUs to use for BLAST (default: 8)')
    return parser.parse_args()

# Define directory names
input_dir = 'tmp'
blast_output_dir = 'blast_tmp'
output_dir_begin = 'begin'
output_dir_end = 'end'
rescued_dir = 'false_positive_chimeras'
db = 'database/EUK'
header = "qseqid stitle qlen slen qstart qend sstart send evalue length nident mismatch gapopen gaps sstrand qcovs pident"
report_file = 'rescue_report.txt'

# Step 1: Clean and create necessary directories
def clean_directory(directory):
    if os.path.exists(directory):
        shutil.rmtree(directory)  # Delete the directory and its contents
    os.makedirs(directory, exist_ok=True)  # Recreate the directory

# Clean and create directories
clean_directory(output_dir_begin)
clean_directory(output_dir_end)
clean_directory(rescued_dir)

# Step 2: Check if the FASTA and BLAST directories contain the same files
def check_file_pairs(fasta_dir, blast_dir):
    fasta_files = set(os.path.splitext(f)[0] for f in os.listdir(fasta_dir) if f.endswith(".fasta"))
    blast_files = set(os.path.splitext(f)[0] for f in os.listdir(blast_dir) if f.endswith(".txt"))

    missing_blasts = fasta_files - blast_files
    missing_fastas = blast_files - fasta_files

    if missing_blasts or missing_fastas:
        logging.warning(f"Missing BLAST results for the following FASTA files: {', '.join(missing_blasts)}")
        logging.warning(f"Missing FASTA files for the following BLAST results: {', '.join(missing_fastas)}")
        
        # Stop execution and raise an error
        raise ValueError("Mismatched FASTA and BLAST files detected. Stopping execution.")

    # Return only the files that have pairs
    valid_files = fasta_files & blast_files
    return valid_files

# Step 3: Filter sequences based on BLAST results with error handling
def filter_sequences(blast_file_path, qcov_threshold=99, pident_threshold=99):
    nonchimeric_sequences = []
    chimeric_candidates = []
    
    with open(blast_file_path, 'r') as file:
        lines = file.readlines()
        if not lines:
            logging.warning(f"No BLAST hits found in file: {blast_file_path}")
            return nonchimeric_sequences, chimeric_candidates

        for line in lines:
            # Skip header lines
            if line.startswith("qseqid"):
                continue
            
            parts = line.strip().split('+')
            if len(parts) < 17:
                logging.warning(f"Skipping line due to unexpected format (less than 17 fields): {line.strip()}")
                continue

            qseqid = parts[0]
            try:
                qcov = float(parts[15])
                pident = float(parts[16])
            except ValueError as e:
                logging.error(f"Error parsing qcov or pident in line: {line.strip()} - {e}")
                continue
            
            if qcov >= qcov_threshold and pident >= pident_threshold:
                nonchimeric_sequences.append(qseqid)
            else:
                chimeric_candidates.append(qseqid)  # Mark these sequences for further analysis
    
    return nonchimeric_sequences, chimeric_candidates

# Step 4: Save nonchimeric sequences to the "false_positive_chimeras" folder without line breaks
def save_nonchimeric_sequences(input_file, nonchimeric_sequences, output_dir):
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    nonchimeric_file = os.path.join(output_dir, f"{base_name}_rescued.fasta")

    with open(input_file) as handle:
        records = [record for record in SeqIO.parse(handle, "fasta") if record.id in nonchimeric_sequences]

    if records:  # Only write if there are nonchimeric sequences
        with open(nonchimeric_file, 'w') as output_handle:
            for record in records:
                output_handle.write(f">{record.id}\n{str(record.seq)}\n")
    
    return len(records)  # Return the number of rescued sequences

# Step 5: Trim sequences and save to separate files
def trim_sequences(input_file, output_dir_begin, output_dir_end, filtered_sequences):
    base_name = os.path.splitext(os.path.basename(input_file))[0]

    with open(input_file) as handle:
        records = [record for record in SeqIO.parse(handle, "fasta") if record.id in filtered_sequences]

    trimmed_begin = []
    trimmed_end = []

    for record in records:
        if len(record.seq) < 2:
            continue  # Skip sequences that are too short to process

        # Trim the first 100 bases
        trimmed_begin.append(record[:100] if len(record.seq) >= 100 else record)

        # Trim the last 100 bases
        trimmed_end.append(record[-100:] if len(record.seq) >= 100 else record)

    if trimmed_begin:
        SeqIO.write(trimmed_begin, os.path.join(output_dir_begin, f"{base_name}_begin.fasta"), "fasta")
    if trimmed_end:
        SeqIO.write(trimmed_end, os.path.join(output_dir_end, f"{base_name}_end.fasta"), "fasta")

# Step 6: Run BLAST on trimmed sequences and clean up chunks
def run_blast(fasta_file, db, header, num_cpus):
    output_dir = os.path.dirname(fasta_file)
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    total_seqs = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))

    if total_seqs == 0:
        logging.warning(f"BLAST skipped: No sequences found in {fasta_file}")
        return

    num_chunks = num_cpus  # Set number of chunks to match the number of CPUs
    seq_per_chunk = total_seqs // num_chunks if total_seqs >= num_chunks else 1

    # Split the FASTA file into chunks
    chunk_files = []
    with open(fasta_file) as f:
        records = list(SeqIO.parse(f, "fasta"))
        for i in range(num_chunks):
            chunk_records = records[i * seq_per_chunk:(i + 1) * seq_per_chunk]
            chunk_file = f"{output_dir}/{base_name}_chunk_{i + 1}.fasta"
            SeqIO.write(chunk_records, chunk_file, "fasta")
            chunk_files.append(chunk_file)
    
    # Run BLAST for each chunk in parallel
    processes = []
    for i in range(num_chunks):
        chunk_file = chunk_files[i]
        blast_output = f"{chunk_file}_blast_results.txt"

        # Ensure the chunk file has sequences
        if not os.path.exists(chunk_file) or os.path.getsize(chunk_file) == 0:
            logging.warning(f"Skipping BLAST for empty chunk file: {chunk_file}")
            continue

        p = subprocess.Popen([
            'blastn',
            '-query', chunk_file,
            '-db', db,
            '-word_size', '7',
            '-task', 'blastn',
            '-num_threads', '1',  # Use 1 thread per chunk
            '-outfmt', f'6 delim=+ {header}',
            '-evalue', '0.001',
            '-strand', 'both',
            '-max_target_seqs', '10',
            '-max_hsps', '1',
            '-out', blast_output
        ])
        processes.append(p)
    
    # Wait for all processes to complete
    for p in processes:
        p.wait()

    # Combine and deduplicate the results
    combined_output = f"{output_dir}/combined_blast_top10hit.txt"
    with open(combined_output, 'w') as outfile:
        for i in range(num_chunks):
            chunk_file = chunk_files[i]
            chunk_blast_output = f"{chunk_file}_blast_results.txt"
            if os.path.exists(chunk_blast_output):
                with open(chunk_blast_output) as infile:
                    outfile.write(infile.read())
                os.remove(chunk_blast_output)  # Remove the individual BLAST output chunk files

    dedup_output = f"{output_dir}/{base_name}.txt"
    if os.path.exists(combined_output):
        with open(combined_output) as infile, open(dedup_output, 'w') as outfile:
            seen = set()
            for line in infile:
                if line.split('+')[0] not in seen:
                    seen.add(line.split('+')[0])
                    outfile.write(line)
        os.remove(combined_output)  # Remove the combined output file

    # Cleanup: Remove chunk files after processing
    for chunk_file in chunk_files:
        if os.path.exists(chunk_file):
            os.remove(chunk_file)

# Step 7: Check database identifiers match between begin and end BLAST results
def check_database_identifier(begin_blast_line, end_blast_line):
    begin_parts = begin_blast_line.split('+')
    end_parts = end_blast_line.split('+')
    
    # Extract database identifiers
    begin_db_id = begin_parts[1]
    end_db_id = end_parts[1]
    
    # Check if they match
    return begin_db_id == end_db_id

# Step 8: Criteria function to determine if a sequence is non-chimeric
def check_criteria(qlen, slen, qstart, qend, sstart, send, qcov, pident):
    if qcov >= 95 and pident >= 95:
        return True  # Non-chimeric
    else:
        return False  # Chimeric

# Step 9: Parse the BLAST output file for a specific sequence
def parse_blast_file(blast_file_path, qseqid):
    with open(blast_file_path, 'r') as file:
        lines = file.readlines()
        if not lines:
            logging.warning(f"No BLAST hits found for {qseqid} in file: {blast_file_path}")
            return None  # No BLAST hits, consider as chimeric
        
        for line in lines:
            if line.startswith("qseqid"):
                continue
            
            if line.startswith(qseqid):
                parts = line.strip().split('+')
                if len(parts) < 17:
                    logging.warning(f"Skipping line due to unexpected format (less than 17 fields): {line.strip()}")
                    continue
                
                qlen = int(parts[2])
                slen = int(parts[3])
                qstart = int(parts[4])
                qend = int(parts[5])
                sstart = int(parts[6])
                send = int(parts[7])
                qcov = float(parts[15])
                pident = float(parts[16])

                return qlen, slen, qstart, qend, sstart, send, qcov, pident

    return None  # Return None if no valid data is found or no BLAST hits

# Step 10: Process each qseqid and determine chimeric/non-chimeric status
def process_qseqid(qseqid):
    if isinstance(qseqid, list):
        qseqid = qseqid[0]  # Extract the first element if it's a list
    
    base_name = qseqid.split('.')[0]
    blast_file_begin = os.path.join(output_dir_begin, f"{base_name}_begin.txt")
    blast_file_end = os.path.join(output_dir_end, f"{base_name}_end.txt")

    begin_data = parse_blast_file(blast_file_begin, qseqid)
    end_data = parse_blast_file(blast_file_end, qseqid)

    if not begin_data or not end_data:
        return None  # Skip this sequence if either begin or end data is missing

    begin_line = next((line for line in open(blast_file_begin) if line.startswith(qseqid)), None)
    end_line = next((line for line in open(blast_file_end) if line.startswith(qseqid)), None)
    
    if begin_line and end_line and not check_database_identifier(begin_line, end_line):
        return None  # Skip this sequence if database identifiers do not match

    non_chimeric_b = check_criteria(*begin_data)
    non_chimeric_e = check_criteria(*end_data)

    if non_chimeric_b and non_chimeric_e:
        return qseqid  # Return the sequence ID if it is non-chimeric
    else:
        return None  # Otherwise, return None

# Step 11: Main function to run the script
def main():
    start_time = time.time()  # Start the timer
    args = parse_args()

    print("Step 1: Cleaning directories...")
    clean_directory(output_dir_begin)
    clean_directory(output_dir_end)
    clean_directory(rescued_dir)

    print("Step 2: Checking for matching FASTA and BLAST result files...")
    valid_files = check_file_pairs(input_dir, blast_output_dir)
    
    rescued_sequences_count = []

    print("Step 3: Filtering sequences based on BLAST results...")
    for base_name in valid_files:
        blast_file = os.path.join(blast_output_dir, f"{base_name}.txt")
        input_file = os.path.join(input_dir, f"{base_name}.fasta")
        
        # Step 3: Filtering sequences based on BLAST results
        nonchimeric_sequences, chimeric_candidates = filter_sequences(blast_file)

        # Step 4: Save non-chimeric sequences immediately
        if nonchimeric_sequences:
            rescued_count = save_nonchimeric_sequences(input_file, nonchimeric_sequences, rescued_dir)
            rescued_sequences_count.append((base_name, rescued_count))

        # Step 5: Process chimeric candidates for further analysis
        if chimeric_candidates:
            trim_sequences(input_file, output_dir_begin, output_dir_end, chimeric_candidates)
    
    print("Step 6: Running BLAST for begin and end directories...")
    for file_name in sorted(os.listdir(output_dir_begin)):  # Sort the files alphabetically
        if file_name.endswith(".fasta"):
            print(f"Running BLAST for {file_name} in begin directory...")
            run_blast(os.path.join(output_dir_begin, file_name), db, header, args.cpus)
    
    for file_name in sorted(os.listdir(output_dir_end)):  # Sort the files alphabetically
        if file_name.endswith(".fasta"):
            print(f"Running BLAST for {file_name} in end directory...")
            run_blast(os.path.join(output_dir_end, file_name), db, header, args.cpus)
    
    print("Step 7: Processing and identifying chimeric sequences...")
    all_nonchimeric_sequences = []

    for base_name in sorted(valid_files):  # Sort the files alphabetically
        blast_file = os.path.join(blast_output_dir, f"{base_name}.txt")
        input_file = os.path.join(input_dir, f"{base_name}.fasta")
        nonchimeric_sequences, _ = filter_sequences(blast_file)

        for qseqid in nonchimeric_sequences:
            result = process_qseqid(qseqid)
            if result:
                all_nonchimeric_sequences.append(result)

    # Generate and save the detailed report
    with open(report_file, 'w') as report:
        report.write("Base Name\tRescued Sequences\n")
        for base_name, count in rescued_sequences_count:
            report.write(f"{base_name}.fasta\t{count}\n")
        total_rescued = sum(count for _, count in rescued_sequences_count)
        report.write(f"Total Rescued Sequences: {total_rescued}\n")

    # Report the total time taken
    elapsed_time = time.time() - start_time
    print(f"Non-chimeric sequences have been rescued and saved in {rescued_dir}")
    print(f"Detailed report saved as {report_file}")
    print(f"Script completed in {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    main()
