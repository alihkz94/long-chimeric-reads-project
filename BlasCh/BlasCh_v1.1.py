"""
Chimera Detection and Recovery Module for Metabarcoding and Environmental DNA (eDNA) Analysis

Description:
    This script processes BLAST XML results to identify, classify, and recover False positive chimeric sequences 
    from metabarcoding or eDNA datasets. Sequences are categorized into three groups: non-chimeric, absolute 
    chimeras, and borderline sequences, based on identity and coverage thresholds.

Usage:
    python chimera_recovery.py

Requirements:
    - Python 3.6+
    - BioPython
    - psutil

Dependencies:
    - os: for file and directory operations
    - shutil: to copy and remove directories and files
    - multiprocessing: for parallel processing to speed up analysis
    - time: to track execution time
    - psutil: to monitor CPU and memory usage
    - Bio.Blast.NCBIXML (from BioPython): to parse BLAST XML outputs
    - Bio.SeqIO (from BioPython): to load and process FASTA files
    - collections (defaultdict): for efficient grouping of taxonomy data
    - logging: to provide process monitoring and error reporting
    - csv: to write detailed sequence classifications to CSV files

Input:
    - FASTA files containing sequences to be evaluated are located in the specified input directory.
    - BLAST XML files with alignment results are in the working directory (or a specified directory).

Output:
    - **FASTA Files**: Classified sequences (non-chimeric, borderline) stored in the output directory.
    - **CSV Files**: Sequence details, including identity, query coverage, classification, and taxonomy.
    - **Text Report**: A summary report in the output directory, showing the overall sequence classification.
    - **Log Files**: Process logs containing system resource usage, file processing details, and error messages.

Configuration:
    The following variables can be adjusted to customize the script’s behavior:
    - `input_dir`: Directory containing the input FASTA files.
    - `directory`: Directory containing BLAST XML files (default: current working directory).
    - `temp_dir`: Directory to store temporary files.
    - `temp_2_dir`: Directory to store intermediate files (removed after execution).
    - `output_dir`: Directory for final output files.
    - `NUM_PROCESSES`: Number of CPU cores to use (default: all available CPUs minus one).
    - `HIGH_IDENTITY_THRESHOLD`: Threshold for high identity percentage (default: 99.0).
    - `HIGH_COVERAGE_THRESHOLD`: Threshold for high query coverage (default: 99.0).

Author: ALI HAKIMZADEH  
Version: 1.1  
Date: 2024-10-23
"""

#load libraries
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
NUM_PROCESSES = max(1, multiprocessing.cpu_count() - 1)  

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

def is_self_hit(hit_def):
    """
    Determine if a hit is a self-hit based on the presence of ';size=' in hit_def.
    """
    return ';size=' in hit_def

def extract_taxonomy(hit_def):
    """
    Extract taxonomy information from Hit_def.
    Assumes taxonomy starts after the first semicolon.
    """
    try:
        taxonomy = hit_def.split(';', 1)[1]
    except IndexError:
        taxonomy = "Unclassified"
    return taxonomy

def analyze_blast_hits(blast_record, query_id):
    """Analyze BLAST hits for a given blast record by processing all HSPs."""
    hits_info = []
    
    # Check if there are any alignments
    if not blast_record.alignments:
        return hits_info
        
    # Process first hit and check if it's a self-hit
    first_alignment = blast_record.alignments[0]
    first_hit_id = extract_query_id(first_alignment.hit_def)
    
    # Determine which hit to check for multiple HSPs
    hit_to_check = first_alignment
    if first_hit_id == query_id and len(blast_record.alignments) > 1:
        # If first hit is self-hit and there's a second hit, use that instead
        hit_to_check = blast_record.alignments[1]
        hit_id = extract_query_id(hit_to_check.hit_def)
    else:
        hit_id = first_hit_id

    # Check for multiple HSPs in the selected hit
    if hit_id != query_id and len(hit_to_check.hsps) > 1:
        return [{
            "hit_id": hit_id,
            "identity": 100,  # Placeholder values
            "coverage": 100,
            "is_same_sample": False,
            "taxonomy": extract_taxonomy(hit_to_check.hit_def),
            "force_chimeric": True,
            "multiple_hsps": True
        }]
    
    # If no multiple HSPs, process all hits normally
    for alignment in blast_record.alignments:
        if not alignment.hsps:
            continue
            
        hit_def = alignment.hit_def
        hit_id = extract_query_id(hit_def)

        # Skip self-hits
        if hit_id == query_id:
            continue

        is_same_sample = is_self_hit(hit_def)
        taxonomy = extract_taxonomy(hit_def) if not is_same_sample else "Self"

        # Get the highest coverage HSP for this hit
        best_hsp = max(alignment.hsps, key=lambda hsp: (hsp.align_length / blast_record.query_length) * 100)
        
        identity_percentage = (best_hsp.identities / best_hsp.align_length) * 100 if best_hsp.align_length > 0 else 0
        coverage_percentage = min((best_hsp.align_length / blast_record.query_length) * 100, 100)

        hits_info.append({
            "hit_id": hit_id,
            "identity": identity_percentage,
            "coverage": coverage_percentage,
            "is_same_sample": is_same_sample,
            "taxonomy": taxonomy,
            "force_chimeric": False,
            "multiple_hsps": False
        })

    return hits_info

def classify_sequence(hits_info):
    """
    Classify sequence based on hit information and updated criteria.
    Returns a tuple: (classification, reason)
    """
    if not hits_info:
        return "non_chimeric", "No significant non-self hits"

    # Check if any hit is forced to be chimeric (multiple HSPs in first non-self hit)
    if any(hit.get("force_chimeric", False) for hit in hits_info):
        return "chimeric", "Multiple alignments in first non-self hit"

    # Separate database hits and self hits
    database_hits = [hit for hit in hits_info if not hit["is_same_sample"]]
    self_hits = [hit for hit in hits_info if hit["is_same_sample"]]

    if not database_hits and self_hits:
        return "chimeric", "Only self-hits, no database hits"

    if database_hits:
        # Check for high-quality matches first
        high_quality_db_hits = [
            hit for hit in database_hits
            if hit["identity"] >= HIGH_IDENTITY_THRESHOLD
            and hit["coverage"] >= HIGH_COVERAGE_THRESHOLD
        ]
        if high_quality_db_hits:
            return "non_chimeric", "High-quality match against database"

        # Check for database hits with high coverage (≥89%)
        high_coverage_db_hits = [
            hit for hit in database_hits
            if hit["coverage"] >= 89.0
        ]
        
        # If any database hit has coverage ≥89%, classify as non-chimeric
        if high_coverage_db_hits:
            max_coverage_hit = max(high_coverage_db_hits, key=lambda x: x["coverage"])
            return "non_chimeric", f"Database hit with high coverage ({max_coverage_hit['coverage']:.2f}%)"

        # Group remaining hits by taxonomy
        taxonomy_groups = defaultdict(list)
        for hit in database_hits:
            taxonomy_groups[hit["taxonomy"]].append(hit)

        unique_taxa = len(taxonomy_groups)

        if unique_taxa == 1:
            return "borderline", "Multiple database hits pointing to the same taxa, none high-quality"
        else:
            return "chimeric", "Multiple database hits pointing to different taxa, none high-quality"

    return "borderline", "Ambiguous classification"

def write_sequences_to_file(seq_ids, fasta_sequences, output_file):
    """Write sequences to a FASTA file only if seq_ids is not empty."""
    if not seq_ids:
        logging.info(f"No sequences to write for {output_file}. Skipping file creation.")
        return
    
    try:
        with open(output_file, 'w', buffering=8192) as out_file:
            for seq_id in sorted(seq_ids):
                if seq_id in fasta_sequences:  # Added safety check
                    out_file.write(f">{seq_id}\n{fasta_sequences[seq_id]}\n")
        logging.info(f"Wrote {len(seq_ids)} sequences to {output_file}")
    except Exception as e:
        logging.error(f"Error writing to {output_file}: {e}")

def write_sequence_details(details, temp_2_dir, fasta_file):
    """
    Write sequence details to a CSV file.
    
    Args:
        details (list): List of sequence detail rows
        temp_2_dir (str): Directory path for output
        fasta_file (str): Original FASTA file name used to generate output filename
    """
    if not details:
        logging.debug("No sequence details to write. Skipping.")
        return
    
    try:
        # Generate output filename based on input FASTA file
        base_filename = os.path.basename(fasta_file).replace(".chimeras.fasta", "")
        csv_file_path = os.path.join(temp_2_dir, f"{base_filename}_sequence_details.csv")
        
        # Check if file exists to determine if we need to write headers
        file_exists = os.path.isfile(csv_file_path)
        
        # Open file in append mode with buffering for better performance
        with open(csv_file_path, 'a', newline='', buffering=8192) as csvfile:
            csvwriter = csv.writer(csvfile)
            
            # Write headers if file is new
            if not file_exists:
                csvwriter.writerow([
                    "Sequence ID",
                    "Query Coverage (%)",
                    "Identity Percentage (%)",
                    "Classification",
                    "Hit Type",
                    "Hit Origin",
                    "Taxonomy"
                ])
            
            # Write all sequence details
            csvwriter.writerows(details)
            
        logging.debug(f"Wrote {len(details)} sequence details to {csv_file_path}")
        
    except Exception as e:
        logging.error(f"Error writing sequence details to CSV: {e}")

def clean_directory(dir_path):
    """Remove all contents within a directory without deleting the directory itself."""
    if not os.path.exists(dir_path):
        os.makedirs(dir_path, exist_ok=True)
        return
    for filename in os.listdir(dir_path):
        file_path = os.path.join(dir_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            logging.error(f'Failed to delete {file_path}. Reason: {e}')

def parse_blast_results(args):
    """Parse BLAST XML results and classify sequences into different categories."""
    xml_file, temp_dir, temp_2_dir, fasta_file = args
    logging.info(f"Starting processing for {xml_file}")
    
    try:
        fasta_sequences = load_fasta_sequences(fasta_file)
    except FileNotFoundError as e:
        logging.error(e)
        return set(), set(), set(), set()
    except Exception as e:
        logging.error(f"Error loading FASTA file {fasta_file}: {e}")
        return set(), set(), set(), set()

    non_chimeric_sequences = set()
    chimeric_sequences = set()
    borderline_sequences = set()
    multiple_alignment_sequences = set()

    sequence_details = []

    try:
        with open(xml_file) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for blast_record in blast_records:
                query_id = extract_query_id(blast_record.query)
                fasta_seq = fasta_sequences.get(query_id)
                if fasta_seq is None:
                    continue

                hits_info = analyze_blast_hits(blast_record, query_id)
                
                # Check if this is a multiple alignment case
                is_multiple_alignment = any(hit.get("multiple_hsps", False) for hit in hits_info)
                
                if is_multiple_alignment:
                    multiple_alignment_sequences.add(query_id)
                else:
                    # Only classify if not multiple alignment
                    classification, reason = classify_sequence(hits_info)
                    if classification == "non_chimeric":
                        non_chimeric_sequences.add(query_id)
                    elif classification == "chimeric":
                        chimeric_sequences.add(query_id)
                    else:
                        borderline_sequences.add(query_id)
                
                # Add details for all hits
                for i, hit in enumerate(hits_info, 1):
                    classification = "multiple_alignment" if is_multiple_alignment else classification
                    sequence_details.append([
                        query_id, 
                        f"{hit['coverage']:.2f}", 
                        f"{hit['identity']:.2f}", 
                        classification,
                        f"Hit {i}",
                        "Same sample" if hit["is_same_sample"] else "Database",
                        hit["taxonomy"]
                    ])

    except Exception as e:
        logging.error(f"Error processing {xml_file}: {e}")
        return set(), set(), set(), set()

    # Write sequences to respective files
    base_filename = os.path.basename(fasta_file).replace(".chimeras.fasta", "")
    
    # Write multiple alignment sequences to temp_2_dir
    if multiple_alignment_sequences:
        write_sequences_to_file(multiple_alignment_sequences, fasta_sequences, 
                              os.path.join(temp_2_dir, f"{base_filename}_multiple_alignments.fasta"))
    
    # Write other categories
    if non_chimeric_sequences:
        write_sequences_to_file(non_chimeric_sequences, fasta_sequences, 
                              os.path.join(temp_dir, f"{base_filename}_non_chimeric.fasta"))
    if borderline_sequences:
        write_sequences_to_file(borderline_sequences, fasta_sequences, 
                              os.path.join(temp_dir, f"{base_filename}_borderline.fasta"))
    if chimeric_sequences:
        write_sequences_to_file(chimeric_sequences, fasta_sequences, 
                              os.path.join(temp_2_dir, f"{base_filename}_chimeric.fasta"))

    # Write sequence details to temp_2_dir
    if sequence_details:
        write_sequence_details(sequence_details, temp_2_dir, fasta_file)

    logging.info(f"Completed processing for {xml_file}")
    return non_chimeric_sequences, chimeric_sequences, borderline_sequences, multiple_alignment_sequences

def generate_report(all_non_chimeric_sequences, all_chimeric_sequences, 
                   all_borderline_chimeras, all_multiple_alignments, 
                   file_results, output_dir):
    """Generate a detailed report summarizing the results, sorted by file name."""
    try:
        report = {
            "Non-Chimeric Sequences": len(all_non_chimeric_sequences),
            "Chimeric Sequences": len(all_chimeric_sequences),
            "Borderline Sequences": len(all_borderline_chimeras),
            "Multiple Alignment Sequences": len(all_multiple_alignments)
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
        clean_directory(dir_path)
    
    args_list = []
    file_results = defaultdict(lambda: defaultdict(int))
    
    all_non_chimeric_sequences = set()
    all_chimeric_sequences = set()
    all_borderline_chimeras = set()
    all_multiple_alignments = set()

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

    if not args_list:
        logging.error("No valid XML and FASTA file pairs found. Exiting.")
        return

    logging.info(f"Using {NUM_PROCESSES} processes for multiprocessing.")
    
    with multiprocessing.Pool(processes=NUM_PROCESSES) as pool:
        results = pool.map(parse_blast_results, args_list, chunksize=len(args_list)//NUM_PROCESSES)

    # Merge all results
    for i, result in enumerate(results):
        non_chimeric_sequences, chimeric_sequences, borderline_sequences, multiple_alignments = result
        
        # Add sequences to their respective sets
        all_non_chimeric_sequences.update(non_chimeric_sequences)
        all_chimeric_sequences.update(chimeric_sequences)
        all_borderline_chimeras.update(borderline_sequences)
        all_multiple_alignments.update(multiple_alignments)
        
        xml_file = args_list[i][0]
        # Store results for reporting, keeping multiple alignments separate
        file_results[xml_file] = {
            "Non-Chimeric Sequences": len(non_chimeric_sequences),
            "Chimeric Sequences": len(chimeric_sequences),
            "Borderline Sequences": len(borderline_sequences),
            "Multiple Alignment Sequences": len(multiple_alignments)
        }

    # Copy only non-chimeric and borderline sequences to rescued_reads
    for filename in os.listdir(input_dir):
        if filename.endswith(".chimeras.fasta"):
            base_filename = os.path.basename(filename).replace(".chimeras.fasta", "")
            non_chimeric_file = os.path.join(temp_dir, f"{base_filename}_non_chimeric.fasta")
            borderline_file = os.path.join(temp_dir, f"{base_filename}_borderline.fasta")
            
            # Copy only non-chimeric and borderline to rescued_reads
            if os.path.exists(non_chimeric_file):
                shutil.copy(non_chimeric_file, os.path.join(output_dir, f"{base_filename}_non_chimeric.fasta"))
            if os.path.exists(borderline_file):
                shutil.copy(borderline_file, os.path.join(output_dir, f"{base_filename}_borderline.fasta"))

    # Generate report
    generate_report(all_non_chimeric_sequences, all_chimeric_sequences, 
                   all_borderline_chimeras, all_multiple_alignments, 
                   file_results, output_dir)

    log_system_usage()

    # Only remove temp_dir, keep temp_2_dir
    if os.path.exists(temp_dir):
        try:
            shutil.rmtree(temp_dir)
            logging.info(f"Removed {temp_dir} directory.")
        except Exception as e:
            logging.error(f"Failed to remove {temp_dir}. Reason: {e}")

    end_time = time.time()
    elapsed_time = end_time - start_time
    logging.info(f"Total time taken: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    # Define directories
    input_dir = "./input"          # Directory containing FASTA files
    directory = "."                # Directory containing XML files
    temp_dir = "./temp"
    temp_2_dir = "./temp_2"
    output_dir = "./rescued_reads"

    # Ensure output directories exist
    for dir_path in [input_dir, directory]:
        if not os.path.exists(dir_path):
            logging.error(f"Directory does not exist: {dir_path}")
            exit(1)

    try:
        process_all_xml_files(directory, temp_dir, temp_2_dir, output_dir, input_dir)

        # Remove only the temp directory after processing
        if os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
                logging.info(f"Removed {temp_dir} directory.")
            except Exception as e:
                logging.error(f"Failed to remove {temp_dir}. Reason: {e}")

        print("Processing complete. Check the rescued_reads folder for final results.")
    except Exception as e:
        logging.error(f"An error occurred during processing: {e}")
