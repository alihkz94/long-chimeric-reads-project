"""
Improved Multiprocessing Chimera Detection for Long-Read Sequencing

UPDATES:
- Implemented multiprocessing for faster processing of multiple BLAST XML files.
- Added detailed error handling and logging.
- Improved report generation with per-file statistics.
- Optimized memory usage and performance.

Description:
    This script processes BLAST XML results to identify and classify chimeric sequences in long-read sequencing data.
    It categorizes sequences as false positive chimeras, absolute chimeras, uncertain chimeras, or non-chimeric sequences
    based on specified alignment criteria. The script utilizes multiprocessing for efficient handling of large datasets.

Usage:
    python improved_chimera_detection.py

Requirements:
    - Python 3.6+
    - BioPython

Dependencies:
    - os
    - multiprocessing
    - Bio.Blast.NCBIXML (from BioPython)
    - Bio.SeqIO (from BioPython)
    - collections (Counter, defaultdict)

Input:
    - FASTA files in the specified input directory
    - BLAST XML result files in the current directory

Output:
    - Classified sequences in separate FASTA files within the output directory
    - Detailed report file summarizing overall and per-file results

Configuration:
    Modify the following variables at the beginning of the script to customize paths and thresholds:
    - input_dir: Directory containing input FASTA files
    - directory: Directory containing BLAST XML files (current directory by default)
    - output_dir: Directory for storing output files
    - NUM_PROCESSES: Number of processes to use for multiprocessing
    - HIGH_IDENTITY_THRESHOLD: Threshold for high identity percentage
    - HIGH_COVERAGE_THRESHOLD: Threshold for high query coverage
    - SIGNIFICANT_COVERAGE_THRESHOLD: Threshold for significant query coverage
    - SIGNIFICANT_IDENTITY_THRESHOLD: Threshold for significant identity percentage

Author: Ali Hakimzadeh
Version: 0.3.0
Date: 2024-09-11
"""
                    
import os
import multiprocessing
from Bio.Blast import NCBIXML
from Bio import SeqIO
from collections import Counter, defaultdict

# Constants for classification criteria
HIGH_IDENTITY_THRESHOLD = 99.0
HIGH_COVERAGE_THRESHOLD = 99.0
SIGNIFICANT_COVERAGE_THRESHOLD = 80.0
SIGNIFICANT_IDENTITY_THRESHOLD = 80.0

# Set the number of processes for multiprocessing
NUM_PROCESSES = 4  # Adjust this based on your system's capabilities

def load_fasta_sequences(fasta_file):
    """
    Load sequences from a FASTA file into a dictionary.
    
    Args:
        fasta_file (str): Path to the FASTA file.
    
    Returns:
        dict: A dictionary with sequence IDs as keys and sequences as values.
    
    Raises:
        FileNotFoundError: If the FASTA file does not exist.
    """
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
    
    fasta_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_sequences[record.id] = str(record.seq)
    print(f"Loaded {len(fasta_sequences)} sequences from {fasta_file}")
    return fasta_sequences

def extract_query_id(blast_query_def):
    """
    Extract the query ID from the BLAST query definition.
    
    Args:
        blast_query_def (str): BLAST query definition.
    
    Returns:
        str: Extracted query ID.
    """
    return blast_query_def

def parse_blast_results(args):
    """
    Parse BLAST XML results and classify sequences into different categories.
    
    Args:
        args (tuple): A tuple containing the XML file path, output directory path, and FASTA file path.
    
    Returns:
        tuple: A tuple containing sets of false positive chimeras, absolute chimeras, uncertain chimeras, and non-chimeric sequences.
    
    Raises:
        FileNotFoundError: If the XML file or FASTA file does not exist.
        ValueError: If no records are found in the BLAST XML file.
    """
    xml_file, output_dir, fasta_file = args
    
    if not os.path.exists(xml_file):
        raise FileNotFoundError(f"XML file not found: {xml_file}")
    
    try:
        fasta_sequences = load_fasta_sequences(fasta_file)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return set(), set(), set(), set()
    
    false_positive_chimeras = set()
    absolute_chimeras = set()
    uncertain_chimeras = set()
    non_chimeric_sequences = set()

    try:
        with open(xml_file) as result_handle:
            blast_records = list(NCBIXML.parse(result_handle))
            
            if not blast_records:
                raise ValueError(f"No records found in BLAST XML: {xml_file}")

            for blast_record in blast_records:
                query_id = extract_query_id(blast_record.query)
                print(f"Processing Query ID: {query_id}")

                if query_id not in fasta_sequences:
                    print(f"Warning: Query ID {query_id} not found in FASTA file {fasta_file}")
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
                    # Check if the high identity hit is not from the same sample
                    if all(extract_query_id(align[0].hit_def).split('_')[0] != query_id.split('_')[0] for align in high_identity_alignments):
                        false_positive_chimeras.add(query_id)
                        print(f"False positive chimera: {query_id}")
                    else:
                        uncertain_chimeras.add(query_id)
                        print(f"Uncertain chimera (high identity, same sample): {query_id}")
                elif len(significant_alignments) > 1:
                    absolute_chimeras.add(query_id)
                    print(f"Absolute chimera: {query_id}")
                elif significant_alignments:
                    uncertain_chimeras.add(query_id)
                    print(f"Uncertain chimera: {query_id}")
                else:
                    non_chimeric_sequences.add(query_id)
                    print(f"Non-chimeric sequence: {query_id}")

            # Write results to files
            base_filename = os.path.basename(fasta_file).replace(".chimeras.fasta", "")
            write_sequences_to_file(false_positive_chimeras, fasta_sequences, os.path.join(output_dir, f"{base_filename}_false_positive_chimeras.fasta"))
            write_sequences_to_file(absolute_chimeras, fasta_sequences, os.path.join(output_dir, f"{base_filename}_absolute_chimeras.fasta"))
            write_sequences_to_file(uncertain_chimeras, fasta_sequences, os.path.join(output_dir, f"{base_filename}_uncertain_chimeras.fasta"))
            write_sequences_to_file(non_chimeric_sequences, fasta_sequences, os.path.join(output_dir, f"{base_filename}_non_chimeric.fasta"))

    except Exception as e:
        print(f"Error processing BLAST XML file {xml_file}: {e}")

    return false_positive_chimeras, absolute_chimeras, uncertain_chimeras, non_chimeric_sequences

def write_sequences_to_file(seq_ids, fasta_sequences, output_file):
    """
    Write sequences to a FASTA file.
    
    Args:
        seq_ids (set): Set of sequence IDs to write.
        fasta_sequences (dict): Dictionary of all sequences.
        output_file (str): Path to the output file.
    """
    with open(output_file, 'w') as out_file:
        for seq_id in sorted(seq_ids):
            out_file.write(f">{seq_id}\n{fasta_sequences[seq_id]}\n")

def process_all_xml_files(directory, output_dir, input_dir):
    """
    Process all BLAST XML files in a directory using multiprocessing.
    
    Args:
        directory (str): Directory containing BLAST XML files.
        output_dir (str): Directory for output files.
        input_dir (str): Directory containing input FASTA files.
    
    Returns:
        tuple: A tuple containing sets of all false positive chimeras, absolute chimeras, uncertain chimeras, and non-chimeric sequences.
    """
    all_false_positive_chimeras = set()
    all_absolute_chimeras = set()
    all_uncertain_chimeras = set()
    all_non_chimeric_sequences = set()
    file_results = defaultdict(lambda: defaultdict(int))

    args_list = []
    for filename in os.listdir(directory):
        if filename.endswith("_blast_results.xml"):
            xml_file = os.path.join(directory, filename)
            fasta_file = os.path.join(input_dir, filename.replace("_blast_results.xml", ".chimeras.fasta"))
            args_list.append((xml_file, output_dir, fasta_file))

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

    return all_false_positive_chimeras, all_absolute_chimeras, all_uncertain_chimeras, all_non_chimeric_sequences, file_results

def generate_report(all_false_positive_chimeras, all_absolute_chimeras, all_uncertain_chimeras, all_non_chimeric_sequences, file_results, output_dir):
    """
    Generate a detailed report summarizing the results.
    
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
        for file, results in file_results.items():
            report_file.write(f"\nFile: {file}\n")
            file_total = sum(results.values())
            for category, count in results.items():
                percentage = (count / file_total) * 100 if file_total > 0 else 0
                report_file.write(f"  {category}: {count} ({percentage:.2f}%)\n")
            report_file.write(f"  Total: {file_total}\n")

# Main script execution
if __name__ == "__main__":
    input_dir = "./input"
    directory = "."
    output_dir = "./output"
    os.makedirs(output_dir, exist_ok=True)

    try:
        all_false_positive_chimeras, all_absolute_chimeras, all_uncertain_chimeras, all_non_chimeric_sequences, file_results = process_all_xml_files(directory, output_dir, input_dir)
        generate_report(all_false_positive_chimeras, all_absolute_chimeras, all_uncertain_chimeras, all_non_chimeric_sequences, file_results, output_dir)
        print("Processing complete. Check the output directory for results and report.")
    except Exception as e:
        print(f"An error occurred during processing: {e}")
