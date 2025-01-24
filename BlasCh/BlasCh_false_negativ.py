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

# Classification thresholds
HIGH_IDENTITY_THRESHOLD = 99.0
HIGH_COVERAGE_THRESHOLD = 99.0
BORDERLINE_COVERAGE_THRESHOLD = 89.0
MULTIPLE_ALIGNMENT_COVERAGE_LIMIT = 85.0

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

def analyze_blast_hits(blast_record, query_id):
    """
    Analyze BLAST hits for a given blast record by processing all HSPs.
    Return:
      hits_info: list of dicts with info about each hit
      best_coverage_first_hit: coverage of the best HSP in the first non-self alignment (for multiple HSP logic)
    """
    hits_info = []
    if not blast_record.alignments:
        return hits_info, 0.0  # No hits at all
    
    # Identify the first alignment and see if it is a self-hit
    first_alignment = blast_record.alignments[0]
    first_hit_id = extract_query_id(first_alignment.hit_def)
    
    # Decide which alignment is the first non-self
    hit_to_check = first_alignment
    if first_hit_id == query_id and len(blast_record.alignments) > 1:
        hit_to_check = blast_record.alignments[1]
        hit_id = extract_query_id(hit_to_check.hit_def)
    else:
        hit_id = first_hit_id

    multiple_hsps_detected = False
    best_coverage_first_hit = 0.0

    # If the first non-self has multiple HSPs, compute the best coverage among them
    if hit_id != query_id and len(hit_to_check.hsps) > 1:
        multiple_hsps_detected = True
        best_hsp = max(hit_to_check.hsps, key=lambda hsp: (hsp.align_length / blast_record.query_length) * 100)
        best_coverage_first_hit = min((best_hsp.align_length / blast_record.query_length) * 100, 100)

    # Gather hit info for each alignment
    for alignment in blast_record.alignments:
        if not alignment.hsps:
            continue
        
        curr_hit_def = alignment.hit_def
        curr_hit_id = extract_query_id(curr_hit_def)
        
        # Skip self-hits if curr_hit_id == query_id
        if curr_hit_id == query_id:
            continue
        
        same_sample = is_self_hit(curr_hit_def)
        
        # Use the entire <Hit_def> string in the taxonomy field
        # This way, you get exactly what's in <Hit_def> in the output CSV
        taxonomy_str = curr_hit_def

        best_hsp = max(alignment.hsps, key=lambda hsp: (hsp.align_length / blast_record.query_length) * 100)
        identity_percentage = (best_hsp.identities / best_hsp.align_length) * 100 if best_hsp.align_length > 0 else 0
        coverage_percentage = min((best_hsp.align_length / blast_record.query_length) * 100, 100)
        
        hits_info.append({
            "hit_id": curr_hit_id,
            "identity": identity_percentage,
            "coverage": coverage_percentage,
            "is_same_sample": same_sample,
            "taxonomy": taxonomy_str,
            "multiple_hsps": multiple_hsps_detected
        })
    
    return hits_info, best_coverage_first_hit

def classify_based_on_coverage(hits_info):
    """
    Classify the set of hits into 'nonchimeric' or 'borderline' based on coverage rules.
    - If no database hits => nonchimeric
    - If any database hit has coverage >= 89 => nonchimeric
    - Else => borderline
    """
    database_hits = [hit for hit in hits_info if not hit["is_same_sample"]]
    if not database_hits:
        return "nonchimeric"
    
    best_db_hit = max(database_hits, key=lambda x: x["coverage"])
    if best_db_hit["coverage"] >= BORDERLINE_COVERAGE_THRESHOLD:
        return "nonchimeric"
    else:
        return "borderline"

def parse_blast_results(args):
    """
    Parse BLAST XML results and classify sequences into:
      1) multiple_alignment
      2) nonchimeric
      3) borderline
    Always returns 3 sets: (nonchimeric, borderline, multiple_alignment).
    """
    xml_file, temp_dir, temp_2_dir, fasta_file = args
    logging.info(f"Starting processing for {xml_file}")

    # Load FASTA sequences
    try:
        fasta_sequences = load_fasta_sequences(fasta_file)
    except FileNotFoundError as e:
        logging.error(e)
        return set(), set(), set()
    except Exception as e:
        logging.error(f"Error loading FASTA file {fasta_file}: {e}")
        return set(), set(), set()

    nonchimeric_sequences = set()
    borderline_sequences = set()
    multiple_alignment_sequences = set()

    sequence_details = []

    try:
        with open(xml_file) as result_handle:
            blast_records = NCBIXML.parse(result_handle)

            for blast_record in blast_records:
                query_id = extract_query_id(blast_record.query)
                if query_id not in fasta_sequences:
                    # Skip if query not in FASTA
                    continue
                
                hits_info, best_cov_first = analyze_blast_hits(blast_record, query_id)

                if not hits_info:
                    # No hits => nonchimeric
                    classification = "nonchimeric"
                else:
                    # Check if multiple HSPs was detected for the first non-self alignment
                    if hits_info[0]["multiple_hsps"]:
                        # If best coverage of that first non-self alignment > 85 => nonchimeric
                        if best_cov_first > MULTIPLE_ALIGNMENT_COVERAGE_LIMIT:
                            classification = "nonchimeric"
                        else:
                            classification = "multiple_alignment"
                    else:
                        # Otherwise, classify by coverage
                        classification = classify_based_on_coverage(hits_info)

                # Assign to sets
                if classification == "multiple_alignment":
                    multiple_alignment_sequences.add(query_id)
                elif classification == "borderline":
                    borderline_sequences.add(query_id)
                else:
                    nonchimeric_sequences.add(query_id)

                # Add details for CSV
                for i, hit in enumerate(hits_info, start=1):
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
        return set(), set(), set()

    # Write sequences to files
    base_filename = os.path.basename(fasta_file).replace(".chimeras.fasta", "")

    if multiple_alignment_sequences:
        write_sequences_to_file(
            multiple_alignment_sequences,
            fasta_sequences,
            os.path.join(temp_2_dir, f"{base_filename}_multiple_alignments.fasta")
        )
    if nonchimeric_sequences:
        write_sequences_to_file(
            nonchimeric_sequences,
            fasta_sequences,
            os.path.join(temp_dir, f"{base_filename}_nonchimeric.fasta")
        )
    if borderline_sequences:
        write_sequences_to_file(
            borderline_sequences,
            fasta_sequences,
            os.path.join(temp_dir, f"{base_filename}_borderline.fasta")
        )

    # Write CSV details
    if sequence_details:
        write_sequence_details(sequence_details, temp_2_dir, fasta_file)

    logging.info(f"Completed processing for {xml_file}")
    return nonchimeric_sequences, borderline_sequences, multiple_alignment_sequences

def write_sequences_to_file(seq_ids, fasta_sequences, output_file):
    """Write sequences to a FASTA file only if seq_ids is not empty."""
    if not seq_ids:
        logging.info(f"No sequences to write for {output_file}. Skipping.")
        return
    try:
        with open(output_file, 'w', buffering=8192) as out_file:
            for seq_id in sorted(seq_ids):
                if seq_id in fasta_sequences:
                    out_file.write(f">{seq_id}\n{fasta_sequences[seq_id]}\n")
        logging.info(f"Wrote {len(seq_ids)} sequences to {output_file}")
    except Exception as e:
        logging.error(f"Error writing to {output_file}: {e}")

def write_sequence_details(details, temp_2_dir, fasta_file):
    """
    Write sequence details to a CSV file.
    """
    if not details:
        logging.debug("No sequence details to write. Skipping.")
        return
    try:
        base_filename = os.path.basename(fasta_file).replace(".chimeras.fasta", "")
        csv_file_path = os.path.join(temp_2_dir, f"{base_filename}_sequence_details.csv")
        
        file_exists = os.path.isfile(csv_file_path)
        
        with open(csv_file_path, 'a', newline='', buffering=8192) as csvfile:
            csvwriter = csv.writer(csvfile)
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

def generate_report(all_nonchimeric, all_borderline, all_multiple, file_results, output_dir):
    """
    Generate a detailed text report summarizing the results.
    Only three categories: Non-Chimeric, Borderline, Multiple Alignment.
    """
    try:
        report = {
            "Non-Chimeric Sequences": len(all_nonchimeric),
            "Borderline Sequences": len(all_borderline),
            "Multiple Alignment Sequences": len(all_multiple)
        }

        total_sequences = sum(report.values())
        report_path = os.path.join(output_dir, "chimera_detection_report.txt")

        with open(report_path, 'w') as report_file:
            report_file.write("Classification Report\n")
            report_file.write("=====================\n\n")
            report_file.write("Overall Summary:\n")
            report_file.write("----------------\n")
            for category, count in report.items():
                percentage = (count / total_sequences) * 100 if total_sequences else 0
                report_file.write(f"{category}: {count} ({percentage:.2f}%)\n")
            report_file.write(f"\nTotal Sequences Processed: {total_sequences}\n\n")
            
            report_file.write("Detailed Results by File:\n")
            report_file.write("-------------------------\n")
            
            sorted_file_results = dict(sorted(file_results.items()))
            for file_path, results in sorted_file_results.items():
                report_file.write(f"\nFile: {os.path.basename(file_path)}\n")
                file_total = sum(results.values())
                for category, count in results.items():
                    percentage = (count / file_total) * 100 if file_total else 0
                    report_file.write(f"  {category}: {count} ({percentage:.2f}%)\n")
                report_file.write(f"  Total: {file_total}\n")
        logging.info(f"Report generated: {report_path}")
    except Exception as e:
        logging.error(f"Error generating report: {e}")

def process_all_xml_files(directory, temp_dir, temp_2_dir, output_dir, input_dir):
    """Process all BLAST XML files and generate final outputs into 3 categories."""
    for dir_path in [temp_dir, temp_2_dir, output_dir]:
        clean_directory(dir_path)
    
    args_list = []
    file_results = defaultdict(lambda: defaultdict(int))
    all_nonchimeric = set()
    all_borderline = set()
    all_multiple = set()

    start_time = time.time()
    log_system_usage()

    # Gather all relevant XML + FASTA pairs
    for filename in os.listdir(directory):
        if filename.endswith("_blast_results.xml"):
            xml_file = os.path.join(directory, filename)
            # Match with same base name .chimeras.fasta
            fasta_file = os.path.join(input_dir, filename.replace("_blast_results.xml", ".chimeras.fasta"))
            if os.path.exists(fasta_file):
                args_list.append((xml_file, temp_dir, temp_2_dir, fasta_file))
            else:
                logging.warning(f"FASTA file for {xml_file} not found. Skipping.")

    if not args_list:
        logging.error("No valid XML and FASTA file pairs found. Exiting.")
        return

    # Safely calculate chunk size
    if len(args_list) < NUM_PROCESSES:
        safe_chunk_size = 1
    else:
        safe_chunk_size = len(args_list) // NUM_PROCESSES
        if safe_chunk_size == 0:
            safe_chunk_size = 1

    logging.info(f"Using {NUM_PROCESSES} processes for multiprocessing with chunk size {safe_chunk_size}.")

    with multiprocessing.Pool(processes=NUM_PROCESSES) as pool:
        results = pool.map(parse_blast_results, args_list, chunksize=safe_chunk_size)

    # Combine results
    for i, result in enumerate(results):
        if result is None:
            # If parse_blast_results returned None due to an unhandled error, skip
            logging.error("A multiprocessing worker returned None. Check parse_blast_results for uncaught exceptions.")
            continue

        try:
            nonchimeric_sequences, borderline_sequences, multiple_alignment_sequences = result
        except Exception as e:
            logging.error(f"Error unpacking result #{i}: {e}")
            continue

        all_nonchimeric.update(nonchimeric_sequences)
        all_borderline.update(borderline_sequences)
        all_multiple.update(multiple_alignment_sequences)
        
        xml_file = args_list[i][0]
        file_results[xml_file] = {
            "Non-Chimeric Sequences": len(nonchimeric_sequences),
            "Borderline Sequences": len(borderline_sequences),
            "Multiple Alignment Sequences": len(multiple_alignment_sequences)
        }

    # Copy nonchimeric & borderline sequences into final output
    for filename in os.listdir(input_dir):
        if filename.endswith(".chimeras.fasta"):
            base_filename = os.path.basename(filename).replace(".chimeras.fasta", "")
            nonchimeric_file = os.path.join(temp_dir, f"{base_filename}_nonchimeric.fasta")
            borderline_file = os.path.join(temp_dir, f"{base_filename}_borderline.fasta")
            
            if os.path.exists(nonchimeric_file):
                shutil.copy(nonchimeric_file, os.path.join(output_dir, f"{base_filename}_nonchimeric.fasta"))
            if os.path.exists(borderline_file):
                shutil.copy(borderline_file, os.path.join(output_dir, f"{base_filename}_borderline.fasta"))

    # Generate the final report
    generate_report(all_nonchimeric, all_borderline, all_multiple, file_results, output_dir)

    log_system_usage()

    # Clean up temp_dir, keep temp_2_dir
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
    # Example usage
    input_dir = "./input"          # Directory containing FASTA files
    directory = "."                # Directory containing XML files
    temp_dir = "./temp"
    temp_2_dir = "./temp_2"
    output_dir = "./rescued_reads"

    for dir_path in [input_dir, directory]:
        if not os.path.exists(dir_path):
            logging.error(f"Directory does not exist: {dir_path}")
            exit(1)

    try:
        process_all_xml_files(directory, temp_dir, temp_2_dir, output_dir, input_dir)
        # Attempt final removal of temp_dir if still present
        if os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
                logging.info(f"Removed {temp_dir} directory.")
            except Exception as e:
                logging.error(f"Failed to remove {temp_dir}. Reason: {e}")

        print("Processing complete. Check the rescued_reads folder for final results.")
    except Exception as e:
        logging.error(f"An error occurred during processing: {e}")
