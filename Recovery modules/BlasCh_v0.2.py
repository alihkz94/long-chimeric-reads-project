import os
import multiprocessing
from Bio.Blast import NCBIXML
from Bio import SeqIO
from collections import Counter

# Function to load sequences from a FASTA file into a dictionary
def load_fasta_sequences(fasta_file):
    fasta_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_sequences[record.id] = str(record.seq)
    print(f"Loaded {len(fasta_sequences)} sequences from {fasta_file}")
    return fasta_sequences

# Function to extract the query ID from the BLAST query definition
def extract_query_id(blast_query_def):
    return blast_query_def

# Function to parse BLAST XML results and classify sequences into different categories
def parse_blast_results(args):
    """
    Parses BLAST XML results and classifies sequences into different categories based on alignment criteria.
    Args:
        args (tuple): A tuple containing the XML file path, output directory path, and FASTA file path.
    Returns:
        tuple: A tuple containing sets of false positive chimeras, absolute chimeras, uncertain chimeras, and non-chimeric sequences.
    Raises:
        Exception: If there is an error processing the BLAST XML file.
    """
    xml_file, output_dir, fasta_file = args
    fasta_sequences = load_fasta_sequences(fasta_file)
    
    false_positive_chimeras = set()
    absolute_chimeras = set()
    uncertain_chimeras = set()
    non_chimeric_sequences = set()

    try:
        with open(xml_file) as result_handle:
            blast_records = list(NCBIXML.parse(result_handle))
            
            if not blast_records:
                print(f"No records found in BLAST XML: {xml_file}")
                return false_positive_chimeras, absolute_chimeras, uncertain_chimeras, non_chimeric_sequences

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
                            
                            if identity_percentage >= 99.0 and query_coverage >= 99.0:
                                high_identity_alignments.append((alignment, hsp))
                            elif query_coverage >= 80.0 and identity_percentage >= 80.0:
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

# Function to write sequences to a FASTA file
def write_sequences_to_file(seq_ids, fasta_sequences, output_file):
    with open(output_file, 'w') as out_file:
        for seq_id in sorted(seq_ids):
            out_file.write(f">{seq_id}\n{fasta_sequences[seq_id]}\n")

# Function to process all BLAST XML files in a directory using multiprocessing
def process_all_xml_files(directory, output_dir, input_dir):
    all_false_positive_chimeras = set()
    all_absolute_chimeras = set()
    all_uncertain_chimeras = set()
    all_non_chimeric_sequences = set()

    args_list = []
    for filename in os.listdir(directory):
        if filename.endswith("_blast_results.xml"):
            xml_file = os.path.join(directory, filename)
            fasta_file = os.path.join(input_dir, filename.replace("_blast_results.xml", ".chimeras.fasta"))
            args_list.append((xml_file, output_dir, fasta_file))

    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.map(parse_blast_results, args_list)

    for result in results:
        false_positive_chimeras, absolute_chimeras, uncertain_chimeras, non_chimeric_sequences = result
        all_false_positive_chimeras.update(false_positive_chimeras)
        all_absolute_chimeras.update(absolute_chimeras)
        all_uncertain_chimeras.update(uncertain_chimeras)
        all_non_chimeric_sequences.update(non_chimeric_sequences)

    return all_false_positive_chimeras, all_absolute_chimeras, all_uncertain_chimeras, all_non_chimeric_sequences

# Function to generate a report summarizing the results
def generate_report(all_false_positive_chimeras, all_absolute_chimeras, all_uncertain_chimeras, all_non_chimeric_sequences, output_dir):
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
        for category, count in report.items():
            percentage = (count / total_sequences) * 100
            report_file.write(f"{category}: {count} ({percentage:.2f}%)\n")
        report_file.write(f"\nTotal Sequences Processed: {total_sequences}\n")

# Main script execution
input_dir = "./input"
directory = "."
output_dir = "./output"
os.makedirs(output_dir, exist_ok=True)

all_false_positive_chimeras, all_absolute_chimeras, all_uncertain_chimeras, all_non_chimeric_sequences = process_all_xml_files(directory, output_dir, input_dir)

generate_report(all_false_positive_chimeras, all_absolute_chimeras, all_uncertain_chimeras, all_non_chimeric_sequences, output_dir)

print("Processing complete. Check the output directory for results and report.")
