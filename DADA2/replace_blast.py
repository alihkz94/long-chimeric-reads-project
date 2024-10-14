import os
import multiprocessing
import logging
import hashlib
import time
from Bio import SeqIO
from functools import partial
import xml.etree.ElementTree as ET
import itertools

# Configure logging to display process name and messages
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(processName)s - %(levelname)s - %(message)s'
)

def compute_md5(sequence):
    """Compute MD5 hash for a given sequence."""
    return hashlib.md5(sequence.encode('utf-8')).hexdigest()

def create_id_mapping(original_dir, chimeras_dir):
    """
    Create a mapping from chimeric IDs to original IDs based on sequence matching.
    """
    id_mapping = {}
    chimera_files = [f for f in os.listdir(chimeras_dir) if f.endswith(".fasta")]

    for filename in chimera_files:
        original_file = os.path.join(original_dir, filename)
        chimera_file = os.path.join(chimeras_dir, filename)

        if not os.path.exists(original_file):
            logging.warning(f"Original file not found for {filename}. Skipping.")
            continue

        try:
            # Load original sequences and create a hash-to-ID mapping
            original_seq_dict = {}
            with open(original_file, 'r') as f:
                for record in SeqIO.parse(f, "fasta"):
                    seq_str = str(record.seq).upper()
                    seq_hash = compute_md5(seq_str)
                    original_seq_dict[seq_hash] = record.description

            # Load chimera sequences
            with open(chimera_file, 'r') as f:
                for record in SeqIO.parse(f, "fasta"):
                    seq_str = str(record.seq).upper()
                    seq_hash = compute_md5(seq_str)
                    chimera_id = record.description

                    if seq_hash in original_seq_dict:
                        id_mapping[chimera_id] = original_seq_dict[seq_hash]
                    else:
                        logging.warning(f"Sequence in {chimera_file} not found in original. ID mapping not created for {chimera_id}.")
        except Exception as e:
            logging.error(f"Error processing {filename}: {e}")

    return id_mapping

def modify_blast_xml_chunk(blast_dir, output_dir, id_mapping, filename, chunk_size=1000):
    """
    Modify the BLAST XML file by replacing chimeric IDs with original IDs in <Iteration_query-def> elements.
    Process the file in chunks to reduce memory usage.
    """
    blast_file = os.path.join(blast_dir, filename)
    output_file = os.path.join(output_dir, filename)

    try:
        # Initialize the output file
        with open(output_file, 'w', encoding='utf-8') as out:
            out.write('<?xml version="1.0"?>\n<BlastOutput>\n')

        # Process the file in chunks
        context = ET.iterparse(blast_file, events=('start', 'end'))
        root = None
        for event, elem in context:
            if event == 'start' and root is None:
                root = elem
            elif event == 'end' and elem.tag == 'Iteration':
                # Process each Iteration element
                query_def = elem.find('Iteration_query-def')
                if query_def is not None:
                    original_query_def = id_mapping.get(query_def.text, query_def.text)
                    query_def.text = original_query_def

                # Write the processed Iteration to the output file
                with open(output_file, 'a', encoding='utf-8') as out:
                    out.write(ET.tostring(elem, encoding='unicode'))

                # Clear the element to free memory
                root.clear()

        # Write the closing tag
        with open(output_file, 'a', encoding='utf-8') as out:
            out.write('</BlastOutput>')

        logging.info(f"Modified BLAST XML file saved: {output_file}")

    except Exception as e:
        logging.error(f"Error processing {filename}: {e}")

def main():
    # Define directories
    original_dir = "./original"
    chimeras_dir = "./nonchimeras"
    blast_dir = "./blast"
    output_dir = "./Blast_modified"

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Create ID mapping from chimeric IDs to original IDs
    logging.info("Creating ID mapping from chimeric IDs to original IDs...")
    id_mapping = create_id_mapping(original_dir, chimeras_dir)
    logging.info(f"ID mapping created with {len(id_mapping)} entries.")

    # Get list of BLAST XML files to process
    blast_files = [f for f in os.listdir(blast_dir) if f.endswith("_blast_results.xml")]

    if not blast_files:
        logging.error("No BLAST XML files found in the blast directory.")
        return

    logging.info(f"Found {len(blast_files)} BLAST XML files to process.")

    # Determine number of processes to use
    num_processes = multiprocessing.cpu_count()
    logging.info(f"Using {num_processes} processes for multiprocessing.")

    # Partial function with fixed arguments
    process_func = partial(
        modify_blast_xml_chunk,
        blast_dir,
        output_dir,
        id_mapping
    )

    start_time = time.time()

    # Use multiprocessing Pool to process files in parallel
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(process_func, blast_files)

    end_time = time.time()
    elapsed_time = end_time - start_time
    logging.info(f"All BLAST XML files processed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()
