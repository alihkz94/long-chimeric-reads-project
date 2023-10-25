import os
import glob
from Bio import SeqIO
import concurrent.futures

# Function to handle duplicate keys (from your previous script)
def handle_duplicate_keys(record_iter):
    counts = {}
    for record in record_iter:
        key = record.id
        if key in counts:
            counts[key] += 1
            key = f"{key}_{counts[key]}"
        else:
            counts[key] = 0
        record.id = key
        yield record

# Create a directory to store the new FASTA files
if not os.path.exists('fasta_files'):
    os.makedirs('fasta_files')

# Read metabar FASTA file once, handling duplicate keys
metabar_dict = SeqIO.to_dict(handle_duplicate_keys(SeqIO.parse("metabar_uchime_input.fasta", "fasta")))

# Function to process each removed chimeric sequences file
def process_removed_file(filepath):
    # Initialize a list to store new FASTA records
    new_fasta_records = []
    
    # Read the sequences in the current removed_chimeric_sequences file
    with open(filepath, 'r') as f:
        removed_seqs = set(f.read().splitlines())
    
    # Map each removed sequence to its original header and sequence in the metabar FASTA file
    for seq_id, sequence in metabar_dict.items():
        if sequence.seq in removed_seqs:
            new_fasta_records.append((seq_id, str(sequence.seq)))
    
    # Write the new FASTA records to a file in the 'fasta_files' directory
    output_filepath = os.path.join('fasta_files', os.path.basename(filepath).replace('.txt', '.fasta'))
    
    with open(output_filepath, 'w') as f:
        for seq_id, seq in new_fasta_records:
            f.write(f">{seq_id}\n{seq}\n")

# Parallel processing of removed chimeric sequences files
with concurrent.futures.ThreadPoolExecutor() as executor:
    executor.map(process_removed_file, glob.glob("removed_chimeric_sequences_*.txt"))



####################ONLY CHIMERIC FASTA HEADERS ONES#########################
import os
import glob
from Bio import SeqIO
import concurrent.futures

# Function to handle duplicate keys (from your previous script)
def handle_duplicate_keys(record_iter):
    counts = {}
    for record in record_iter:
        key = record.id
        if key in counts:
            counts[key] += 1
            key = f"{key}_{counts[key]}"
        else:
            counts[key] = 0
        record.id = key
        yield record

# Create a directory to store the new FASTA files
if not os.path.exists('fasta_files'):
    os.makedirs('fasta_files')

# Read metabar FASTA file once, handling duplicate keys
metabar_dict = SeqIO.to_dict(handle_duplicate_keys(SeqIO.parse("metabar_uchime_input.fasta", "fasta")))

# Function to process each removed chimeric sequences file
def process_removed_file(filepath):
    # Initialize a list to store new FASTA records
    new_fasta_records = []
    
    # Read the sequences in the current removed_chimeric_sequences file
    with open(filepath, 'r') as f:
        removed_seqs = set(f.read().splitlines())
    
    # Map each removed sequence to its original header and sequence in the metabar FASTA file
    # and filter only those that start with "chimera_"
    for seq_id, sequence in metabar_dict.items():
        if sequence.seq in removed_seqs and seq_id.startswith("chimera_"):
            new_fasta_records.append((seq_id, str(sequence.seq)))
    
    # Write the new FASTA records to a file in the 'fasta_files' directory
    output_filepath = os.path.join('fasta_files', os.path.basename(filepath).replace('.txt', '_chimeric.fasta'))
    
    with open(output_filepath, 'w') as f:
        for seq_id, seq in new_fasta_records:
            f.write(f">{seq_id}\n{seq}\n")

# Parallel processing of removed chimeric sequences files
with concurrent.futures.ThreadPoolExecutor() as executor:
    executor.map(process_removed_file, glob.glob("removed_chimeric_sequences_*.txt"))
