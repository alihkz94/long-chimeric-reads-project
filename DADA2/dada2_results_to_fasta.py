# Import the os module for file operations
import os

# Initialize dictionaries and lists
fasta_dict = {}
seq_ids_to_extract = []

# File paths
fasta_path = 'chimera.fasta'  # Replace with the path to your chimera.fasta file
seq_ids_path = 'chimeric_sequences.txt'  # Replace with the path to your chimeric_sequences.txt file

# Manually parse the chimera.fasta file and populate the fasta_dict
with open(fasta_path, 'r') as f:
    current_header = ""
    current_sequence = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_header:
                fasta_dict[current_header] = current_sequence
            current_header = line[1:].split(" ")[0]  # Take only the first part of the header (identifier)
            current_sequence = ""
        else:
            current_sequence += line
    # Add the last sequence
    if current_header:
        fasta_dict[current_header] = current_sequence

# Manually parse the chimeric_sequences.txt file and populate the seq_ids_to_extract list
with open(seq_ids_path, 'r') as f:
    for line in f:
        line = line.strip()
        seq_ids_to_extract.append(line)

# Match based on sequence content
matching_sequences = {key: value for key, value in fasta_dict.items() if value in seq_ids_to_extract}

# Generate a new FASTA file using the matching_sequences dictionary
output_fasta_content = ""
for seq_id, seq_content in matching_sequences.items():
    output_fasta_content += f">{seq_id}\n{seq_content}\n"

# Save the new FASTA file as filtered_chimera.fasta
output_fasta_path = 'filtered_chimera_pooled.fasta'  # Replace with the desired output file path
with open(output_fasta_path, 'w') as f:
    f.write(output_fasta_content)

# Verify if the new file has been generated successfully
file_exists = os.path.exists(output_fasta_path)
file_size = os.path.getsize(output_fasta_path) if file_exists else 0

print(f"File exists: {file_exists}, File size: {file_size} bytes")
