import pandas as pd

# Define the path to the input FASTA file
input_path = "/path/to/your/dereplicated_metabar_v2.fasta"

# Dictionary to hold sequences and their frequencies
sequences = {}

with open(input_path, "r") as fasta_file:
    seq = ''
    for line in fasta_file:
        line = line.strip()
        if line.startswith(">"):  # Header line
            if seq:  # If there's a sequence from the previous header, save it
                sequences[seq] = freq
                seq = ''
            # Extract frequency from the header
            size_marker = ";size="
            if size_marker in line:
                size_start = line.index(size_marker) + len(size_marker)
                size_end = line[size_start:].find(';') if ';' in line[size_start:] else None
                freq = int(line[size_start:] if size_end is None else line[size_start:size_start+size_end])
        else:
            seq += line
    # Handle the last sequence
    if seq:
        sequences[seq] = freq

# Convert the dictionary to a dataframe
sequence_table = pd.DataFrame(sequences.items(), columns=["Sequence", "Frequency"])

# Save the dataframe to a CSV file
output_path = "/path/to/save/sequence_table_for_dada2.csv"
sequence_table.to_csv(output_path, index=False)
