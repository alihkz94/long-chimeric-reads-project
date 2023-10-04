# Import required modules
from collections import Counter
import pandas as pd

# Custom FASTA parser function
def parse_fasta(file_path):
    sequences = []
    with open(file_path, "r") as f:
        sequence = ""
        for line in f:
            if line.startswith(">"):
                if sequence:
                    sequences.append(sequence)
                sequence = ""
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence)
    return sequences

# Read sequences using custom FASTA parser
fasta_file_path = "chimeric_input.fasta"
sequences = parse_fasta(fasta_file_path)

# Count the frequency of each unique sequence
sequence_counter = Counter(sequences)

# Create a DataFrame from the sequence counter
seqtab_df = pd.DataFrame.from_dict(sequence_counter, orient='index', columns=['Frequency']).reset_index()
seqtab_df.columns = ['Sequence', 'Frequency']

# Save the DataFrame to a CSV file
csv_file_path = "chimeric_sequence_table.csv"
seqtab_df.to_csv(csv_file_path, index=False)
