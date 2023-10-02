import csv

# Load and parse files
def load_fasta(file_path):
    with open(file_path, 'r') as file:
        sequences = {}
        for line in file:
            if line.startswith('>'):
                header = line.strip()[1:]
                sequences[header] = next(file).strip()
        return sequences

def load_csv(file_path):
    with open(file_path, 'r') as file:
        return list(csv.DictReader(file))

def load_tsv(file_path):
    with open(file_path, 'r') as file:
        return list(csv.DictReader(file, delimiter='\t'))

# Load data
original_sequences = load_fasta('/home/alihkz/dada2_rate/chimeric_input.fasta')
filtered_chimeras = load_fasta('/home/alihkz/dada2_rate/filtered_chimera.fasta')
chimera_info = load_tsv('/home/alihkz/dada2_rate/chimera_info.tsv')

# Identify sequence ID key
seq_id_key = list(chimera_info[0].keys())[0]

# Extract chimera IDs
chimera_ids = set(entry[seq_id_key] for entry in chimera_info)

# Initialize counters
TP, FP, TN, FN = 0, 0, 0, 0

# Debugging: Verify data integrity
print(f"Total original sequences: {len(original_sequences)}")
print(f"Total filtered chimeras: {len(filtered_chimeras)}")
print(f"Total chimera info records: {len(chimera_info)}")

# Compute metrics
for header, sequence in original_sequences.items():
    is_chimera = header in chimera_ids
    is_filtered = header in filtered_chimeras

    if is_chimera and is_filtered:
        TP += 1
    elif is_chimera and not is_filtered:
        FN += 1
    elif not is_chimera and is_filtered:
        FP += 1
    elif not is_chimera and not is_filtered:
        TN += 1

# Debugging: Verify calculations
print(f"Debugging Info: TP+FP+TN+FN = {TP+FP+TN+FN}")

# Calculate Sensitivity and Specificity
sensitivity = TP / (TP + FN) if TP + FN > 0 else 'Undefined'
specificity = TN / (TN + FP) if TN + FP > 0 else 'Undefined'

# Print Results
result_dict = {
    'True Positives': TP,
    'False Positives': FP,
    'True Negatives': TN,
    'False Negatives': FN,
    'Sensitivity': sensitivity,
    'Specificity': specificity
}
print("Calculated Metrics:", result_dict)
