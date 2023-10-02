import csv

# Function to parse CSV or TSV files
def parse_delimited_file(file_path, delimiter=','):
    with open(file_path, 'r') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        return list(reader)

# Function to parse FASTA files
def parse_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                sequences[header] = next(f).strip()
    return sequences

# Load data (Replace these paths with your actual file paths)
original_sequences = parse_fasta('chimeric_input.fasta')
vsearch_filtered_chimeras = parse_fasta('vsearch_long_chimeric.fasta')
chimera_info = parse_delimited_file('chimera_info.tsv', delimiter='\t')

# Identify sequence ID key
seq_id_key = list(chimera_info[0].keys())[0]

# Extract chimera IDs
chimera_ids = set(entry[seq_id_key] for entry in chimera_info)

# Initialize counters
TP, FP, TN, FN = 0, 0, 0, 0

# Compute metrics
for header in original_sequences.keys():
    is_chimera = header in chimera_ids
    is_filtered = header in vsearch_filtered_chimeras

    if is_chimera and is_filtered:
        TP += 1
    elif is_chimera and not is_filtered:
        FN += 1
    elif not is_chimera and is_filtered:
        FP += 1
    elif not is_chimera and not is_filtered:
        TN += 1

# Calculate Sensitivity and Specificity
sensitivity = TP / (TP + FN) if TP + FN > 0 else 'Undefined'
specificity = TN / (TN + FP) if TN + FP > 0 else 'Undefined'

# Results
result_dict = {
    'True Positives': TP,
    'False Positives': FP,
    'True Negatives': TN,
    'False Negatives': FN,
    'Sensitivity': sensitivity,
    'Specificity': specificity
}

print("Calculated Metrics:", result_dict)
