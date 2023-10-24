from Bio import SeqIO
from sklearn.metrics import precision_recall_fscore_support
import pandas as pd
import glob
import concurrent.futures

# Function to handle duplicate keys
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

# Initialize DataFrame to store results
df = pd.DataFrame(columns=['File', 'True Positives', 'False Positives', 'True Negatives', 'False Negatives', 'Sensitivity', 'Specificity', 'Precision', 'Recall', 'F1 Score'])

# Read metabar FASTA file once
metabar_dict = SeqIO.to_dict(handle_duplicate_keys(SeqIO.parse("metabar_uchime_input.fasta", "fasta")))

# Function to process each file
def process_file(filepath):
    # Initialize counts
    TP, FP, TN, FN = 0, 0, 0, 0
    
    # Read removed chimeric sequences from current file and convert to set for fast look-up
    with open(filepath, 'r') as f:
        removed_seqs = set(f.read().splitlines())
    
    # Calculate TP, FP, TN, FN
    for seq_id, sequence in metabar_dict.items():
        if sequence.seq in removed_seqs:
            if "chimera" in seq_id:
                TP += 1
            else:
                FP += 1
        else:
            if "chimera" in seq_id:
                FN += 1
            else:
                TN += 1
    
    # Compute metrics
    sensitivity = TP / (TP + FN)
    specificity = TN / (TN + FP)
    y_true = [1]*TP + [0]*FP + [1]*FN + [0]*TN
    y_pred = [1]*(TP + FP) + [0]*(FN + TN)
    precision, recall, f1_score, _ = precision_recall_fscore_support(y_true, y_pred, average='binary')
    
    return {
        'File': filepath,
        'True Positives': TP,
        'False Positives': FP,
        'True Negatives': TN,
        'False Negatives': FN,
        'Sensitivity': sensitivity,
        'Specificity': specificity,
        'Precision': precision,
        'Recall': recall,
        'F1 Score': f1_score
    }

# Parallel processing of files
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_file, glob.glob("removed_chimeric_sequences_*.txt")))

# Add results to DataFrame
for result in results:
    df = df.append(result, ignore_index=True)

# Save DataFrame as a CSV file for the final report
df.to_csv("final_report.csv", index=False)
