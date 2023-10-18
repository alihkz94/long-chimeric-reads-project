from Bio import SeqIO
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_fscore_support

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

# Read metabar FASTA file
metabar_dict = SeqIO.to_dict(handle_duplicate_keys(SeqIO.parse("metabar_uchime_input.fasta", "fasta")))

# Read removed chimeric sequences
with open("removed_chimeric_sequences_8.txt", 'r') as f:
    removed_seqs = f.read().splitlines()

# Initialize counts
TP = 0  # True Positives
FP = 0  # False Positives
TN = 0  # True Negatives
FN = 0  # False Negatives

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

# Calculate Sensitivity and Specificity
sensitivity = TP / (TP + FN)
specificity = TN / (TN + FP)

print(f"True Positives: {TP}")
print(f"False Positives: {FP}")
print(f"True Negatives: {TN}")
print(f"False Negatives: {FN}")
print(f"Sensitivity: {sensitivity}")
print(f"Specificity: {specificity}")

# Use scikit-learn to calculate precision, recall, and F1-score
y_true = [1]*TP + [0]*FP + [1]*FN + [0]*TN
y_pred = [1]*(TP + FP) + [0]*(FN + TN)
precision, recall, f1_score, _ = precision_recall_fscore_support(y_true, y_pred, average='binary')

print(f"Precision: {precision}")
print(f"Recall: {recall}")
print(f"F1 Score: {f1_score}")

# Confusion Matrix
cm = confusion_matrix(y_true, y_pred)
print("Confusion Matrix:")
print(cm)
