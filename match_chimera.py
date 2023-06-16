from Bio import SeqIO

# Function to parse the chimera_info_file and store it in a dictionary
def parse_chimera_info(chimera_info_file):
    chimera_dict = {}
    with open(chimera_info_file, "r") as f:
        next(f)  # skip header
        for line in f:
            chimera_id, seq1_id, seq2_id, breakpoint = line.strip().split("\t")
            chimera_dict[chimera_id] = (seq1_id, seq2_id, int(breakpoint))
    return chimera_dict

# Function to parse the FASTA file of detected chimeric sequences and store the IDs in a set
def parse_detected_chimeras(fasta_file):
    detected_chimeras = set()
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            detected_chimeras.add(record.id.split(";")[0])  # Remove ";size=2" from the ID
    return detected_chimeras

chimera_info_file = (
    "/home/ali/Documents/simulated_data/chimera/filtering/"
    "chimera_Filtered_out/chimeras/chimera_info.txt"
)

detected_chimeras_fasta = (
    "/home/ali/Documents/simulated_data/chimera/filtering/"
    "chimera_Filtered_out/chimeras/ch_coi20.denovo.chimeras.fasta"
)

# Parse the chimera_info_file and detected_chimeras_fasta
chimera_dict = parse_chimera_info(chimera_info_file)
detected_chimeras = parse_detected_chimeras(detected_chimeras_fasta)

# Compare the generated and detected chimeras
match_count = 0
for chimera_id in chimera_dict:
    if chimera_id in detected_chimeras:
        match_count += 1

# Calculate the proportion of generated chimeras that were successfully detected
proportion_detected = match_count / len(chimera_dict)

print(f"Generated chimeras: {len(chimera_dict)}")
print(f"Detected chimeras: {match_count}")
print(f"Proportion of generated chimeras detected: {proportion_detected:.2%}")


#### Report generation ####

total_non_chimera_reads = 1765  # Set this value to the number of non-chimera reads in your dataset

# Function to calculate true positives, false positives, true negatives, and false negatives
def calculate_rates(chimera_dict, detected_chimeras, total_non_chimera_reads):
    true_positives = len(detected_chimeras.intersection(chimera_dict.keys()))
    false_positives = len(detected_chimeras.difference(chimera_dict.keys()))
    false_negatives = len(chimera_dict.keys()) - true_positives
    true_negatives = total_non_chimera_reads - false_positives

    return true_positives, false_positives, true_negatives, false_negatives

# Calculate rates
tp, fp, tn, fn = calculate_rates(chimera_dict, detected_chimeras, total_non_chimera_reads)

print(f"True Positives: {tp}")
print(f"False Positives: {fp}")
print(f"True Negatives: {tn}")
print(f"False Negatives: {fn}")
