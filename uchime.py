import subprocess
import matplotlib.pyplot as plt
from Bio import SeqIO

# Define the parameters
parameters = {
    '--dn': [1.4, 1.5, 1.6, 1.7, 1.8],
    '--mindiffs': [1, 2, 3],
    '--minh': [0.1, 0.2, 0.28],
    '--xn': [6.0, 7.0, 8.0, 9.0],
    '--mindiv': [0.6, 0.7, 0.8, 0.9]
}

# Define the input and output files
input_file = 'input.fasta'
output_file = 'output.fasta'

# Define the base command
base_cmd = f'vsearch --uchime_denovo {input_file} --abskew 1.5 --nonchimeras {output_file}'

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

# Function to calculate true positives, false positives, true negatives, and false negatives
def calculate_rates(chimera_dict, detected_chimeras, total_non_chimera_reads):
    true_positives = len(detected_chimeras.intersection(chimera_dict.keys()))
    false_positives = len(detected_chimeras.difference(chimera_dict.keys()))
    false_negatives = len(chimera_dict.keys()) - true_positives
    true_negatives = total_non_chimera_reads - false_positives

    return true_positives, false_positives, true_negatives, false_negatives

# Loop over the parameters
results = {}
for param, values in parameters.items():
    results[param] = []
    for value in values:
        # Construct the command
        cmd = f'{base_cmd} {param} {value}'
        
        # Run the command
        process = subprocess.run(cmd, shell=True, capture_output=True)
        
        # Parse the output to get the number of chimeras
        output = process.stdout.decode()
        num_chimeras = int(output.split('\n')[-2].split()[0])
        
        # Store the result
        results[param].append(num_chimeras)

        # Parse the chimera_info_file and detected_chimeras_fasta
        chimera_dict = parse_chimera_info('chimera_info.txt')
        detected_chimeras = parse_detected_chimeras(output_file)

        # Calculate rates
        tp, fp, tn, fn = calculate_rates(chimera_dict, detected_chimeras, total_non_chimera_reads)

        print(f"For {param} = {value}:")
        print(f"True Positives: {tp}")
        print(f"False Positives: {fp}")
        print(f"True Negatives: {tn}")
        print(f"False Negatives: {fn}")

# Plot the results
for param, values in parameters.items():
    plt.plot(values, results[param], label=param)
plt.xlabel('Parameter Value')
plt.ylabel('Number of Chimeras')
plt.legend()
plt.title('Effect of VSEARCH Parameters on Chimera Detection')
plt.savefig('vsearch_results.png')
