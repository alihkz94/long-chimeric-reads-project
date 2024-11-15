import os
import subprocess
import pandas as pd

# Create a temporary directory if it doesn't exist
if not os.path.exists('tempdir'):
    os.makedirs('tempdir')

# Define the parameters for --chimeras_denovo
parameters = {
    '--abskew': [2,3,4],
    '--chimeras_diff_pct': [0.5, 0.6, 0.7, 0.8, 0.9],
    '--chimeras_parents_max': [2,3],
    '--chimeras_length_min' : [10,20,30],
    '--chimeras_parts': [2,3,4]
}

# Define input file
input_file = 'metabar_uchime_input.fasta'  # Modify with your actual input file
def parse_chimera_info(filename):
    # Define column names for the chimera output
    column_names = ['chimera_id', 'parent_seq1', 'parent_seq2', 'breakpoint', 'support', 'other_info']
    data = pd.read_csv(filename, sep="\t", header=None, names=column_names)
    return data.to_dict(orient="records")

# Dereplicate sequences
dereplicate_cmd = f'vsearch --derep_fulllength {input_file} --sizein --sizeout --fasta_width 0 --output tempdir/{input_file}.derep.fasta'

# Run the dereplication command
subprocess.run(dereplicate_cmd, shell=True)

# Loop over the parameters
results = {}
for param, values in parameters.items():
    results[param] = []
    for value in values:
        # Construct the vsearch --chimeras_denovo command
        chimeras_file = f'tempdir/chimeras_{param}_{value}.fasta'
        nonchimeras_file = f'tempdir/nonchimeras_{param}_{value}.fasta'
        vsearch_cmd = f'vsearch --chimeras_denovo tempdir/{input_file}.derep.fasta {param} {value} --sizein --sizeout --fasta_width 0 \
            --chimeras {chimeras_file} --nonchimeras {nonchimeras_file}'

        # Run the vsearch command
        subprocess.run(vsearch_cmd, shell=True)

        # Parse the output to get the number of chimeras
        if os.path.exists(chimeras_file):
            chimera_dict = parse_chimera_info(chimeras_file)
            num_chimeras = len(chimera_dict)
        else:
            num_chimeras = 0  # If the file does not exist, assume no chimeras were found

        # Store the results
        results[param].append((value, num_chimeras))
