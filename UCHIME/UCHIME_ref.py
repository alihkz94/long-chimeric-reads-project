import os
import subprocess
import matplotlib.pyplot as plt
import pandas as pd

# Define the parameters
parameters = {
    '--dn': [1.4, 1.5, 1.6, 1.7, 1.8],
    '--mindiffs': [1, 2, 3],
    '--minh': [0.1, 0.2, 0.28],
    '--xn': [6.0, 7.0, 8.0, 9.0],
    '--mindiv': [0.6, 0.7, 0.8, 0.9],
    '--abskew': [1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
}

# Define the input and output files
input_file = 'metabar_uchime_input.fasta'
db_file = 'ITS.fasta'

# Define the base command
base_cmd = f'vsearch --uchime_ref {input_file} --threads 8 --db {db_file}'

def parse_chimera_info(filename):
    # Define the column names
    column_names = ["chimera_id", "seq1_id", "seq1_file", "seq2_id", "seq2_file", "breakpoint", \
                    "reversed_status", "ratio", "chimera_seq_length"]

    # Load the data
    data = pd.read_csv(filename, sep="\t", header=None, names=column_names)

    # Convert the DataFrame to a dictionary
    data_dict = data.to_dict(orient="records")

    return data_dict

# Loop over the parameters
results = {}
for param, values in parameters.items():
    results[param] = []
    for value in values:
        # Construct the command
        output_chimera = f'chimeric_{param}_{value}.fasta'
        output_nonchimera = f'non_chimer_{param}_{value}.fasta'
        cmd = f'{base_cmd} {param} {value} --chimeras {output_chimera} --nonchimeras {output_nonchimera}'

        # Run the command
        process = subprocess.run(cmd, shell=True)

        # Parse the output to get the number of chimeras
        chimera_info_file = os.path.join('output_directory', f"chimera_info_{param}_{value}.tsv")
        if os.path.exists(chimera_info_file):
            chimera_dict = parse_chimera_info(chimera_info_file)
            num_chimeras = len(chimera_dict)
        else:
            num_chimeras = 0  # If the file does not exist, assume no chimeras were found

        # Store the result
        results[param].append(num_chimeras)
