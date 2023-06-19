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
    '--mindiv': [0.6, 0.7, 0.8, 0.9]
}

# Define the input and output files
input_file = 'metabar_uchime_input.fasta'

# Define the base command
base_cmd = f'vsearch --uchime_denovo {input_file} --abskew 2'

def parse_chimera_info(filename):
    # Define the column names
    column_names = ["chimera_id", "seq1_id", "seq1_file", "seq2_id", "seq2_file", "breakpoint", "reversed_status", "ratio", "chimera_seq_length"]

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
        output_file = f'output_{param}_{value}.fasta'
        cmd = f'{base_cmd} {param} {value} --nonchimeras {output_file}'

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

# Plot the results
for param, values in parameters.items():
    plt.plot(values, results[param], label=param)
plt.xlabel('Parameter Value')
plt.ylabel('Number of Chimeras')
plt.legend()
plt.title('Effect of VSEARCH Parameters on Chimera Detection')
plt.savefig('vsearch_results.png')
