import os
import subprocess
import matplotlib.pyplot as plt

# Define the parameters
parameters = {
    '--dn': [1.4, 1.5, 1.6, 1.7, 1.8],
    '--mindiffs': [1, 2, 3],
    '--minh': [0.1, 0.2, 0.28],
    '--xn': [6.0, 7.0, 8.0, 9.0],
    '--mindiv': [0.6, 0.7, 0.8, 0.9]
}

# Define the input file
input_file = 'metabar_uchime_input.fasta'

# Define the base command
base_cmd = f'vsearch --uchime_denovo {input_file} --abskew 1.5 --nonchimeras'

# Loop over the parameters
results = {}
for param, values in parameters.items():
    results[param] = []
    for value in values:
        # Construct the command
        output_dir = f"output_{param.replace('--','')}_value_{value}"  # use parameter name and value to create distinct directory names
        os.makedirs(output_dir, exist_ok=True)  # create the output directory
        output_file = os.path.join(output_dir, 'output.fasta')
        param_file = os.path.join(output_dir, 'params.txt')
        cmd = f'{base_cmd} {output_file} {param} {value}'
        
        # Run the command
        print(f"Running command: {cmd}")
        process = subprocess.run(cmd, shell=True, capture_output=True)
        
        # Parse the output to get the number of chimeras
        output = process.stdout.decode()
        try:
            num_chimeras = int(output.split('\n')[-2].split()[0])
        except IndexError:
            num_chimeras = 0  # assign a default value if the expected output is not found
        
        # Store the result
        results[param].append(num_chimeras)
        
        # Save the used parameters to a file
        with open(param_file, "w") as f:
            f.write(f"{param} {value}\n")
