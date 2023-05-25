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

# Define the input and output files
input_file = 'input.fasta'
output_file = 'output.fasta'

# Define the base command
base_cmd = f'vsearch --uchime_denovo {input_file} --abskew 1.5 --nonchimeras {output_file}'

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

# Plot the results
for param, values in parameters.items():
    plt.plot(values, results[param], label=param)
plt.xlabel('Parameter Value')
plt.ylabel('Number of Chimeras')
plt.legend()
plt.title('Effect of VSEARCH Parameters on Chimera Detection')
plt.savefig('vsearch_results.png')
