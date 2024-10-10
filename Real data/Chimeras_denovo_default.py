import subprocess
import os

# Directory containing the merged FASTA files
input_dir = '/mnt/d/vsearch_default'

# Directories to store the output files
chimeras_dir = '/mnt/d/vsearch_default/chimeras'
nonchimeras_dir = '/mnt/d/vsearch_default/nonchimeras'

# Create output directories if they do not exist
if not os.path.exists(chimeras_dir):
    os.makedirs(chimeras_dir)
if not os.path.exists(nonchimeras_dir):
    os.makedirs(nonchimeras_dir)

# Iterate over FASTA files in the input directory
for fasta_file in os.listdir(input_dir):
    if fasta_file.endswith('.fasta'):
        input_file_path = os.path.join(input_dir, fasta_file)
        
        # Construct output file paths
        chimeras_file_path = os.path.join(chimeras_dir, fasta_file.replace('.fasta', '_chimeras.fasta'))
        nonchimeras_file_path = os.path.join(nonchimeras_dir, fasta_file.replace('.fasta', '_nonchimeras.fasta'))
        
        # Construct the vsearch command
        vsearch_cmd = [
            'vsearch',
            '--chimeras_denovo', input_file_path,
            '--sizein', '--sizeout',
            '--fasta_width', '0',
            '--chimeras', chimeras_file_path,
            '--nonchimeras', nonchimeras_file_path
        ]

        # Execute the command
        subprocess.run(vsearch_cmd)

print("Chimera filtering completed. Check the 'chimeras' and 'nonchimeras' directories for the output.")
