import os
import subprocess
import pandas as pd

# Create tempdir if it doesn't exist
if not os.path.exists('tempdir'):
    os.makedirs('tempdir')

# Define the parameters
parameters = {
    '--minh': [0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01]
}

# Define the input and output files
input_file = 'metabar_uchime_input.fasta'
id_value = 0.97  # Modify as per your requirements

def parse_chimera_info(filename):
    column_names = ["chimera_id", "seq1_id", "seq1_file", "seq2_id", "seq2_file", "breakpoint", "reversed_status", "ratio", "chimera_seq_length"]
    data = pd.read_csv(filename, sep="\t", header=None, names=column_names)
    return data.to_dict(orient="records")

# Dereplicate and pre-cluster sequences
dereplicate_cmd = f'vsearch --derep_fulllength {input_file} --sizein --sizeout --fasta_width 0 --output tempdir/{input_file}.derep.fasta'
precluster_cmd = f'vsearch --cluster_size tempdir/{input_file}.derep.fasta --id {id_value} --sizein --sizeout --fasta_width 0 --centroids tempdir/{input_file}.preclustered.fasta'

subprocess.run(dereplicate_cmd, shell=True)
subprocess.run(precluster_cmd, shell=True)

# Loop over the parameters
results = {}
for param, values in parameters.items():
    results[param] = []
    for value in values:
        output_file = f'output_{param}_{value}.fasta'
        uchime_cmd = f'vsearch --uchime_denovo tempdir/{input_file}.preclustered.fasta {param} {value} --sizein --sizeout --fasta_width 0 \
            --mindiv 0.4 --dn 1.6 --chimeras tempdir/{output_file} --uchimeout tempdir/uchime_out_{param}_{value}.tsv --borderline tempdir/borderline_{param}_{value}.fasta'

        subprocess.run(uchime_cmd, shell=True)

        # Combine chimeric and borderline reads
        combined_file = f'tempdir/combined_{param}_{value}.fasta'
        with open(combined_file, 'w') as outfile:
            for fname in [f'tempdir/{output_file}', f'tempdir/borderline_{param}_{value}.fasta']:
                with open(fname) as infile:
                    outfile.write(infile.read())

        # Parse UCHIME output
        chimera_info_file = f'tempdir/uchime_out_{param}_{value}.tsv'
        if os.path.exists(chimera_info_file):
            chimera_dict = parse_chimera_info(chimera_info_file)
            num_chimeras = len(chimera_dict)
        else:
            num_chimeras = 0

        results[param].append(num_chimeras)
