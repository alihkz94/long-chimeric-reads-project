import os
import sys
from Bio import SeqIO
import numpy as np
import pandas as pd
import csv


#load the fasta file
fasta_sequences = list(SeqIO.parse(open("ITS.fasta"),'fasta'))

# get lengths of sequences and species names for later use
lengths = []
species_names = []
values = []

for i, fasta in enumerate(fasta_sequences):
    name, sequence = fasta.id, str(fasta.seq)
    lengths.append(len(sequence))
    
    # get species name from header
    species_name = "_".join(name.split("|")[1:]).replace(" ", "_")
    species_names.append(f"{species_name}_{i}")

    # generate value for simlord -n parameter
    value = int(1000 + ((9000 / (len(fasta_sequences) - 1)) * i))
    values.append(value)

# sort lengths, species names and values based on the lengths
sorted_indices = np.argsort(lengths)
lengths = np.array(lengths)[sorted_indices]
species_names = np.array(species_names)[sorted_indices]
values = np.array(values)[sorted_indices]

# create a dataframe to store the results
df = pd.DataFrame(columns=['File', 'Number of Sequences'])

# iterate over each sequence and run the simlord command
for i in range(len(fasta_sequences)):
    temp_file_name = f"temp_{species_names[i]}.fasta"
    new_file_name = species_names[i]
    value = values[i]
    
    # generate temp fasta file
    with open(temp_file_name, "w") as temp_file:
        SeqIO.write(fasta_sequences[sorted_indices[i]], temp_file, "fasta")
    
    # run the simlord command
    os.system(f"simlord -n {value} --read-reference {temp_file_name} {new_file_name} --max-passes 30")

    # delete the temporary fasta file
    os.remove(temp_file_name)


#To check the number of sequences and generate a table for it:
 
folder_path = '/home/ali/Documents/simulated_data/fasta_generation/simluate_fasta'
output_csv = '/home/ali/Documents/simulated_data/fasta_generation/num_seq_ITS.csv'

# List all files in the folder
file_list = os.listdir(folder_path)

# Initialize the list to store file names and sequence counts
data = [['File Name', 'Sequence Count']]

# Process each file
for file_name in file_list:
    # Check if the file is a FASTQ file
    if file_name.endswith('.fastq'):
        file_path = os.path.join(folder_path, file_name)
        
        # Count the number of sequences in the file
        seq_count = sum(1 for _ in SeqIO.parse(file_path, 'fastq'))
        
        # Add the file name and sequence count to the data list
        data.append([file_name, seq_count])

# Write the data to a CSV file
with open(output_csv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(data)
print('CSV table generated successfully.')
