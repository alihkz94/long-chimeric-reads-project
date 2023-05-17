import os
from Bio import SeqIO
import numpy as np

# load the fasta file
fasta_sequences = list(SeqIO.parse(open("fungal_species.fasta"),'fasta'))

# get lengths of sequences and species names for later use
lengths = []
species_names = []

for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    lengths.append(len(sequence))
    
    # get species name from header
    species_name = name.split("|")[-1].replace(" ", "_")
    species_names.append(species_name)

# sort lengths and species names based on the lengths
sorted_indices = np.argsort(lengths)
lengths = np.array(lengths)[sorted_indices]
species_names = np.array(species_names)[sorted_indices]

# generate values for simlord -n parameter
values = np.linspace(1000, 10000, num=len(lengths), dtype=int)

# iterate over each sequence and run the simlord command
for i in range(len(lengths)):
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
