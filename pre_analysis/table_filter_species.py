# FILTER THE TABLE FOR GETTING AN OUTPUT WITH SPECIES

import pandas as pd

# Read the CSV file into a pandas DataFrame without column names
df = pd.read_csv('chimera.csv', header=None)

# Drop rows with any empty cell in the DataFrame
df.dropna(inplace=True)

# Filter rows where all 10 columns have non-null values
filtered_df = df[df.notnull().all(axis=1)]

# Save the filtered DataFrame to a new CSV file
filtered_df.to_csv('filtered_file.csv', index=False, header=False)




##Generate FASTA files with species names on it: 

import csv

fasta_file_path = "fungi.fasta"
csv_file_path = "fasta.csv"
output_file_path = "fungal_species.fasta"

# Load the CSV table into a dictionary
species_dict = {}
with open(csv_file_path, "r") as csv_file:
    reader = csv.reader(csv_file)
    for row in reader:
        fasta_header = row[0]
        species_name = row[1]
        species_dict[fasta_header] = species_name

# Process the FASTA file and create the new file with updated headers
with open(fasta_file_path, "r") as fasta_file, open(output_file_path, "w") as output_file:
    current_header = None
    for line in fasta_file:
        if line.startswith(">"):
            fasta_header = line.strip()[1:]
            if fasta_header in species_dict:
                current_header = fasta_header
                species_name = species_dict[fasta_header]
                new_header = f">{fasta_header}|{species_name}\n"
                output_file.write(new_header)
            else:
                current_header = None
        elif current_header is not None:
            output_file.write(line)
