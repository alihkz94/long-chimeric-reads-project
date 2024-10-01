#!/bin/bash

# Directory containing the FASTA files
FASTA_DIR="."

# Directory to store the BLAST databases
DB_DIR="./databases"

# Create BLAST databases
for fasta_file in ${FASTA_DIR}/*.fasta; do
    # Get the base name of the file (without extension)
    base_name=$(basename ${fasta_file} .fasta)
    
    # Create a directory for this database
    mkdir -p "${DB_DIR}/${base_name}"
    
    # Create the BLAST database
    makeblastdb -in "${fasta_file}" -dbtype nucl -out "${DB_DIR}/${base_name}/${base_name}"
    
    echo "Created BLAST database for ${base_name}"
done

echo "All BLAST databases have been created."