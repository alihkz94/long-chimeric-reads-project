#!/bin/bash

# Create the "dereplicated" folder if it doesn't exist
mkdir -p dereplicated

# Loop through all files ending with ".fasta" in the current directory
for file in *.fasta; do
  # Extract the filename without extension
  filename="${file%.*}"

  # Run vsearch with the desired options
  vsearch --derep_fulllength "$file" --sizein --sizeout --fasta_width 0 --output "dereplicated/$filename.derep.fasta"

  # Print a success message
  echo "Dereplicated '$file' to 'dereplicated/$filename.derep.fasta'"
done

echo "All FASTA files dereplicated!"
