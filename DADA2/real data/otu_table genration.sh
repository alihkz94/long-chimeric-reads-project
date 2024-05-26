#!/bin/bash

# Combine all FASTA files into one
cat *.fasta > combined.fasta

# Perform global dereplication 
vsearch --derep_fulllength combined.fasta --output dereplicated.fasta --sizeout --minuniquesize 1

# Use the dereplicated sequences as a database
vsearch --search_exact combined.fasta --db dereplicated.fasta --otutabout otu_table.txt --threads 8

echo "OTU table generated with exact sequence matches."

# Clean up intermediate files if necessary
rm combined.fasta

echo "Script completed successfully."

#rereplicate the file to be suitable as an input for the DADA2 script: 
vsearch --rereplicate dereplicated.fasta --output replicated.fasta
