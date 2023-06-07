#!/bin/bash

output_file="metabar_uchime_input.fasta"

# Remove the output file if it already exists
if [ -f "$output_file" ]; then
    rm "$output_file"
fi

# Loop through all the Fasta files in the folder
for input_file in *.fasta; do
    awk '/^>/ {if (NR>1) printf("\n"); printf("%s\n", $0); next;} {printf("%s", $0);} END {printf("\n");}' "$input_file" >> "$output_file"
done
