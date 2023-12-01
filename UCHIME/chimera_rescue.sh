#!/bin/bash

input_dir="/home/ali/Documents/simulated_data/analysis/Blast/treat_chimera/chimeras/combined_chimeras"
output_dir="/home/ali/Documents/simulated_data/analysis/Blast/treat_chimera/chimeras/combined_chimeras/rescued_chimeras"
min_occurrence=2

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Temporary files
temp_file=$(mktemp)
headers_file=$(mktemp)

# Process each FASTA file
for file in "$input_dir"/*.fasta; do
    sequence=""
    header=""
    while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ $line == ">"* ]]; then
            if [ -n "$sequence" ]; then
                # Save the previous sequence and header
                echo "$sequence" >> "$temp_file"
                echo "$header###$sequence" >> "$headers_file"
            fi
            header=$line
            sequence=""
        else
            sequence+=$line
        fi
    done < "$file"
    # Save the last sequence in the file
    if [ -n "$sequence" ]; then
        echo "$sequence" >> "$temp_file"
        echo "$header###$sequence" >> "$headers_file"
    fi
done

# Sort and count sequences
sort "$temp_file" | uniq -c | while read -r count sequence; do
    if [ "$count" -ge "$min_occurrence" ]; then
        # Find all headers for this sequence
        grep "###$sequence$" "$headers_file" | sed 's/###.*$//' | while read -r header; do
            echo -e "$header\n$sequence" >> "$output_dir/RescuedChimera.fasta"
        done
    fi
done

# Clean up
rm "$temp_file" "$headers_file"

echo "Processing complete. Output saved to $output_dir/RescuedChimera.fasta"
