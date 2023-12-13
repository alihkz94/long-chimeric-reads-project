#!/bin/bash

# Define directories
input_dir="/media/ali/data_store/test_chim/chimeras"
non_chimeric_dir="/media/ali/data_store/test_chim"
new_non_chimeric_dir="/media/ali/data_store/test_chim/new_non_chimeric"

# Create new directory for modified non-chimeric files
mkdir -p "$new_non_chimeric_dir"

# Report file
report_file="report.txt"

# Minimum sequence occurrence
min_occurrence=2

# Initialize report file
echo "Filename,Before,After,Rescued" > "$report_file"
total_rescued=0

# Process each chimeric file
for chimera_file in "$input_dir"/*.fasta; do
    basename=$(basename "$chimera_file" .chimeras.fasta)
    non_chimeric_file="$non_chimeric_dir/$basename.fasta"
    new_non_chimeric_file="$new_non_chimeric_dir/$basename.fasta"

    # Copy the original non-chimeric file to the new directory
    cp "$non_chimeric_file" "$new_non_chimeric_file"

    # Count sequences in the original non-chimeric file
    count_before=$(grep -c "^>" "$non_chimeric_file")

    # Initialize counter for rescued sequences
    rescued=0

    # Process chimeric sequences
    if [ -f "$chimera_file" ] && [ -s "$chimera_file" ]; then
        # Temporary file to store unique sequences from the chimeric file
        temp_chimera_file="/tmp/temp_chimera_$basename.fasta"
        touch "$temp_chimera_file"

        grep -v "^>" "$chimera_file" | sort | uniq -c | while read count sequence; do
            if [ "$count" -ge "$min_occurrence" ]; then
                # Find the header for the sequence
                header=$(grep -B 1 "$sequence" "$chimera_file" | grep "^>")
                # Add sequence to the new non-chimeric file
                echo -e "$header\n$sequence" >> "$new_non_chimeric_file"
                ((rescued++))
            fi
        done

        # Remove temporary file
        rm "$temp_chimera_file"
    fi

    # Count sequences in the enhanced non-chimeric file
    count_after=$(grep -c "^>" "$new_non_chimeric_file")
    # Calculate the number of rescued sequences
    rescued=$((count_after - count_before))

    # Write to report file
    echo "$basename,$count_before,$count_after,$rescued" >> "$report_file"
    total_rescued=$((total_rescued + rescued))
done

# Report total rescued sequences
echo "Total Rescued Sequences: $total_rescued" >> "$report_file"

echo "Processing complete. Chimeric sequences processed into non-chimeric files."
