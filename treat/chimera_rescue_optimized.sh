## The script is designed to add True nonchimeric reads from the chimeric files into the non-chimeric ones. For the input, three directories are needed. One for the chimeric reads which needs the name of the 
## fasta files be in ".chimeras.fasta" format. The other is for the non_chimeric directory with ".fasta" format files, and later on, the "Rechime_non_chimeric" directory will be created for storing the new 
## non_chimeric ones. The report.txt file will be created at the end which holds the number of the sequences related to each file in three states of "before, after, rescued". The file "rechimed.pool.fasta"
## will be created to store all the new nonchimeric sequences added to the different samples.
## usage:  ./chimera_rescue.sh -c ./chimeras -n . 
## "-c" option is for the chimeric reads directory and the "-n" option is for the nonchimeric reads directory.

#!/bin/bash

# Function to check if a directory exists
check_dir() {
    if [ ! -d "$1" ]; then
        echo "Error: Directory $1 does not exist."
        exit 1
    fi
}

# Function to validate FASTA format
validate_fasta() {
    if [ ! -s "$1" ] || ! grep -q "^>" "$1"; then
        echo "Error: File $1 does not appear to be in FASTA format."
        exit 1
    fi
}

# Initialize variables for directory paths
chimeric_dir=""
non_chimeric_dir=""

# Parse command-line options
while getopts ":c:n:" opt; do
  case $opt in
    c) chimeric_dir="$OPTARG"
    ;;
    n) non_chimeric_dir="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done

# Check if directories were provided
if [ -z "$chimeric_dir" ] || [ -z "$non_chimeric_dir" ]; then
    echo "Usage: $0 -c <chimeric_reads_dir> -n <non_chimeric_dir>"
    exit 1
fi

# Automatically set the new_non_chimeric_dir
new_non_chimeric_dir="${non_chimeric_dir}/Rechime_non_chimeric"

# Create a subfolder for rescued sequences
seq_rescued_dir="${new_non_chimeric_dir}/seq_rescued"

# Path for the pooled rescued sequences file
rechimed_pool_file="${seq_rescued_dir}/rechimed.pool.fasta"

# Remove the existing new_non_chimeric_dir if it exists and create a new one, including the seq_rescued subdirectory
if [ -d "$new_non_chimeric_dir" ]; then
    rm -rf "$new_non_chimeric_dir"
fi
mkdir -p "$seq_rescued_dir"

# Check if directories exist
check_dir "$chimeric_dir"
check_dir "$non_chimeric_dir"

# Preprocess all FASTA files to ensure sequences are not wrapped and parsed correctly
for fasta_file in "$chimeric_dir"/*.fasta "$non_chimeric_dir"/*.fasta; do
    seqkit seq -w 0 "$fasta_file" > "${fasta_file}.tmp" && mv "${fasta_file}.tmp" "$fasta_file"
done

# Report file path
report_file="${new_non_chimeric_dir}/report.txt"

# Minimum sequence occurrence
min_occurrence=2

# Initialize report file
echo "Filename,Before,After,Rescued" > "$report_file"
total_rescued=0

# Main analysis logic
for chimera_file in "$chimeric_dir"/*.chimeras.fasta; do
    # Skip validation for empty chimeric files
    if [ ! -s "$chimera_file" ]; then
        echo "Notice: Skipping empty file $chimera_file."
        continue
    fi

    # Validate FASTA format
    validate_fasta "$chimera_file"

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
    declare -A seq_header_map
    while IFS= read -r line; do
        if [[ $line == ">"* ]]; then
            current_header=$line
        else
            seq_header_map["$line"]+="${current_header}\n"
        fi
    done < "$chimera_file"

    for sequence in "${!seq_header_map[@]}"; do
        count=$(echo -e "${seq_header_map[$sequence]}" | grep -c "^>")
        if [ "$count" -ge "$min_occurrence" ]; then
            header=$(echo -e "${seq_header_map[$sequence]}" | head -n 1)
            echo -e "$header\n$sequence" >> "$new_non_chimeric_file"
            echo -e "$header\n$sequence" >> "$rechimed_pool_file"
            ((rescued++))
        fi
    done

    # Count sequences in the enhanced non-chimeric file and report
    count_after=$(grep -c "^>" "$new_non_chimeric_file")
    echo "$basename,$count_before,$count_after,$rescued" >> "$report_file"
    total_rescued=$((total_rescued + rescued))
done

# Final report
echo "Total Rescued Sequences: $total_rescued" >> "$report_file"
echo "Processing complete. Chimeric sequences processed into non-chimeric files."
