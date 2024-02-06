## The script is designed to add True nonchimeric reads from the chimeric files into the non-chimeric ones. For the input, three directories are needed. One for the chimeric reads which needs the name of the 
## fasta files be in ".chimeras.fasta" format. The other is for the non_chimeric directory with ".fasta" format files, and later on, the "new_non_chimeric" directory will be created for storing the new 
## non_chimeric ones. The report.txt file will be created at the end which holds the number of the sequences related to each file in three states of "before, after, rescued". 

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
    if ! grep -q "^>" "$1"; then
        echo "Error: File $1 does not appear to be in FASTA format."
        exit 1
    fi
}

# Function to validate nucleotide sequences in a FASTA file
validate_sequence() {
    if ! grep -v "^>" "$1" | grep -qvE "^[ACGTNacgtn]+$"; then
        echo "Error: Invalid nucleotide sequence detected in file $1."
        exit 1
    fi
}

# Function to clean up temporary files
cleanup() {
    echo "Cleaning up temporary files..."
    rm -f /tmp/temp_chimera_*.fasta
}

# Trap to execute the cleanup function on script exit
trap cleanup EXIT

# Define directories
input_dir="/media/ali/data_store/test_chim/chimeras"
non_chimeric_dir="/media/ali/data_store/test_chim"
new_non_chimeric_dir="/media/ali/data_store/test_chim/new_non_chimeric"

# Create a new directory for modified non-chimeric files
mkdir -p "$new_non_chimeric_dir"

# Check if directories exist
check_dir "$input_dir"
check_dir "$non_chimeric_dir"
check_dir "$new_non_chimeric_dir"

# Report file
report_file="report.txt"

# Minimum sequence occurrence
min_occurrence=2

# Initialize report file
echo "Filename,Before,After,Rescued" > "$report_file"
total_rescued=0

# Process each chimeric file
for chimera_file in "$input_dir"/*.fasta; do
    # Validate FASTA format and sequence integrity
    validate_fasta "$chimera_file"
    validate_sequence "$chimera_file"

    basename=$(basename "$chimera_file" .chimeras.fasta)
    non_chimeric_file="$non_chimeric_dir/$basename.fasta"
    new_non_chimeric_file="$new_non_chimeric_dir/$basename.fasta"

    # Check for successful copy
    if ! cp "$non_chimeric_file" "$new_non_chimeric_file"; then
        echo "Error: Failed to copy file $non_chimeric_file."
        exit 1
    fi

    # Count sequences in the original non-chimeric file
    count_before=$(grep -c "^>" "$non_chimeric_file")
    if [ $? -ne 0 ]; then
        echo "Error: Failed to count sequences in $non_chimeric_file."
        exit 1
    fi

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
                if [ $? -ne 0 ]; then
                    echo "Error: Failed to find header for a sequence in $chimera_file."
                    continue
                fi
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
    if [ $? -ne 0 ]; then
        echo "Error: Failed to count sequences in $new_non_chimeric_file."
        exit 1
    fi

    # Calculate the number of rescued sequences
    rescued=$((count_after - count_before))

    # Write to report file
    echo "$basename,$count_before,$count_after,$rescued" >> "$report_file"
    total_rescued=$((total_rescued + rescued))
done

# Report total rescued sequences
echo "Total Rescued Sequences: $total_rescued" >> "$report_file"

echo "Processing complete. Chimeric sequences processed into non-chimeric files."
