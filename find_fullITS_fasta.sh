#!/bin/bash

# Define the source and destination directories
source_dir="/media/ali/data_store/itsx_out"
destination_dir="/media/ali/data_store/full_ITS"

# Create the destination directory if it doesn't exist
mkdir -p "$destination_dir"

# Find all FASTA files in subfolders of the source directory with ".full_and_partial" in their names
# but excluding "ITS1.full_and_partial.fasta" and "ITS2.full_and_partial.fasta"
find "$source_dir" -type f -regex ".*/[^/]*\.full_and_partial\.fasta" \
  ! -regex ".*/[^/]*ITS1\.full_and_partial\.fasta" \
  ! -regex ".*/[^/]*ITS2\.full_and_partial\.fasta" | while read -r fasta_file; do
  # Copy the file to the destination directory
  cp "$fasta_file" "$destination_dir"
done

echo "Files copied to $destination_dir."
