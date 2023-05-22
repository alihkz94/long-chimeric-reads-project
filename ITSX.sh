#!/bin/bash
#SBATCH --job-name="ITSX"
#SBATCH --cpus-per-task=128
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --partition amd
#SBATCH --time=150:00:00


# Create the main output directory if it doesn't exist
mkdir -p itsx_out

# Loop over all FASTA files in the current directory
for fasta in *.fasta
do
  # Create a subdirectory for the current file
  subdir="itsx_out/${fasta%.*}"
  mkdir -p "$subdir"

  # Run ITSx with the specified settings and save the output in the subdirectory
  ITSx \
  --cpu 128 \
  -i "$fasta" \
  -o "$subdir/${fasta%.*}" \
  -nhmmer T \
  --complement T \
  --only_full T \
  -domains 2 \
  -score 0 \
  --partial 50 
done
