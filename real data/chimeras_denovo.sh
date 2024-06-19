#!/bin/bash

# Directory containing your FASTA files
FASTA_DIR="/home/ali/Documents/simulated_data/analysis/tagJump_filt/original"

# Create folders for each step
mkdir -p dereplicated chimeras
mkdir -p chimeras/combined_chimeras  # Create subfolder for combined chimeras

# Loop through each FASTA file in the directory
for FASTA_FILE in "$FASTA_DIR"/*.fasta; do
    BASENAME=$(basename "$FASTA_FILE" .fasta)

    # 1. Dereplication
    echo "Dereplication for $BASENAME..."
    vsearch --derep_fulllength "$FASTA_FILE" --sizein --sizeout --fasta_width 0 --output "dereplicated/${BASENAME}_dereplicated.fasta"

    # 2. Chimera Filtering with UCHIME Denovo
    echo "Chimera Filtering with UCHIME Denovo for $BASENAME..."
    vsearch --chimeras_denovo "dereplicated/${BASENAME}_dereplicated.fasta" --abskew 3 \
    --nonchimeras "chimeras/${BASENAME}_nonchimeras.fasta" --chimeras "chimeras/${BASENAME}_chimeras.fasta" 

echo "Script completed."
