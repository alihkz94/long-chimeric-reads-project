###### WITH DEREPLICATION ########

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
    vsearch --uchime_denovo "dereplicated/${BASENAME}_dereplicated.fasta" --mindiv 0.4 --dn 1.6 --minh 0.08 --sizein --sizeout --fasta_width 0 \
    --nonchimeras "chimeras/${BASENAME}_denovo_nonchimeras.fasta" --chimeras "chimeras/${BASENAME}_denovo_chimeras.fasta" --borderline "chimeras/${BASENAME}_denovo_borderline.fasta"

    # Combining chimeric and borderline sequences for each file
    echo "Combining chimeric and borderline sequences for $BASENAME..."
    cat "chimeras/${BASENAME}_denovo_chimeras.fasta" "chimeras/${BASENAME}_denovo_borderline.fasta" > "chimeras/combined_chimeras/${BASENAME}_combined_chimeras.fasta"
done

echo "Script completed."



###### WITHOUT DEREPLICATION ########

#!/bin/bash

# Directory containing your FASTA files
FASTA_DIR="/home/ali/Documents/simulated_data/analysis/tagJump_filt/original"

# Create folders for each step
mkdir -p chimeras
mkdir -p chimeras/combined_chimeras  # Create subfolder for combined chimeras

# Loop through each FASTA file in the directory
for FASTA_FILE in "$FASTA_DIR"/*.fasta; do
    BASENAME=$(basename "$FASTA_FILE" .fasta)

    # Chimera Filtering with UCHIME Denovo
    echo "Chimera Filtering with UCHIME Denovo for $BASENAME..."
    vsearch --uchime_denovo "$FASTA_FILE" --mindiv 0.4 --dn 1.6 --minh 0.08 --sizein --sizeout --fasta_width 0 \
    --nonchimeras "chimeras/${BASENAME}_denovo_nonchimeras.fasta" --chimeras "chimeras/${BASENAME}_denovo_chimeras.fasta" --borderline "chimeras/${BASENAME}_denovo_borderline.fasta"

    # Combining chimeric and borderline sequences for each file
    echo "Combining chimeric and borderline sequences for $BASENAME..."
    cat "chimeras/${BASENAME}_denovo_chimeras.fasta" "chimeras/${BASENAME}_denovo_borderline.fasta" > "chimeras/combined_chimeras/${BASENAME}_combined_chimeras.fasta"
done

echo "Script completed."
