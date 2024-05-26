#!/bin/bash
#SBATCH --job-name="chimeras"
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --partition amd
#SBATCH --time=02:00:00

DIR="/gpfs/space/home/alihakim/homepc/usearch/chimeras_denovo_default"
DB="/gpfs/space/home/alihakim/homepc/usearch/utax.fasta"
REPORT_FILE="$DIR/vsearch_summary_report.txt"

# Initialize report file
echo -e "file_name\t matching_sequences\t total_sequences\t percentage" > $REPORT_FILE

# Process each FASTA file in the directory
for fasta_file in "$DIR"/*.fasta; do
    # Extract the basename of the fasta file
    BASENAME=$(basename "$fasta_file" .fasta)
    
    # Prepare output directory and files
    OUTPUT_DIR="$DIR/$BASENAME"
    mkdir -p "$OUTPUT_DIR"
    MATCHED_OUTPUT="$OUTPUT_DIR/${BASENAME}_matched_output.txt"
    UNMATCHED_OUTPUT="$OUTPUT_DIR/${BASENAME}_unmatched_output.fasta"
    VSEARCH_LOG="$OUTPUT_DIR/${BASENAME}_vsearch.log"

    # Run vsearch --usearch_global and capture the log output
    vsearch --usearch_global $fasta_file \
            --db $DB \
            --id 0.8 \
            --blast6out $MATCHED_OUTPUT \
            --strand both \
            --threads 32 \
            --maxaccepts 1 \
            --maxrejects 32 \
            --query_cov 0.8 \
            --output_no_hits > $VSEARCH_LOG 2>&1

    # Extract sequences that did not match (query coverage < 0.8)
    awk '$2 == "*" {print $1}' $MATCHED_OUTPUT > $OUTPUT_DIR/unmatched_queries.txt

    # Create a file with unmatched sequences
    seqtk subseq $fasta_file $OUTPUT_DIR/unmatched_queries.txt > $UNMATCHED_OUTPUT

    # Extract matching sequences statistics from the log
    MATCHED_SEQUENCES=$(grep -oP 'Matching unique query sequences: \K\d+' $VSEARCH_LOG)
    TOTAL_SEQUENCES=$(grep -oP 'Matching unique query sequences: \d+ of \K\d+' $VSEARCH_LOG)
    PERCENTAGE=$(grep -oP 'Matching unique query sequences: \d+ of \d+ \(\K[\d.]+(?=%)' $VSEARCH_LOG)

    # Append the statistics to the report file
    echo -e "$BASENAME.fasta\t$MATCHED_SEQUENCES\t$TOTAL_SEQUENCES\t$PERCENTAGE" >> $REPORT_FILE
done