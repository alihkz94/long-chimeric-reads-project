#!/bin/bash
set -e # Exit on error

# Input: FASTA file for chimera filtering
INPUT_FILE="input_sequences.fasta"
REFERENCE_DB="path_to_reference_db.fasta"  # Modify this to the path of your reference database for BLAST

# 1. Pre-dereplication
echo "Pre-dereplication..."
vsearch --derep_fulllength $INPUT_FILE --sizein --sizeout --fasta_width 0 --output dereplicated.fasta

# 2. Pre-clustering
echo "Pre-clustering..."
vsearch --cluster_size dereplicated.fasta --id 0.97 --sizein --sizeout --fasta_width 0 --centroids preclustered.fasta

# 3. Chimera Filtering with UCHIME Denovo
echo "Chimera Filtering with UCHIME Denovo..."
vsearch --uchime_denovo preclustered.fasta --mindiv 0.5 --dn 1.6 --threads 8 --uchimeout denovo.uchime \
--nonchimeras denovo_nonchimeras.fasta --chimeras denovo_chimeras.fasta

# 4. Chimera Filtering with UCHIME Ref
echo "Chimera Filtering with UCHIME Ref..."
vsearch --uchime_ref denovo_nonchimeras.fasta --mindiv 0.5 --dn 1.6 --threads 8 \
--db $REFERENCE_DB \
--nonchimeras ref_nonchimeras.fasta --chimeras ref_chimeras.fasta

# 5. BLAST Analysis for Flagged Chimeric Sequences
echo "BLAST analysis for flagged chimeric sequences..."

# Combine chimeric sequences from both UCHIME Denovo and UCHIME Ref
cat denovo_chimeras.fasta ref_chimeras.fasta > combined_chimeras.fasta

blastn -query combined_chimeras.fasta \
-db $REFERENCE_DB \
-task blastn \
-word_size 7 \
-num_threads 12 \
-outfmt "6 delim=+ qseqid stitle qlen slen qstart qend sstart send evalue length nident mismatch gapopen gaps sstrand qcovs pident" \
-evalue 0.001 \
-strand both \
-max_target_seqs 10 \
-max_hsps 1 \
-out combined_chimeras_blast_results.txt

# Define a threshold for percentage identity. Sequences above this threshold are considered true non-chimeric.
PIDENT_THRESHOLD=90  # Adjust based on your dataset and requirements.

awk -F'+' -v threshold="$PIDENT_THRESHOLD" '{ if ($15 >= threshold) print $1 }' combined_chimeras_blast_results.txt > true_non_chimeras.txt

# Extract true non-chimeric sequences from the flagged chimeric ones
seqkit grep -f true_non_chimeras.txt combined_chimeras.fasta > retrieved_non_chimeras.fasta

# Combine the retrieved non-chimeric sequences with the initial non-chimeric sequences
cat ref_nonchimeras.fasta retrieved_non_chimeras.fasta > final_non_chimeric_sequences.fasta

echo "Chimera filtering module completed!"
