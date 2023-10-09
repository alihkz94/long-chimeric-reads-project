#!/bin/bash
#SBATCH --job-name="blast_treat"
#SBATCH --cpus-per-task=128
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --partition amd
#SBATCH --time=10:00:00
#SBATCH --output=blast_treat.out

set -e # Exit on error

# Input: FASTA file for chimera filtering
INPUT_FILE="/gpfs/space/home/alihakim/analysis/treat_chimera/metabar_uchime_input.fasta"
REFERENCE_DB="/gpfs/space/home/alihakim/analysis/treat_chimera/ITS.fasta"

# 1. Pre-dereplication
echo "Pre-dereplication..."
vsearch --derep_fulllength $INPUT_FILE --sizein --sizeout --fasta_width 0 --output dereplicated.fasta

# 2. Pre-clustering
echo "Pre-clustering..."
vsearch --cluster_size dereplicated.fasta --id 0.97 --sizein --sizeout --fasta_width 0 --centroids preclustered.fasta

# 3. Chimera Filtering with UCHIME Denovo
echo "Chimera Filtering with UCHIME Denovo..."
vsearch --uchime_denovo preclustered.fasta --mindiv 0.4 --dn 1.6 --minh 0.01 --threads 128 --uchimeout denovo.uchime \
--nonchimeras denovo_nonchimeras.fasta --chimeras denovo_chimeras.fasta

# 4. Chimera Filtering with UCHIME Ref
echo "Chimera Filtering with UCHIME Ref..."
vsearch --uchime_ref denovo_nonchimeras.fasta --mindiv 0.4 --dn 1.6 --minh 0.01 --threads 128 \
--db $REFERENCE_DB \
--nonchimeras ref_nonchimeras.fasta --chimeras ref_chimeras.fasta

# 5. Combine chimeric sequences
echo "Combining chimeric sequences..."
cat denovo_chimeras.fasta ref_chimeras.fasta > combined_chimeras.fasta

# 6. Reverse complement all sequences
echo "Reverse complementing sequences..."
seqkit seq -r combined_chimeras.fasta > combined_chimeras_reversed.fasta

# 7. BLAST Analysis for Flagged Chimeric Sequences (Original and Reversed)
echo "BLAST analysis for flagged chimeric sequences..."
blastn -query combined_chimeras.fasta \
-db /gpfs/space/home/alihakim/analysis/databse/UNITE  \
-word_size 7 \
-task blastn \
-num_threads 128 \
-outfmt "6 delim=+ qseqid stitle qlen slen qstart qend sstart send evalue length nident mismatch gapopen gaps sstrand qcovs pident" \
-evalue 0.001 \
-strand both \
-max_target_seqs 10 \
-max_hsps 1 \
-out combined_chimeras_blast_results.txt

blastn -query combined_chimeras_reversed.fasta \
-db /gpfs/space/home/alihakim/analysis/databse/UNITE  \
-word_size 7 \
-task blastn \
-num_threads 128 \
-outfmt "6 delim=+ qseqid stitle qlen slen qstart qend sstart send evalue length nident mismatch gapopen gaps sstrand qcovs pident" \
-evalue 0.001 \
-strand both \
-max_target_seqs 10 \
-max_hsps 1 \
-out combined_chimeras_reversed_blast_results.txt

# 8. Additional Filtering Steps
PIDENT_THRESHOLD=90
QC_THRESHOLD=90

awk -F'+' -v pident="$PIDENT_THRESHOLD" -v qc="$QC_THRESHOLD" '{ if ($15 >= pident && $14 >= qc) print $1 }' combined_chimeras_blast_results.txt > true_non_chimeras.txt
awk -F'+' -v pident="$PIDENT_THRESHOLD" -v qc="$QC_THRESHOLD" '{ if ($15 >= pident && $14 >= qc) print $1 }' combined_chimeras_reversed_blast_results.txt > true_non_chimeras_reversed.txt

# 9. Extract true non-chimeric sequences
seqkit grep -f true_non_chimeras.txt combined_chimeras.fasta > retrieved_non_chimeras.fasta
seqkit grep -f true_non_chimeras_reversed.txt combined_chimeras_reversed.fasta > retrieved_non_chimeras_reversed.fasta

# 10. Combine the retrieved non-chimeric sequences
cat ref_nonchimeras.fasta retrieved_non_chimeras.fasta retrieved_non_chimeras_reversed.fasta > final_non_chimeric_sequences.fasta

echo "Chimera filtering module completed!"
