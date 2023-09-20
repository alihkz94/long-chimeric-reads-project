# After UCHIME Ref chimera filtering

# 6.5. BLAST Analysis for Flagged Chimeric Sequences
echo "BLAST analysis for flagged chimeric sequences..."

# Combine chimeric sequences from both UCHIME Denovo and UCHIME Ref
cat uchime_denovo_out/*.uchime.chimeras uchime_ref_out/*.uchime.chimeras > combined_chimeras.fasta

# BLAST the combined chimeric sequences
blastn -query combined_chimeras.fasta \
-db /home/ali/Documents/simulated_data/analysis/pipe_test/database/UNITE \
-word_size 7 \
-num_threads 8 \
-outfmt "6 delim=+ qseqid stitle qlen slen qstart qend sstart send evalue length nident mismatch gapopen gaps sstrand qcovs pident" \
-evalue 0.001 \
-strand both \
-max_target_seqs 10 \
-max_hsps 2 \
-dust no \
-soft_masking true \
-penalty -1 \
-reward 1 \
-gapopen 1 \
-gapextend 2 \
-out combined_chimeras_blast_results.txt

# Define a threshold for percentage identity. Sequences above this threshold are considered true non-chimeric.
PIDENT_THRESHOLD=90  # Adjust based on your dataset and requirements.

awk -F'+' -v threshold="$PIDENT_THRESHOLD" '{ if ($15 >= threshold) print $1 }' combined_chimeras_blast_results.txt > true_non_chimeras.txt

# Extract true non-chimeric sequences from the flagged chimeric ones
seqkit grep -f true_non_chimeras.txt combined_chimeras.fasta > retrieved_non_chimeras.fasta

# Combine the retrieved non-chimeric sequences with the initial non-chimeric sequences
cat uchime_ref_out/*.fasta retrieved_non_chimeras.fasta > final_non_chimeric_sequences.fasta
