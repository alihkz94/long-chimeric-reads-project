#!/bin/bash
#SBATCH --job-name="d_blast"
#SBATCH --cpus-per-task=128
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --partition amd
#SBATCH --time=20:00:00

mkdir -p blast_out

# Splitting the chimeric_sequences file into chunks
TOTAL_SEQS=$(grep -c '^>' chimeric_sequences.fasta)
NUM_CHUNKS=10
SEQ_PER_CHUNK=$((TOTAL_SEQS / NUM_CHUNKS))
awk -v prefix="chunk_" -v chunk_size="$SEQ_PER_CHUNK" '/^>/{n++;if(n%chunk_size==1) {m++; close(f); f=prefix m ".fasta"}} {print > f}' chimeric_sequences.fasta

# Running BLAST on each chunk in parallel
for file in chunk_*.fasta; do
    (
        blastn -query $file \
        -db /gpfs/space/home/alihakim/analysis/database/UNITE  \
        -word_size 7 \
        -task blastn \
        -num_threads 12 \
        -outfmt "6 delim=+ qseqid stitle qlen slen qstart qend sstart send evalue length nident mismatch gapopen gaps sstrand qcovs pident" \
        -evalue 0.001 \
        -strand both \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -out blast_out/${file%.fasta}_blast_results.txt

        rm $file
    ) &
done
wait

# Concatenating the BLAST results into one file
cat blast_out/*_blast_results.txt > blast_out/chimeric_sequences_blast_results.txt
rm blast_out/*_blast_results.txt
