#!/bin/bash
#SBATCH --job-name="blast"
#SBATCH --cpus-per-task=128
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --partition amd
#SBATCH --time=10:00:00
#SBATCH --output=blast_new.out

 Directory containing your fastq files
FASTQ_DIR="/gpfs/space/home/alihakim/vlad/filtered"

# Convert FastQ to FastA
echo "Converting FastQ to FastA..."
mkdir -p fasta_files
for file in $FASTQ_DIR/*.fastq; do
    seqkit fq2fa "$file" > "fasta_files/$(basename "${file%.fastq}.fasta")"
done

# Directory containing your fasta files
FASTA_DIR="fasta_files"

# Process each fasta file
for file in $FASTA_DIR/*.fasta; do

    BASENAME=$(basename "$file" .fasta)

    # Create folders for each step
    mkdir -p dereplicated preclustered chimeras

    # 1. Pre-dereplication
    echo "Pre-dereplication for $BASENAME..."
    vsearch --derep_fulllength "$file" --sizein --sizeout --fasta_width 0 --output "dereplicated/${BASENAME}_dereplicated.fasta"

done

# Process dereplicated files
for file in dereplicated/*.fasta; do

    BASENAME=$(basename "$file" .fasta)

    # 2. Pre-clustering
    echo "Pre-clustering for $BASENAME..."
    vsearch --cluster_size "$file" --id 0.97 --sizein --sizeout --fasta_width 0 --centroids "preclustered/${BASENAME}_preclustered.fasta"

    # 3. Chimera Filtering with UCHIME Denovo
    echo "Chimera Filtering with UCHIME Denovo for $BASENAME..."
    vsearch --uchime_denovo "preclustered/${BASENAME}_preclustered.fasta" --mindiv 0.4 --dn 1.6 --minh 0.06 --threads 128 --uchimeout "chimeras/${BASENAME}_denovo.uchime" \
    --nonchimeras "chimeras/${BASENAME}_denovo_nonchimeras.fasta" --chimeras "chimeras/${BASENAME}_denovo_chimeras.fasta" --borderline "chimeras/${BASENAME}_denovo_borderline.fasta"

done

# Combine chimeric and borderline sequences from all files
echo "Combining chimeric and borderline sequences..."
cat chimeras/*_denovo_chimeras.fasta chimeras/*_denovo_borderline.fasta > combined_chimeras.fasta

# BLAST Analysis for Flagged Chimeric Sequences (Original and Reversed)

# Create folder for BLAST output
mkdir -p blast_out

# Splitting the combined chimeras file into chunks
TOTAL_SEQS=$(grep -c '^>' combined_chimeras.fasta)
NUM_CHUNKS=10
SEQ_PER_CHUNK=$((TOTAL_SEQS / NUM_CHUNKS))
awk -v prefix="chunk_" -v chunk_size="$SEQ_PER_CHUNK" '/^>/{n++;if(n%chunk_size==1) {m++; close(f); f=prefix m ".fasta"}} {print > f}' combined_chimeras.fasta

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
cat blast_out/*_blast_results.txt > blast_out/combined_chimeras_blast_results.txt
rm blast_out/*_blast_results.txt
