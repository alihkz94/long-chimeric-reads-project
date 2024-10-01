#!/bin/bash
#SBATCH --job-name="xml_dada2"
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --partition amd
#SBATCH --time=180:00:00

# Directory containing the FASTA files and where results will be stored
DIR="."

# Directory containing the BLAST databases we created
DB_DIR="$DIR/blast_dbs"

# Number of CPUs to use for BLAST
CPUS=32

# Path to the existing BLAST database
BLAST_DB="/gpfs/space/home/alihakim/database/euk/EUK"

# Run BLAST analysis only for FASTA files
for chimera_file in $DIR/*.fasta; do
    # Get the base name of the file (without .fasta extension)
    BASENAME=$(basename "$chimera_file" .fasta)
    
    # Path to the corresponding BLAST database we created
    ERR_DB="${DB_DIR}/${BASENAME}/${BASENAME}"
    
    # Run BLASTn
    blastn -query "$chimera_file" \
           -db "$BLAST_DB $ERR_DB" \
           -word_size 7 \
           -task blastn \
           -num_threads "$CPUS" \
           -outfmt 5 \
           -evalue 0.001 \
           -strand both \
           -max_target_seqs 10 \
           -max_hsps 1 \
           -out "$DIR/${BASENAME}_blast_results.xml"
    
    echo "Completed BLAST analysis for ${BASENAME}"
done

echo "All BLAST analyses completed."



###### script for the HPC to run seeveral times in the meantime in arrays#####
#!/bin/bash

#SBATCH --job-name=3dada
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --partition=amd
#SBATCH --time=180:00:00
#SBATCH --array=1-5%5  # Adjust this based on the number of FASTA files

# Directory containing the FASTA files and where results will be stored
DIR="."
# Directory containing the BLAST databases we created
DB_DIR="$DIR/blast_dbs"
# Number of CPUs to use for BLAST
CPUS=32
# Path to the existing BLAST database
BLAST_DB="/gpfs/space/home/alihakim/database/euk/EUK"

# Get the list of FASTA files
FASTA_FILES=($(ls $DIR/*.fasta))

# Get the current FASTA file based on the array task ID
CURRENT_FILE=${FASTA_FILES[$SLURM_ARRAY_TASK_ID-1]}

# Get the base name of the file (without .fasta extension)
BASENAME=$(basename "$CURRENT_FILE" .fasta)

# Path to the corresponding BLAST database we created
ERR_DB="${DB_DIR}/${BASENAME}/${BASENAME}"

# Run BLASTn
blastn -query "$CURRENT_FILE" \
       -db "$BLAST_DB $ERR_DB" \
       -word_size 7 \
       -task blastn \
       -num_threads "$CPUS" \
       -outfmt 5 \
       -evalue 0.001 \
       -strand both \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -out "$DIR/${BASENAME}_blast_results.xml"

echo "Completed BLAST analysis for ${BASENAME}"
