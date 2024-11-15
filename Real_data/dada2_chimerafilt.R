## DADA2 Pipeline for Chimera Filtering and Sequence Processing ##

# load necessary libraries
library(Biostrings)
library(data.table)
library(dplyr)
library(tidyr)
library(dada2)
library(digest)
library(fs)
library(Hmisc)
library(tidyverse)
# Setting directory and reading sequences
setwd("/gpfs/space/home/alihakim/mypc/dada2_all")
sequences <- readDNAStringSet("replicated.fasta")

# Extract sample names from sequence headers
# Assuming the sample name is the part before the first dot in the FASTA header
sample_ids <- sapply(names(sequences), function(x) {
  s <- strsplit(x, "\\.")[[1]]
  return(s[1])
})

# Create a data frame of sequences and sample IDs
seq_data <- data.frame(sample_id = sample_ids, sequence = as.character(sequences), stringsAsFactors = FALSE)

# Display some of the sample names to confirm correct extraction
print("Sample IDs extracted (sample):")
print(unique(sample_ids))
# Count sequence occurrences per sample
seq_counts <- seq_data %>%
  group_by(sample_id, sequence) %>%
  summarise(count = n(), .groups = 'drop')

# Pivot the data to create a matrix-like table with samples as rows and sequences as columns
seqtab_df <- pivot_wider(seq_counts, names_from = sequence, values_from = count, values_fill = list(count = 0))

# Convert to a matrix and set row names as sample IDs for DADA2
seqtab <- as.matrix(seqtab_df[,-1])  # Excluding sample_id column for the matrix
rownames(seqtab) <- seqtab_df$sample_id

# Print sample matrix data to verify structure
print("Dimension of sequence table:")
print(dim(seqtab))

# Save the sequence table for later use
saveRDS(seqtab, "sequence_table.rds")

# Load and use in DADA2
seqtab <- readRDS("sequence_table.rds")
# Existing code for chimera filtering
nonchimeric_seqtab <- removeBimeraDenovo(seqtab,method = "consensus",
                                         multithread = TRUE, verbose = TRUE)
# Extract sequences
nonchimeric_sequences <- as.character(getSequences(nonchimeric_seqtab))
discarded_chimera <- seqtab[,colnames(seqtab) %nin% colnames(nonchimeric_seqtab)]
chimeric_sequences <- as.character(getSequences(discarded_chimera))

# Ensure directory existence
ensureDirExists <- function(dir_name) {
  if (!dir_exists(dir_name)) {
    dir_create(dir_name)
  }
}

# Write sequences to FASTA with counts, excluding sequences with zero counts
writeFastaWithCounts <- function(sequences, counts, headers, dir_name, file_name) {
  ensureDirExists(dir_name)
  file_path <- file.path(dir_name, file_name)
  
  # Initialize a vector to hold FASTA content
  fasta_lines <- c()
  
  # Loop through sequences and prepare FASTA content
  for (i in seq_along(sequences)) {
    if (counts[i] > 0) {  # Check if count is greater than zero
      # Modified header format to use ";size=" 
      fasta_header <- paste(">", headers[i], ";size=", counts[i], sep = "")
      fasta_sequence <- sequences[i]
      fasta_lines <- c(fasta_lines, fasta_header, fasta_sequence)
    }
  }
  
  # Write to file only if there's content to avoid empty files
  if (length(fasta_lines) > 0) {
    con <- file(file_path, "w")
    writeLines(fasta_lines, con)
    close(con)
  }
}

# Process and organize sequences by sample, respecting counts
processAndWriteSequences <- function(seqtab, isChimeric = FALSE) {
  # Regular expression to match specified file extensions
  pattern <- "\\.fastq\\.gz$|\\.fastq$|\\.fq\\.gz$|\\.fq$|\\.fa\\.gz$|\\.fasta$"
  
  # Modify sample names by removing the matched file extensions
  samples <- unique(gsub(pattern, "", rownames(seqtab)))
  dir_name <- ifelse(isChimeric, "chimeras", "non_chimeric")
  
  # Loop through each sample
  for (sample in samples) {
    # Subset for the current sample; ensure drop = FALSE to keep data frame structure
    sample_seqtab <- seqtab[grep(paste0("^", sample), rownames(seqtab)), , drop = FALSE]
    
    # Skip processing for empty or all-zero-count samples
    if (ncol(sample_seqtab) == 0 || all(colSums(sample_seqtab) == 0)) next
    
    # Extracting sequences and counts
    sequences <- colnames(sample_seqtab)
    counts <- colSums(sample_seqtab)
    valid_indices <- counts > 0  # Indices where count > 0
    headers <- sapply(sequences[valid_indices], function(x) digest(x, algo = "sha1"))
    
    # File naming convention
    file_name <- paste0(sample, ifelse(isChimeric, "_chimeric.fasta", ".fasta"))
    
    # Write sequences with valid counts to FASTA
    writeFastaWithCounts(sequences[valid_indices], counts[valid_indices], headers, dir_name, file_name)
  }
}

processAndWriteSequences(nonchimeric_seqtab, isChimeric = FALSE)
processAndWriteSequences(discarded_chimera, isChimeric = TRUE)
