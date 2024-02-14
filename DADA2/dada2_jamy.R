# Libraries
library(BiocGenerics)
library(Biostrings)
library(Rcpp)
library(dada2)
library(ShortRead)
library(tibble)
library(seqinr)
library(Hmisc)
library(openssl)

# Set working directory
setwd("~/dada2")

# Load pre-processed FASTA files (assumed dereplicated)
fasta_files <- list.files(pattern="*.fasta", full.names=TRUE)

# Initialize an empty list for sequence tables
seqtabs <- list()

for(fasta_file in fasta_files) {
  seqs <- readDNAStringSet(fasta_file)
  seqtabs[[basename(fasta_file)]] <- table(sapply(seqs, as.character))
}

# Convert the list of tables into a single data frame
seqtab <- do.call(rbind, seqtabs)

# Ensure column names are unique sequences
colnames(seqtab) <- make.names(colnames(seqtab), unique=TRUE)

# Optionally, save sequence table for later use
saveRDS(seqtab, "sequence_table.rds")
# Load sequence table if needed
seqtab <- readRDS("sequence_table.rds")

# Chimera removal (if not already done)
nonchimeric_seqtab <- removeBimeraDenovo(seqtab, method = "consensus",
                                         minFoldParentOverAbundance = 1.5, minParentAbundance = 2, 
                                         allowOneOff = TRUE, minOneOffParentDistance = 2, maxShift = 100, 
                                         verbose = TRUE)


# Number of sequences before and after
orig_count <- ncol(seqtab)  
filtered_count <- ncol(nonchimeric_seqtab)
discarded_count <- orig_count - filtered_count

# Extract chimeric sequences
discarded_chimera <- seqtab[,colnames(seqtab) %nin% colnames(nonchimeric_seqtab)]
chimeric_sequences <- getSequences(discarded_chimera)

# Create a DNAStringSet object
chimeric_dna <- DNAStringSet(chimeric_sequences)

# Extract chimeric sequences - already done in your script
chimeric_sequences <- colnames(discarded_chimera)

# Create headers for the FASTA file
chimera_headers <- paste0(">", "chimera_", seq_along(chimeric_sequences))

# Combine headers and sequences into a FASTA format
fasta_lines <- mapply(function(header, sequence) paste(header, sequence, sep = "\n"), chimera_headers, chimeric_sequences, SIMPLIFY = FALSE)
fasta_content <- unlist(fasta_lines)

# Write the FASTA file
writeLines(fasta_content, con = "chimeric_sequences.fasta")
