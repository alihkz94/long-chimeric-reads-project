####Libraries####
library(BiocGenerics)
library(Biostrings)
library(Rcpp)
library(dada2)
library(ShortRead)
library(tibble)

####Directory####
setwd("/gpfs/space/home/alihakim/vlad")


# List all fastq files in the working directory
fnFs <- list.files(pattern = ".fastq", full.names = TRUE)

# Define file paths for filtered files
filtFs <- file.path("filtered", basename(fnFs))

# Quality filtering
out <- filterAndTrim(fnFs, filtFs, truncLen = 0, maxN = 0, maxEE = 2, truncQ = 2, minQ = 3, rm.phix = FALSE, multithread = TRUE)
cat("Filtered", sum(out), "reads from", length(out), "sample(s).\n")


# Check if all filtered files exist
missing_files <- filtFs[!file.exists(filtFs)]
if(length(missing_files) > 0) {
  warning("The following filtered files do not exist:\n", paste(missing_files, collapse = "\n"))
}

# Update filtFs to only include existing files
filtFs <- filtFs[file.exists(filtFs)]

# Determine the number of subsets to process
n_subsets <- 2  # Adjust this number based on your system's capacity

# Determine the number of files per subset
files_per_subset <- ceiling(length(filtFs) / n_subsets)

# Initialize an empty list to store the dereplicated data
derep_list <- vector("list", n_subsets)

# Loop through each subset of files
for(i in seq_len(n_subsets)){
  
  # Determine the file indices for this subset
  idx <- ((i - 1) * files_per_subset + 1):min(i * files_per_subset, length(filtFs))
  
  # Select the files for this subset
  subset_files <- filtFs[idx]
  
  # Dereplicate this subset of files
  derep_list[[i]] <- derepFastq(subset_files, verbose = TRUE)
  
}

# Combine all the dereplicated data into a single list
derepFs <- do.call(c, derep_list)

# Learn error rates
dadaFs <- dada(derepFs, err = NULL, selfConsist = TRUE, HOMOPOLYMER_GAP_PENALTY = 1, multithread = TRUE)

# Create sequence table
seqtab <- makeSequenceTable(dadaFs)
# Inspect sequence table
print(dim(seqtab))
# Save sequence table for later chimera filtering
saveRDS(seqtab, "sequence_table.rds")


# Load the sequence table if needed
seqtab <- readRDS("sequence_table.rds")

# Chimera filtering
nonchimeric_seqtab <- removeBimeraDenovo(seqtab, minSampleFraction = 0.9, ignoreNNegatives = 1, 
                                         minFoldParentOverAbundance = 1.5, minParentAbundance = 2, 
                                         allowOneOff = TRUE, minOneOffParentDistance = 2, maxShift = 300, 
                                         multithread = TRUE, verbose = TRUE)

# Identify chimeric sequences
chimeric_seqtab <- seqtab[!rownames(seqtab) %in% rownames(nonchimeric_seqtab),]

# Saving the results to fasta files
uniquesToFasta(nonchimeric_seqtab, fout="nonchimeric_sequences.fasta")
uniquesToFasta(chimeric_seqtab, fout="chimeric_sequences.fasta")
