# Quality filtering with DADA2 pipeline for real data
# Libraries
library(Rcpp)
library(dada2)

# Directory
setwd("/media/ali/data_store/test_dada2")

# List all fastq files in the working directory
fnFs <- list.files(pattern = ".fastq", full.names = TRUE)

# Read the original FASTQ files
original_reads <- lapply(fnFs, readFastq)

# Convert to DNAStringSet objects
original_dna <- lapply(original_reads, sread)

# Extract the headers
original_headers <- lapply(original_dna, names)

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