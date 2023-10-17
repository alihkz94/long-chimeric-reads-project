##Libraries
library(BiocGenerics)
library(Biostrings)
library(Rcpp)
library(dada2)
library(ShortRead)
library(tibble)
setwd("/gpfs/space/home/alihakim/analysis/dada2")

# Function to read fasta sequences into a named vector
read_fasta_to_vector <- function(fasta_file) {
  fasta_seqs <- readDNAStringSet(fasta_file)
  return(setNames(as.character(fasta_seqs), names(fasta_seqs)))
}

# Read in the sequence table
seqtab <- read.csv("chimeric_sequence_table.csv", stringsAsFactors = FALSE)

# Convert to matrix format similar to your original script
seqtab_matrix <- matrix(seqtab$Frequency, nrow=nrow(seqtab), ncol=1, 
                        dimnames=list(seqtab$Sequence, "Sample1"))

# Read in the table of known chimeric sequences (from fasta)
known_chimeras <- read_fasta_to_vector("metabar_uchime_input.fasta")

# Initialize a ground truth vector
ground_truth <- ifelse(seqtab$Sequence %in% known_chimeras, 1, 0)

# List of options to test, tailored for PacBio
options_list <- list(
    list(minSampleFraction = 0.9,ignoreNNegatives=1,minFoldParentOverAbundance = 1.5, 
    minParentAbundance = 2, allowOneOff = FALSE, minOneOffParentDistance = 2, maxShift =16),
      list(minSampleFraction = 0.9,ignoreNNegatives=1,minFoldParentOverAbundance = 1.5,
      minParentAbundance = 2, allowOneOff = TRUE, minOneOffParentDistance = 2, maxShift =16)
)

# Initialize a data frame to store metrics
metrics_df <- tibble()

# Loop through options and perform chimera removal
for (opts in options_list) {
  cleaned_seqtab_matrix <- removeBimeraDenovo(t(seqtab_matrix), method = "consensus",
                                              minFoldParentOverAbundance = opts$minFoldParentOverAbundance,
                                              minSampleFraction = opts$minSampleFraction,
                                              ignoreNNegatives = opts$ignoreNNegatives,
                                              minParentAbundance = opts$minParentAbundance,
                                              allowOneOff = opts$allowOneOff,
                                              minOneOffParentDistance = opts$minOneOffParentDistance,
                                              maxShift = opts$maxShift,
                                              multithread=TRUE, verbose = TRUE)
  
  # Identify removed chimeric sequences
  removed_seqs <- rownames(seqtab_matrix)[!rownames(seqtab_matrix) %in% rownames(t(cleaned_seqtab_matrix))]
  
  # Create a binary vector where 1 means the sequence was removed
  removed_binary <- ifelse(seqtab$Sequence %in% removed_seqs, 1, 0)
  
  # Calculate TP, FP, TN, FN
  TP <- sum(removed_binary & ground_truth)
  FP <- sum(removed_binary & !ground_truth)
  TN <- sum(!removed_binary & !ground_truth)
  FN <- sum(!removed_binary & ground_truth)
  
  # Calculate Sensitivity and Specificity
  Sensitivity <- TP / (TP + FN)
  Specificity <- TN / (TN + FP)
  
  # Add to metrics data frame
  opts_str <- paste(unlist(opts), collapse = "_")
  metrics_df <- rbind(metrics_df, tibble(Options = opts_str, TP = TP, FP = FP, TN = TN, FN = FN, Sensitivity = Sensitivity, Specificity = Specificity))
}

# Save metrics data frame to a CSV file
write.csv(metrics_df, "DADA2_Sensitivity_Specificity_Metrics.csv")
