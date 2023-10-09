##Libraries
library(BiocGenerics)
library(Biostrings)
library(Rcpp)
library(dada2)
library(ShortRead)
library(tibble)
setwd("~/Documents/simulated_data/analysis/DADA2")

# Read in sequence table
seqtab <- read.csv("chimeric_sequence_table.csv", stringsAsFactors = FALSE)

# Convert to matrix format
seqtab_matrix <- matrix(seqtab$Frequency, nrow=nrow(seqtab), ncol=1, 
                        dimnames=list(seqtab$Sequence, "Sample1"))

# List of options to test, tailored for PacBio
options_list <- list(
  list(minFoldParentOverAbundance = 1.5, minParentAbundance = 5, allowOneOff = TRUE, minOneOffParentDistance = 6, maxShift = 50),
  list(minFoldParentOverAbundance = 2, minParentAbundance = 10, allowOneOff = TRUE, minOneOffParentDistance = 6, maxShift = 100),
  list(minFoldParentOverAbundance = 1, minParentAbundance = 8, allowOneOff = TRUE, minOneOffParentDistance = 4, maxShift = 16)
)

# Initialize a data frame to store metrics
metrics_df <- tibble()

# Loop through options and perform chimera detection
for (opts in options_list) {
  chimera <- isBimeraDenovo(t(seqtab_matrix),
                            minFoldParentOverAbundance = opts$minFoldParentOverAbundance,
                            minParentAbundance = opts$minParentAbundance,
                            allowOneOff = opts$allowOneOff,
                            minOneOffParentDistance = opts$minOneOffParentDistance,
                            maxShift = opts$maxShift,
                            multithread=TRUE, verbose = TRUE)
  
  # Extract chimeric sequences
  Chimeric_sequences <- seqtab[chimera,]
  
  # Save results for further analysis
  output_file <- paste0("Chimeric_sequences_", paste(unlist(opts), collapse = "_"), ".csv")
  write.csv(Chimeric_sequences, output_file)
  
  # Assume you have a ground truth vector 'ground_truth' where 1 means chimera and 0 means non-chimera
  # ground_truth <- c(1, 0, 1, 0, ...)
  
  # Calculate TP, FP, TN, FN
  TP <- sum(chimera & ground_truth)
  FP <- sum(chimera & !ground_truth)
  TN <- sum(!chimera & !ground_truth)
  FN <- sum(!chimera & ground_truth)
  
  # Calculate Sensitivity and Specificity
  Sensitivity <- TP / (TP + FN)
  Specificity <- TN / (TN + FP)
  
  # Add to metrics data frame
  opts_str <- paste(unlist(opts), collapse = "_")
  metrics_df <- rbind(metrics_df, tibble(Options = opts_str, TP = TP, FP = FP, TN = TN, FN = FN, Sensitivity = Sensitivity, Specificity = Specificity))
}

# Save metrics data frame to a CSV file
write.csv(metrics_df, "DADA2_Sensitivity_Specificity_Metrics.csv")
