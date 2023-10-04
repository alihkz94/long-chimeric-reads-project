##Libraries
library(BiocGenerics)
library(Biostrings)
library(Rcpp)
library(dada2)
library(ShortRead)
setwd("~/dada2")

####chimera filtering####
seqtab <- read.csv("chimeric_sequence_table.csv", stringsAsFactors = FALSE)
seqtab_matrix <- matrix(seqtab$Frequency, nrow=nrow(seqtab), ncol=1, 
                        dimnames=list(seqtab$Sequence, "Sample1"))


result <- removeBimeraDenovo(t(seqtab_matrix), method="consensus", 
                                    multithread=TRUE, verbose = TRUE)


seqtab_matrix_transposed <- t(seqtab_matrix)
# Get the names of the sequences before chimera removal
seqs_before <- colnames(seqtab_matrix_transposed)

# Get the names of the sequences after chimera removal
seqs_after <- colnames(result)

# Find chimeric sequences (those present before but not after chimera removal)
chimeric_seqs <- setdiff(seqs_before, seqs_after)
# Concatenate chimeric sequences into a single string with spaces between them
chimeric_seqs_string <- paste(chimeric_seqs, collapse = " ")

# Save chimeric sequences to a file
writeLines(chimeric_seqs, con = "chimeric_sequences.txt")
