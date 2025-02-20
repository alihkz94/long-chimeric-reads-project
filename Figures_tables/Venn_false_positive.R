################################################################################
#Generate a Venn diagram to compare the results of three different methods for #
#identifying false positive chimeric sequences in a simulated dataset.         #
################################################################################

# Load necessary libraries
library(VennDiagram)
library(Biostrings)

setwd("~/Documents/simulated_data/analysis/figures/blast")

# Function to calculate intersections and generate a Venn diagram
generate_venn_diagram <- function(file1, file2, file3, output_file) {
  # Load sequences from each FASTA file into sets
  seqs1 <- unique(as.character(readDNAStringSet(file1)))
  seqs2 <- unique(as.character(readDNAStringSet(file2)))
  seqs3 <- unique(as.character(readDNAStringSet(file3)))
  
  # Calculate intersections and unique sets
  common_all <- intersect(intersect(seqs1, seqs2), seqs3)
  common_1_2 <- setdiff(intersect(seqs1, seqs2), common_all)
  common_1_3 <- setdiff(intersect(seqs1, seqs3), common_all)
  common_2_3 <- setdiff(intersect(seqs2, seqs3), common_all)
  unique_1 <- setdiff(seqs1, union(seqs2, seqs3))
  unique_2 <- setdiff(seqs2, union(seqs1, seqs3))
  unique_3 <- setdiff(seqs3, union(seqs1, seqs2))
  
  # Print the results
  cat("Unique to File 1:", length(unique_1), "\n")
  cat("Unique to File 2:", length(unique_2), "\n")
  cat("Unique to File 3:", length(unique_3), "\n")
  cat("Common to Files 1 and 2:", length(common_1_2), "\n")
  cat("Common to Files 1 and 3:", length(common_1_3), "\n")
  cat("Common to Files 2 and 3:", length(common_2_3), "\n")
  cat("Common to all three files:", length(common_all), "\n")
  
  # Generate a Venn diagram
  venn.plot <- venn.diagram(
    x = list(
      `uchime_denovo` = seqs1,
      `chimeras_denovo` = seqs2,
      `DADA2` = seqs3
    ),
    category.names = c("uchime_denovo", "chimeras_denovo", "removeBimeraDenovo"),
    filename = NULL,
    output = TRUE,
    imagetype = "png",
    height = 2000,  # Adjusted for better fit in Word document
    width = 2000,   # Adjusted for better fit in Word document
    resolution = 300,  # Sufficient for high quality without being overly large
    compression = "lzw",
    col = "black",
    fill = c("#D55E00", "#0072B2", "#CC79A7"),
    alpha = 0.50,
    cex = 1.5,  # Adjusted for better readability in Word document
    cat.cex = 1.5,  # Adjusted for better readability in Word document
    cat.pos = c(-15, 15, 180),  # Adjusted positions
    cat.dist = c(0.05, 0.05, 0.025),  # Adjusted distances
    cat.col = c("black", "black", "black"),
    main = "",
    main.cex = 0,  # No title
    margin = 0.1  # Reduced margin for better fit
  )
  
  # Save the Venn diagram in high-resolution
  jpeg(output_file, width = 2000, height = 2000, res = 300)
  grid.draw(venn.plot)
  dev.off()
}
# Example usage
uchime_denovo <- "/home/ali/Documents/simulated_data/analysis/figures/blast/uchime_denovo.fasta"
chimeras_denovo <- "/home/ali/Documents/simulated_data/analysis/figures/blast/chimeras_denovo.fasta"
DADA2<- "/home/ali/Documents/simulated_data/analysis/figures/blast/dada2.fasta"
output_file <- "/home/ali/Documents/simulated_data/analysis/figures/blast/venn_diagram.png"
generate_venn_diagram(uchime_denovo, chimeras_denovo, DADA2, output_file)
