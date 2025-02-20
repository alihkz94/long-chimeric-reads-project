# Load necessary libraries
library(VennDiagram)
library(Biostrings)

# Directory
setwd("~/Documents/simulated_data/analysis/figures/BlasCh/false_negative_chimeras_venn")

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
      `_uchime_denovo` = seqs1,
      `_chimeras_denovo` = seqs2,
      `_DADA2` = seqs3
    ),
    category.names = c(" uchime_denovo", " chimeras_denovo", " removeBimeraDenovo"),
    filename = NULL,
    output = TRUE,
    imagetype = "jpg",
    height = 2000,  # Adjusted for better fit in Word document
    width = 2000,   # Adjusted for better fit in Word document
    resolution = 300,  # Sufficient for high quality without being overly large
    compression = "lzw",
    col = "black",
    fill = c("#FF0000", "#0000FF", "#00FF00"),  # New colors
    alpha = 0.50,
    cex = 1.7,  # Adjusted for better readability in Word document
    cat.cex = 1.7,  # Adjusted for better readability in Word document
    cat.pos = c(-10, 10, 180),  # Adjusted positions to bring labels closer
    cat.dist = c(0.035, 0.035, 0.025),  # Adjusted distances to bring labels closer
    cat.just = list(c(0.5, 1), c(0.5, 1), c(0.5, 0.5)),  # Center and adjust text position
    cat.default.pos = "outer",  # Ensure the text is outside the circles
    cat.col = c("black", "black", "black"),
    main = "",
    main.cex = 0,  # No title
    margin = 0.1  # Reduced margin for better fit
  )
  
  # Save the Venn diagram in high-resolution
  jpeg(output_file, width = 3000, height = 3000, res = 300)
  grid.draw(venn.plot)
  dev.off()
}

# Example usage
file1 <- "/home/ali/Documents/simulated_data/analysis/figures/BlasCh/false_negative_chimeras_venn/uchime_denovo_nonchimeras.fasta"
file2 <- "/home/ali/Documents/simulated_data/analysis/figures/BlasCh/false_negative_chimeras_venn/chimeras_denovo_nonchimeras.fasta"
file3 <- "/home/ali/Documents/simulated_data/analysis/figures/BlasCh/false_negative_chimeras_venn/DADA2_nonchimeras.fasta"
output_file <- "/home/ali/Documents/simulated_data/analysis/figures/BlasCh/false_negative_chimeras_venn/venn_diagram.jpg"
generate_venn_diagram(file1, file2, file3, output_file)


###########SVG format#########

# Load necessary libraries
library(VennDiagram)
library(Biostrings)

# Directory
setwd("~/Documents/simulated_data/analysis/figures/BlasCh/false_negative_chimeras_venn")

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
      `_uchime_denovo` = seqs1,
      `_chimeras_denovo` = seqs2,
      `_DADA2` = seqs3
    ),
    category.names = c(" UCHIME_denovo", " chimeras_denovo", " DADA2"),
    filename = NULL,
    output = TRUE,
    imagetype = "svg",
    height = 3000,  # Adjusted for better fit in Word document
    width = 3000,   # Adjusted for better fit in Word document
    resolution = 300,  # Sufficient for high quality without being overly large
    compression = "lzw",
    col = "black",
    fill = c("#FF0000", "#0000FF", "#00FF00"),  # New colors
    alpha = 0.50,
    cex = 2,  # Increased for better readability of intersection counts
    cat.cex = 2,  # Increased for better readability of category labels
    cat.pos = c(-10, 10, 180),  # Adjusted positions to bring labels closer
    cat.dist = c(0.035, 0.035, 0.025),  # Adjusted distances to bring labels closer
    cat.just = list(c(0.5, 1), c(0.5, 1), c(0.5, 0.5)),  # Center and adjust text position
    cat.default.pos = "outer",  # Ensure the text is outside the circles
    cat.col = c("black", "black", "black"),
    main = "",
    main.cex = 0,  # No title
    margin = 0.1  # Reduced margin for better fit
  )
  
  # Save the Venn diagram in SVG format
  svg(output_file, width = 10, height = 10)  # Dimensions in inches
  grid.draw(venn.plot)
  dev.off()
}

# Example usage
file1 <- "/home/ali/Documents/simulated_data/analysis/figures/BlasCh/false_negative_chimeras_venn/uchime_denovo_nonchimeras.fasta"
file2 <- "/home/ali/Documents/simulated_data/analysis/figures/BlasCh/false_negative_chimeras_venn/chimeras_denovo_nonchimeras.fasta"
file3 <- "/home/ali/Documents/simulated_data/analysis/figures/BlasCh/false_negative_chimeras_venn/DADA2_nonchimeras.fasta"
output_file <- "/home/ali/Documents/simulated_data/analysis/figures/BlasCh/false_negative_chimeras_venn/venn_diagram.svg"
generate_venn_diagram(file1, file2, file3, output_file)
