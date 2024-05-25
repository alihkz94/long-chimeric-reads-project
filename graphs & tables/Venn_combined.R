# Load necessary libraries
library(VennDiagram)
library(Biostrings)
library(grid)
library(gridExtra)

# Function to calculate intersections and generate a Venn diagram
generate_venn_diagram <- function(file1, file2, file3, labels, colors, title) {
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
  cat(title, "\n")
  cat("Unique to File 1:", length(unique_1), "\n")
  cat("Unique to File 2:", length(unique_2), "\n")
  cat("Unique to File 3:", length(unique_3), "\n")
  cat("Common to Files 1 and 2:", length(common_1_2), "\n")
  cat("Common to Files 1 and 3:", length(common_1_3), "\n")
  cat("Common to Files 2 and 3:", length(common_2_3), "\n")
  cat("Common to all three files:", length(common_all), "\n")
  
  # Generate a Venn diagram
  venn.plot <- venn.diagram(
    x = setNames(list(seqs1, seqs2, seqs3), labels),
    category.names = labels,
    filename = NULL,
    output = FALSE,
    imagetype = "png",
    height = 2000,
    width = 2000,
    resolution = 300,
    compression = "lzw",
    col = "black",
    fill = colors,
    alpha = 0.50,
    cex = 2.5,  # Increase text size for intersection counts
    cat.cex = 2.5,  # Increase text size for category labels
    cat.pos = c(-15, 15, 180),
    cat.dist = c(0.05, 0.05, 0.025),
    cat.col = c("black", "black", "black"),
    main = title,
    main.cex = 0,
    margin = 0.1
  )
  return(gTree(children = venn.plot))
}

# File paths
uchime_denovo <- "/home/ali/Documents/simulated_data/analysis/Venn2/best/uchime_denovo.fasta"
chimeras_denovo <- "/home/ali/Documents/simulated_data/analysis/Venn2/best/chimeras_denovo.fasta"
DADA2 <- "/home/ali/Documents/simulated_data/analysis/Venn2/best/DADA2.fasta"

rechimed_uchime_denovo <- "/home/ali/Documents/simulated_data/analysis/Venn2/best/rechimed/uchime_denovo.fasta"
rechimed_chimeras_denovo <- "/home/ali/Documents/simulated_data/analysis/Venn2/best/rechimed/chimeras_denovo.fasta"
rechimed_DADA2 <- "/home/ali/Documents/simulated_data/analysis/Venn2/best/rechimed/DADA2.fasta"

# Generate both Venn diagrams
venn1 <- generate_venn_diagram(uchime_denovo, chimeras_denovo, DADA2,
                               c("UCHIME_denovo", "Chimeras_denovo", "DADA2"),
                               c("#E69F00", "#56B4E9", "#009E73"),
                               "Original Venn Diagram")

venn2 <- generate_venn_diagram(rechimed_uchime_denovo, rechimed_chimeras_denovo, rechimed_DADA2,
                               c("ReChimed UCHIME_denovo", "ReChimed chimeras_denovo", "ReChimed DADA2"),
                               c("#FF0000", "#0000FF", "#00FF00"),
                               "ReChimed Venn Diagram")

# Combine both Venn diagrams into one SVG with increased text size
svg("/home/ali/Documents/simulated_data/analysis/Venn2/best/combined_venn_diagrams.svg", width = 20, height = 10)
grid.arrange(venn1, venn2, ncol = 2)
dev.off()
