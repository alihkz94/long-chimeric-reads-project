# Set working directory (adjust to your path)
setwd("~/Documents/simulated_data/analysis/figures/upset_plot")

# Load required packages
library(ComplexUpset)
library(tidyverse)
library(Biostrings)
library(scales)  # for comma-formatting large numbers

##############################################################################
# 1) Read metadata: sample name + substrate, no header
##############################################################################
metadata <- read.csv("metadata.csv", header = FALSE, stringsAsFactors = FALSE)
colnames(metadata) <- c("SampleName", "Substrate")

##############################################################################
# 2) Define your method folders
##############################################################################
method_folders <- c(
  "uchime_default",
  "uchime_adjusted",
  "chimeras_denovo_default",
  "chimeras_denovo_adjusted",
  "removeBimeraDenovo"  # or "DADA2" if that's your correct folder
)

##############################################################################
# 3) Sequentially read FASTA files with error handling (to avoid crashes)
##############################################################################
all_data_list <- list()

for (method_folder in method_folders) {
  for (i in seq_len(nrow(metadata))) {
    sample_name <- metadata$SampleName[i]
    substrate   <- metadata$Substrate[i]
    fasta_file  <- file.path(method_folder, paste0(sample_name, ".fasta"))
    
    if (file.exists(fasta_file)) {
      tryCatch({
        seqs <- readDNAStringSet(fasta_file)
        seqs_char <- unique(as.character(seqs))
        
        df <- data.frame(
          SampleName = sample_name,
          Substrate  = substrate,
          Sequence   = seqs_char,
          Method     = method_folder,
          stringsAsFactors = FALSE
        )
        all_data_list[[length(all_data_list) + 1]] <- df
      }, error = function(e) {
        message(sprintf("Error reading file %s: %s", fasta_file, e$message))
      })
    }
  }
}

all_data <- bind_rows(all_data_list)

##############################################################################
# 4) Pivot to wide format (presence/absence)
##############################################################################
all_data <- all_data %>% mutate(Presence = TRUE)
df_wide <- all_data %>%
  pivot_wider(
    id_cols     = c("SampleName", "Substrate", "Sequence"),
    names_from  = "Method",
    values_from = "Presence",
    values_fill = FALSE
  )

##############################################################################
# 5) Produce the UpSet plot with improved aesthetics and SVG output
##############################################################################
methods_for_upset <- c(
  "chimeras_denovo_default",
  "removeBimeraDenovo",
  "chimeras_denovo_adjusted",
  "uchime_adjusted",
  "uchime_default"
)

# Build the UpSet plot with corrected aesthetics
upset_plot <- upset(
  data = df_wide,
  intersect = methods_for_upset,
  
  # Intersection size bar (top) with corrected number labels
  base_annotations = list(
    'intersection_size' = intersection_size(
      counts = TRUE,
      bar_number_threshold = 0,
      text = list(
        angle = 45,    # Rotate text
        vjust = -0.5,  # Position above bars
        hjust = -0.1,  # Adjust horizontal position
        size = 3.5     # Adjust text size
      )
    ) +
      scale_fill_manual(values = "black", guide = "none") +
      scale_y_continuous(
        "Number of (sample–sequence) entries shared",
        labels = comma,
        # Extend y-axis limits to accommodate labels
        expand = expansion(mult = c(0, 0.15))
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, hjust = 0.5)
      ) +
      ggtitle("Chimera Filtering Methods")
  ),
  
  # Set size bar (left) with improved number labels
  set_sizes = (
    upset_set_size(
      geom = geom_bar(fill = "blue")
    ) +
      # Add white number labels centered on bars
      geom_text(
        aes(y = after_stat(count), label = after_stat(count)),
        stat = "count",
        color = "white",
        position = position_stack(vjust = 0.5),
        size = 3.5
      ) +
      scale_y_continuous(
        labels = comma,
        name = "Set size"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  ),
  
  sort_intersections = "descending",
  
  # Matrix theme
  themes = upset_modify_themes(
    list(
      'intersections_matrix' = theme(
        text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    )
  ),
  
  # Simplified matrix parameters
  matrix = intersection_matrix(
    geom = geom_point(size = 3),
    segment = geom_segment(size = 0.8)
  )
)

# Save with adjusted dimensions
ggsave(
  filename = "chimera_upset.svg",
  plot = upset_plot,
  device = "svg",
  width = 15,
  height = 8,
  dpi = 300
)
