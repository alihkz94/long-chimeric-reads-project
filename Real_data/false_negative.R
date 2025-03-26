###############################################################################
# Optimized Visualization Script for ITS Amplicon Chimera Analysis
###############################################################################

# ----------------------------- LIBRARIES --------------------------------------
library(dplyr)
library(ggplot2)
library(viridis)      # For color palettes
library(randomForest)
library(ggfortify)    # For PCA plotting
library(patchwork)    # For arranging multiple plots
library(scales)       # For better axis formatting
library(cowplot)      # For plot grid arrangement
library(tidyr)        # For data reshaping

# ----------------------------- DATA LOADING -----------------------------------
false_neg <- read.csv("false_negatives_merged.csv", stringsAsFactors = FALSE)
inferrnal <- read.csv("results_inferrnal.csv", stringsAsFactors = FALSE)

# ----------------------------- DATA PREPROCESSING -----------------------------
# Clean sequence IDs consistently
false_neg <- false_neg %>% mutate(sequence_id = sub(";size=.*", "", sequence_id))
inferrnal <- inferrnal %>% mutate(sequence_id = sub(";size=.*", "", sequence_id))

# Merge datasets and calculate key metrics
merged_data <- inner_join(false_neg, inferrnal, by = "sequence_id") %>%
  mutate(
    seq_length = nchar(sequence),
    start      = as.numeric(segment1_start),
    end        = as.numeric(segment1_end),
    rel_start  = start / seq_length
  )

# Summarize hits per sequence with improved naming
hit_summary <- merged_data %>%
  group_by(sequence_id, group) %>%
  summarise(
    n_5_8s_hits       = n(),
    avg_infernal_score = mean(score, na.rm = TRUE),
    seq_length         = first(seq_length),
    rel_start_min      = min(rel_start, na.rm = TRUE),
    rel_start_max      = max(rel_start, na.rm = TRUE),
    hit_range          = max(rel_start, na.rm = TRUE) - min(rel_start, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    hit_category = ifelse(n_5_8s_hits > 1, "Multiple", "Single")
  )

# Filter only chimeric sequences (multiple hits)
chimeric_data <- hit_summary %>% filter(hit_category == "Multiple")

# Rename algorithm groups to be more publication-friendly
chimeric_data <- chimeric_data %>%
  mutate(algorithm = case_when(
    group == "chimeras_denovo_adjusted" ~ "ChimeraX-Adjusted",
    group == "chimeras_denovo_default" ~ "ChimeraX-Default",
    group == "removeBimeraDenovo" ~ "DADA2-Chimera",
    group == "uchime_denovo_adjusted" ~ "UCHIME-Adjusted",
    group == "uchime_denovo_default" ~ "UCHIME-Default",
    TRUE ~ group
  ))

# --------------------- IMPROVED SCATTER PLOT ---------------------------------
# Add marginal histograms for better distribution visualization
library(ggExtra)

# Create the main scatter plot with enhanced aesthetics
scatter_plot <- ggplot(chimeric_data, aes(x = seq_length, y = hit_range)) +
  geom_point(aes(color = algorithm), alpha = 0.8, size = 2) +
  scale_y_log10(
    labels = scales::label_scientific(),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    oob = scales::squish, # Handle out-of-bounds values
    limits = c(NA, NA)
  ) +
  scale_color_viridis_d(option = "D", begin = 0.1, end = 0.9, name = "Algorithm") +
  scale_x_continuous(labels = scales::comma, limits = c(0, max(chimeric_data$seq_length) * 1.05)) +
  labs(
    x = "Sequence Length (bp)",
    y = "Range of Relative 5.8S Start Positions (log scale)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# Add marginal histograms
scatter_plot_with_marginals <- ggMarginal(
  scatter_plot, 
  type = "histogram",
  fill = "grey70",
  bins = 30,
  size = 10
)

# --------------------- IMPROVED INFERNAL SCORES PLOT -------------------------
# Create violin plots with overlaid boxplots for better distribution visualization
infernal_plot <- ggplot(chimeric_data, 
                        aes(x = reorder(algorithm, avg_infernal_score, FUN = median), 
                            y = avg_infernal_score)) +
  geom_violin(aes(fill = algorithm), alpha = 0.6, trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7, fill = "white") +
  geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.5) +
  scale_fill_viridis_d(option = "A", begin = 0.1, end = 0.9, guide = "none") +
  labs(
    x = NULL,
    y = "Average Infernal Score"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")
  )

# --------------------- IMPROVED DENSITY PLOT ---------------------------------
# Use a more effective comparison layout with filled density plots
# Convert to long format for better plotting
density_data <- chimeric_data %>%
  dplyr::select(sequence_id, algorithm, rel_start_min) %>%
  filter(!is.na(rel_start_min))

# Calculate kernel density estimates for smoother visualization
density_estimates <- density_data %>%
  group_by(algorithm) %>%
  summarize(
    count = n(),
    .groups = "drop"
  )

# Create improved density plot with overlaid quantile markers
density_plot <- ggplot(density_data, aes(x = rel_start_min, fill = algorithm)) +
  geom_density(alpha = 0.6, adjust = 1.5) +
  geom_rug(aes(color = algorithm), alpha = 0.6, size = 0.5) +
  scale_fill_viridis_d(option = "D", begin = 0.1, end = 0.9) +
  scale_color_viridis_d(option = "D", begin = 0.1, end = 0.9) +
  facet_wrap(~ algorithm, ncol = 1, scales = "free_y") +
  labs(
    x = "Minimum Relative 5.8S Start Position",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )

# ---------------------- IMPROVED PCA VISUALIZATION ---------------------------
# Calculate more interpretable metrics for chimeric severity
z_hit_range <- if(sd(chimeric_data$hit_range, na.rm = TRUE) == 0) {
  0
} else {
  (chimeric_data$hit_range - mean(chimeric_data$hit_range, na.rm = TRUE)) /
    sd(chimeric_data$hit_range, na.rm = TRUE)
}

z_infernal <- if(sd(chimeric_data$avg_infernal_score, na.rm = TRUE) == 0) {
  0
} else {
  (chimeric_data$avg_infernal_score - mean(chimeric_data$avg_infernal_score, na.rm = TRUE)) /
    sd(chimeric_data$avg_infernal_score, na.rm = TRUE)
}

# Calculate chimeric severity metric combining multiple features
chimeric_data <- chimeric_data %>%
  mutate(
    severity_score = (z_hit_range + z_infernal) / 2,
    severity = ifelse(severity_score > median(severity_score, na.rm = TRUE), "High", "Low")
  )

# PCA analysis with selected features
pca_features <- chimeric_data %>%
  select(avg_infernal_score, hit_range, seq_length) %>%
  # Handle NA values properly
  na.omit()

# Remove constant columns
pca_features <- pca_features[, sapply(pca_features, function(x) sd(x, na.rm = TRUE)) > 0]

# Run PCA with proper scaling
pca_result <- prcomp(pca_features, scale. = TRUE, center = TRUE)
pca_data <- as.data.frame(pca_result$x) %>%
  # Add back metadata
  mutate(
    severity = chimeric_data$severity[row.names(.) %in% row.names(chimeric_data)],
    algorithm = chimeric_data$algorithm[row.names(.) %in% row.names(chimeric_data)]
  )

# Calculate variance explained
var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)

# Create improved PCA plot with ellipses
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  # Add ellipses to show group boundaries
  stat_ellipse(aes(fill = severity), geom = "polygon", alpha = 0.2) +
  geom_point(aes(color = severity), size = 2.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", size = 0.5) +
  scale_color_manual(values = c("High" = "#440154", "Low" = "#21908C"), 
                     name = "Chimeric\nSeverity") +
  scale_fill_manual(values = c("High" = "#440154", "Low" = "#21908C"), 
                    name = "Chimeric\nSeverity") +
  labs(
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Add the feature loadings (biplot arrows) with more subtle appearance
loadings <- as.data.frame(pca_result$rotation)
loadings$feature <- rownames(loadings)

# Scale loadings for better visibility
loading_scale <- 3.5
loadings$PC1 <- loadings$PC1 * loading_scale
loadings$PC2 <- loadings$PC2 * loading_scale

pca_plot <- pca_plot +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "darkred", size = 0.6) +
  geom_text(data = loadings, 
            aes(x = PC1 * 1.1, y = PC2 * 1.1, label = feature),
            color = "darkred", size = 3.5, fontface = "bold")

# --------------------- NEW SWARM PLOT FOR ALGORITHM COMPARISON --------------
# Create a swarm plot to better visualize algorithm performance
library(ggbeeswarm)

# Calculate algorithm performance metrics
algorithm_metrics <- chimeric_data %>%
  group_by(algorithm) %>%
  summarize(
    mean_hit_range = mean(hit_range, na.rm = TRUE),
    mean_score = mean(avg_infernal_score, na.rm = TRUE),
    high_severity_pct = 100 * sum(severity == "High", na.rm = TRUE) / n(),
    n_sequences = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_score))

# Before creating the plot, check for groups with insufficient data
algorithm_counts <- table(chimeric_data$algorithm)
print(algorithm_counts)

# Filter out algorithms with insufficient data (if needed)
valid_algorithms <- names(algorithm_counts[algorithm_counts >= 2])
chimeric_data_filtered <- chimeric_data %>% 
  filter(algorithm %in% valid_algorithms)

# Then use chimeric_data_filtered instead of chimeric_data for the algorithm_comparison plot

# Create performance comparison plot
algorithm_comparison <- ggplot(chimeric_data, 
                               aes(x = reorder(algorithm, avg_infernal_score, FUN = median), 
                                   y = hit_range)) +
  # Check for sufficient data points per group before using quasirandom
  geom_point(aes(color = severity), alpha = 0.7, size = 2, position = position_jitter(width = 0.2)) +
  scale_y_log10(labels = scales::label_scientific(), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                oob = scales::squish) + # Handle out-of-bounds values
  scale_color_manual(values = c("High" = "#440154", "Low" = "#21908C")) +
  labs(
    x = NULL,
    y = "Hit Range (log scale)",
    color = "Chimeric Severity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# --------------------- SUMMARY STATISTICS TABLE -----------------------------
# Create a summary statistics table
summary_stats <- chimeric_data %>%
  group_by(algorithm) %>%
  summarize(
    n_sequences = n(),
    median_hit_range = median(hit_range, na.rm = TRUE),
    median_score = median(avg_infernal_score, na.rm = TRUE),
    avg_seq_length = mean(seq_length, na.rm = TRUE),
    high_severity_pct = 100 * sum(severity == "High", na.rm = TRUE) / n(),
    .groups = "drop"
  ) %>%
  arrange(desc(median_score))

# --------------------- COMBINED FIGURE ARRANGEMENT -------------------------
# Create figure labels for publication
fig_a_label <- ggdraw() + draw_label("A", fontface = "bold", size = 16, x = 0, hjust = 0)
fig_b_label <- ggdraw() + draw_label("B", fontface = "bold", size = 16, x = 0, hjust = 0)
fig_c_label <- ggdraw() + draw_label("C", fontface = "bold", size = 16, x = 0, hjust = 0)
fig_d_label <- ggdraw() + draw_label("D", fontface = "bold", size = 16, x = 0, hjust = 0)
fig_e_label <- ggdraw() + draw_label("E", fontface = "bold", size = 16, x = 0, hjust = 0)

# Arrange plots with improved layout - create a 2x2 grid with the new comparison plot
panel_a <- plot_grid(fig_a_label, scatter_plot_with_marginals, ncol = 1, rel_heights = c(0.1, 1))
panel_b <- plot_grid(fig_b_label, infernal_plot, ncol = 1, rel_heights = c(0.1, 1))
panel_c <- plot_grid(fig_c_label, pca_plot, ncol = 1, rel_heights = c(0.1, 1))
panel_d <- plot_grid(fig_d_label, density_plot, ncol = 1, rel_heights = c(0.1, 1))
panel_e <- plot_grid(fig_e_label, algorithm_comparison, ncol = 1, rel_heights = c(0.1, 1))

# Top row contains scatter and infernal plots
top_row <- plot_grid(panel_a, panel_b, ncol = 2, rel_widths = c(1.2, 0.8))

# Middle row contains PCA plot and new algorithm comparison
middle_row <- plot_grid(panel_c, panel_e, ncol = 2)

# Create a more balanced layout
combined_plots <- plot_grid(
  top_row,
  middle_row,
  panel_d,
  ncol = 1,
  rel_heights = c(1, 1, 1.2)
)

# Add overall title
title <- ggdraw() + 
  draw_label(
    "Comparative Analysis of ITS Chimera Detection Algorithms",
    fontface = "bold",
    size = 16,
    x = 0.5,
    hjust = 0.5
  )

final_figure <- plot_grid(title, combined_plots, ncol = 1, rel_heights = c(0.1, 1))

# Save the optimized plots
ggsave("figure1_optimized.pdf", final_figure, width = 12, height = 16, dpi = 300)
ggsave("figure1_optimized.png", final_figure, width = 12, height = 16, dpi = 300)

# Save individual panels for supplementary or alternative arrangements
ggsave("figure1A_sequence_length_vs_hit_range.pdf", panel_a, width = 7, height = 5, dpi = 300)
ggsave("figure1B_infernal_scores.pdf", panel_b, width = 5, height = 5, dpi = 300)
ggsave("figure1C_pca_features.pdf", panel_c, width = 6, height = 5, dpi = 300)
ggsave("figure1D_start_position_density.pdf", panel_d, width = 7, height = 8, dpi = 300)
ggsave("figure1E_algorithm_comparison.pdf", panel_e, width = 6, height = 5, dpi = 300)

# Export summary statistics as CSV
write.csv(summary_stats, "algorithm_summary_statistics.csv", row.names = FALSE)

# Function to print a nicely formatted summary of findings
print_key_findings <- function() {
  cat("Key Findings from ITS Chimera Detection Analysis:\n")
  cat("------------------------------------------------\n")
  
  best_algo <- summary_stats$algorithm[which.max(summary_stats$median_score)]
  worst_algo <- summary_stats$algorithm[which.min(summary_stats$median_score)]
  
  cat("1. ", best_algo, " produced the highest median Infernal scores (", 
      round(max(summary_stats$median_score), 2), ")\n", sep="")
  
  cat("2. ", worst_algo, " produced the lowest median Infernal scores (", 
      round(min(summary_stats$median_score), 2), ")\n", sep="")
  
  cat("3. The PCA analysis indicates that sequence length, hit range, and Infernal scores\n",
      "   collectively explain ", sum(var_explained), "% of the variation in the data\n", sep="")
  
  cat("4. The algorithms identified between ", min(summary_stats$n_sequences), 
      " and ", max(summary_stats$n_sequences), " potential chimeric sequences\n", sep="")
}

# Print summary of key findings
print_key_findings()
