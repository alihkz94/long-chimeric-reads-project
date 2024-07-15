######################################################################################################################################################################
#The scripts are for the box plot generations of chimeric reads and also statical test analysis to see how different methods vary in the detection of chimeric reads## 
######################################################################################################################################################################


#####chimeric reads boxplot for custom settings####
# Load the necessary libraries
library(ggplot2)
library(reshape2)
library(scales)

# Define the data frame
data <- data.frame(
  FILE = c("ERR6454461.fasta", "ERR6454462.fasta", "ERR6454463.fasta", "ERR6454464.fasta", 
           "ERR6454465.fasta", "ERR6454466.fasta", "ERR6454467.fasta", "ERR6454468.fasta", 
           "ERR6454469.fasta", "ERR6454470.fasta", "ERR6454471.fasta", "ERR6454472.fasta", 
           "ERR6454473.fasta", "ERR6454474.fasta", "ERR6454475.fasta", "ERR6454476.fasta", 
           "ERR6454477.fasta", "ERR6454478.fasta"),
  Uchime_denovo = c(1285, 316, 367, 3732, 1099, 884, 1454, 3122, 3004, 10218, 
                    3077, 1482, 25588, 1483, 3012, 4807, 4768, 619),
  Chimeras_denovo = c(6601, 3081, 4979, 7908, 14540, 3350, 6956, 16298, 18552, 110918, 
                      4136, 5715, 24044, 9558, 18119, 10859, 14591, 6818),
  DADA2 = c(20593, 10245, 13502, 26106, 45094, 9758, 19784, 45443, 51898, 303374, 
            15053, 18272, 63023, 29093, 57004, 33134, 43393, 18486)
)

# Melt the data frame to long format for ggplot2
data_long <- melt(data, id.vars = "FILE", variable.name = "Method", value.name = "Count")

# Create the boxplot with log scale and enhanced styling
boxplot <- ggplot(data_long, aes(x = Method, y = Count, fill = Method)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10(labels = comma) +  # Use logarithmic scale and format labels
  labs(title = "Chimeric reads abudances for custom settings",
       x = "Method",
       y = "Counts") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.1, color = "gray"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  scale_fill_brewer(palette = "Set2") +  # Use a color palette from RColorBrewer
  geom_jitter(width = 0.2, size = 1, alpha = 0.6)  # Add jitter for data points

# Display the plot
print(boxplot)


##### chimeric reads box plot for deafault settings####
# Load the necessary libraries
library(ggplot2)
library(reshape2)
library(scales)

# Define the new data frame
new_data <- data.frame(
  FILE = c("ERR6454461.fasta", "ERR6454462.fasta", "ERR6454463.fasta", "ERR6454464.fasta", 
           "ERR6454465.fasta", "ERR6454466.fasta", "ERR6454467.fasta", "ERR6454468.fasta", 
           "ERR6454469.fasta", "ERR6454470.fasta", "ERR6454471.fasta", "ERR6454472.fasta", 
           "ERR6454473.fasta", "ERR6454474.fasta", "ERR6454475.fasta", "ERR6454476.fasta", 
           "ERR6454477.fasta", "ERR6454478.fasta"),
  Uchime_denovo = c(307, 50, 76, 345, 37, 100, 753, 253, 697, 2414, 
                    381, 502, 9213, 136, 560, 1220, 616, 121),
  Chimeras_denovo = c(21856, 7429, 12978, 33667, 50519, 13357, 19077, 45361, 44504, 185475, 
                      13411, 11075, 53684, 19796, 39870, 34560, 37460, 13800),
  DADA2 = c(20593, 10245, 13502, 26106, 45094, 9758, 19784, 45443, 51898, 303374, 
            15053, 18272, 63023, 29093, 57004, 33134, 43393, 18486)
)

# Melt the new data frame to long format for ggplot2
new_data_long <- melt(new_data, id.vars = "FILE", variable.name = "Method", value.name = "Count")

# Create the boxplot with log scale and enhanced styling for the new data
default_boxplot <- ggplot(new_data_long, aes(x = Method, y = Count, fill = Method)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10(labels = comma) +  # Use logarithmic scale and format labels
  labs(title = "Chimeric reads abudances for default settings",
       x = "Method",
       y = "Count") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.1, color = "gray"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  scale_fill_brewer(palette = "Set2") +  # Use a color palette from RColorBrewer
  geom_jitter(width = 0.2, size = 1, alpha = 0.6)  # Add jitter for data points

# Display the plot
print(default_boxplot)


##### compare chimeric reads abudances for default settings ####
# Define the new data frame
data <- data.frame(
  FILE = c("ERR6454461.fasta", "ERR6454462.fasta", "ERR6454463.fasta", "ERR6454464.fasta", 
           "ERR6454465.fasta", "ERR6454466.fasta", "ERR6454467.fasta", "ERR6454468.fasta", 
           "ERR6454469.fasta", "ERR6454470.fasta", "ERR6454471.fasta", "ERR6454472.fasta", 
           "ERR6454473.fasta", "ERR6454474.fasta", "ERR6454475.fasta", "ERR6454476.fasta", 
           "ERR6454477.fasta", "ERR6454478.fasta"),
  Uchime_denovo = c(307, 50, 76, 345, 37, 100, 753, 253, 697, 2414, 
                    381, 502, 9213, 136, 560, 1220, 616, 121),
  Chimeras_denovo = c(21856, 7429, 12978, 33667, 50519, 13357, 19077, 45361, 44504, 185475, 
                      13411, 11075, 53684, 19796, 39870, 34560, 37460, 13800),
  DADA2 = c(20593, 10245, 13502, 26106, 45094, 9758, 19784, 45443, 51898, 303374, 
            15053, 18272, 63023, 29093, 57004, 33134, 43393, 18486)
)

# Reshape data for Kruskal-Wallis test
data_long <- data %>%
  pivot_longer(cols = -FILE, names_to = "Method", values_to = "Values")

# Perform Kruskal-Wallis test for all categories
kruskal_test_result <- kruskal.test(Values ~ Method, data = data_long)

# Print the result of the Kruskal-Wallis test
print(kruskal_test_result)

# Function to perform pairwise Wilcoxon test
pairwise_wilcox_test <- function(data, group_col, value_col) {
  combinations <- combn(unique(data[[group_col]]), 2)
  results <- apply(combinations, 2, function(pair) {
    group1 <- data %>% filter(!!sym(group_col) == pair[1])
    group2 <- data %>% filter(!!sym(group_col) == pair[2])
    test_result <- wilcox.test(group1[[value_col]], group2[[value_col]])
    p_value <- test_result$p.value
    list(pair = paste(pair[1], "vs", pair[2]), p_value = p_value)
  })
  do.call(rbind, lapply(results, as.data.frame))
}

# Perform pairwise Wilcoxon tests
pairwise_results <- pairwise_wilcox_test(data_long, "Method", "Values")
pairwise_results <- as.data.frame(pairwise_results)
pairwise_results$p_value <- as.numeric(as.character(pairwise_results$p_value))
pairwise_results$Significance <- ifelse(pairwise_results$p_value < 0.05, "Significant", "Not Significant")

# Print the pairwise results
print(pairwise_results)

# Interpretation based on the results
if (kruskal_test_result$p.value < 0.05) {
  cat("The differences between the three categories are statistically significant.\n")
} else {
  cat("The differences between the three categories are not statistically significant.\n")
}

# Interpretation for pairwise comparisons
apply(pairwise_results, 1, function(row) {
  cat(paste("Comparison between", row["pair"], "is", row["Significance"], "with p-value:", row["p_value"], "\n"))
})


#####compare chimeric reads abundances for custom settings####
library(tidyverse)

# Create the data frame
data <- data.frame(
  FILE = c("ERR6454461.fasta", "ERR6454462.fasta", "ERR6454463.fasta", "ERR6454464.fasta", 
           "ERR6454465.fasta", "ERR6454466.fasta", "ERR6454467.fasta", "ERR6454468.fasta", 
           "ERR6454469.fasta", "ERR6454470.fasta", "ERR6454471.fasta", "ERR6454472.fasta", 
           "ERR6454473.fasta", "ERR6454474.fasta", "ERR6454475.fasta", "ERR6454476.fasta", 
           "ERR6454477.fasta", "ERR6454478.fasta"),
  Uchime_denovo = c(1285, 316, 367, 3732, 1099, 884, 1454, 3122, 3004, 10218, 
                    3077, 1482, 25588, 1483, 3012, 4807, 4768, 619),
  Chimeras_denovo = c(6601, 3081, 4979, 7908, 14540, 3350, 6956, 16298, 18552, 110918, 
                      4136, 5715, 24044, 9558, 18119, 10859, 14591, 6818),
  DADA2 = c(20593, 10245, 13502, 26106, 45094, 9758, 19784, 45443, 51898, 303374, 
            15053, 18272, 63023, 29093, 57004, 33134, 43393, 18486)
)

# Reshape data for Kruskal-Wallis test
data_long <- data %>%
  pivot_longer(cols = -FILE, names_to = "Method", values_to = "Values")

# Perform Kruskal-Wallis test for all categories
kruskal_test_result <- kruskal.test(Values ~ Method, data = data_long)

# Print the result of the Kruskal-Wallis test
print(kruskal_test_result)

# Function to perform pairwise Wilcoxon test
pairwise_wilcox_test <- function(data, group_col, value_col) {
  combinations <- combn(unique(data[[group_col]]), 2)
  results <- apply(combinations, 2, function(pair) {
    group1 <- data %>% filter(!!sym(group_col) == pair[1])
    group2 <- data %>% filter(!!sym(group_col) == pair[2])
    test_result <- wilcox.test(group1[[value_col]], group2[[value_col]])
    p_value <- test_result$p.value
    list(pair = paste(pair[1], "vs", pair[2]), p_value = p_value)
  })
  do.call(rbind, lapply(results, as.data.frame))
}

# Perform pairwise Wilcoxon tests
pairwise_results <- pairwise_wilcox_test(data_long, "Method", "Values")
pairwise_results <- as.data.frame(pairwise_results)
pairwise_results$p_value <- as.numeric(as.character(pairwise_results$p_value))
pairwise_results$Significance <- ifelse(pairwise_results$p_value < 0.05, "Significant", "Not Significant")

# Print the pairwise results
print(pairwise_results)

# Interpretation based on the results
if (kruskal_test_result$p.value < 0.05) {
  cat("The differences between the three categories are statistically significant.\n")
} else {
  cat("The differences between the three categories are not statistically significant.\n")
}

# Interpretation for pairwise comparisons
apply(pairwise_results, 1, function(row) {
  cat(paste("Comparison between", row["pair"], "is", row["Significance"], "with p-value:", row["p_value"], "\n"))
})


#####compare number of nonchimeric+blastn recovered reads together####

library(tidyverse)

# Updated data frame with new numbers
data <- data.frame(
  FILE = c("ERR6454461.fasta", "ERR6454462.fasta", "ERR6454463.fasta", "ERR6454464.fasta", 
           "ERR6454465.fasta", "ERR6454466.fasta", "ERR6454467.fasta", "ERR6454468.fasta", 
           "ERR6454469.fasta", "ERR6454470.fasta", "ERR6454471.fasta", "ERR6454472.fasta", 
           "ERR6454473.fasta", "ERR6454474.fasta", "ERR6454475.fasta", "ERR6454476.fasta", 
           "ERR6454477.fasta", "ERR6454478.fasta"),
  Uchime_denovo = c(176676, 63706, 65934, 237814, 290092, 104990, 111013, 327232, 255230, 1003766, 
                    153725, 68738, 212383, 211269, 192938, 193873, 216363, 54320),
  Chimeras_denovo = c(174795, 63606, 64185, 237089, 285813, 104162, 109186, 325685, 255447, 939056, 
                      152345, 64341, 218874, 204087, 177562, 194670, 206370, 54016),
  DADA2 = c(263934, 100016, 102723, 384127, 471039, 171681, 179385, 547796, 426004, 1103263, 
            232323, 88413, 320800, 264681, 248409, 318179, 312409, 86740)
)

# Reshape data for Kruskal-Wallis test
data_long <- data %>%
  pivot_longer(cols = -FILE, names_to = "Method", values_to = "Values")

# Perform Kruskal-Wallis test for all categories
kruskal_test_result <- kruskal.test(Values ~ Method, data = data_long)

# Print the result of the Kruskal-Wallis test
print(kruskal_test_result)

# Function to perform pairwise Wilcoxon test
pairwise_wilcox_test <- function(data, group_col, value_col) {
  combinations <- combn(unique(data[[group_col]]), 2)
  results <- apply(combinations, 2, function(pair) {
    group1 <- data %>% filter(!!sym(group_col) == pair[1])
    group2 <- data %>% filter(!!sym(group_col) == pair[2])
    test_result <- wilcox.test(group1[[value_col]], group2[[value_col]])
    p_value <- test_result$p.value
    list(pair = paste(pair[1], "vs", pair[2]), p_value = p_value)
  })
  do.call(rbind, lapply(results, as.data.frame))
}

# Perform pairwise Wilcoxon tests
pairwise_results <- pairwise_wilcox_test(data_long, "Method", "Values")
pairwise_results <- as.data.frame(pairwise_results)
pairwise_results$p_value <- as.numeric(as.character(pairwise_results$p_value))
pairwise_results$Significance <- ifelse(pairwise_results$p_value < 0.05, "Significant", "Not Significant")

# Print the pairwise results
print(pairwise_results)

# Interpretation based on the results
if (kruskal_test_result$p.value < 0.05) {
  cat("The differences between the three categories are statistically significant.\n")
} else {
  cat("The differences between the three categories are not statistically significant.\n")
}

# Interpretation for pairwise comparisons
apply(pairwise_results, 1, function(row) {
  cat(paste("Comparison between", row["pair"], "is", row["Significance"], "with p-value:", row["p_value"], "\n"))
})

