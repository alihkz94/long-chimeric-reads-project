library(dplyr)
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
