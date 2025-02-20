# Load necessary packages
library(ggplot2)
library(dplyr)

# Define the data
data <- data.frame(
  module = factor(c("uchime_denovo default", "removeBimeraDenovo default", "chimeras_denovo default",
                    "uchime_denovo adjusted", "chimeras_denovo adjusted")),
  algorithm = c("uchime_denovo", "removeBimeraDenovo", "chimeras_denovo",
                "uchime_denovo", "chimeras_denovo"),
  F1_score = c(0.87, 0.14, 0.15,
               0.89, 0.71)
)

# Sort data by F1_score in descending order
data <- data %>% arrange(desc(F1_score))

# ***REORDER THE FACTOR LEVELS BASED ON F1_SCORE***
data$module <- factor(data$module, levels = data$module[order(desc(data$F1_score))])

# Define a fixed color palette for each algorithm
algorithm_palette <- c(
  "uchime_denovo" = "red",
  "removeBimeraDenovo" = "green",
  "chimeras_denovo" = "purple"
)

# Create the plot
p <- ggplot(data, aes(x = module, y = F1_score, fill = algorithm)) +
  geom_point(shape = 21, size = 6) +
  labs(title = "F1 Scores of Modules",
       x = "Module",
       y = "F1 Score") +
  scale_fill_manual(values = algorithm_palette) +
  geom_text(aes(label = round(F1_score, 2)), vjust = -1.2, size = 3.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
        plot.margin = margin(12, 12, 12, 12, "mm"),
        plot.title = element_text(hjust = 0.5, vjust = 1)) +
  ylim(0, 1.2)

# Optionally, save as PNG for digital uses
ggsave("F1_Scores_modules.png", plot = p, width = 8, height = 6, dpi = 1200)
