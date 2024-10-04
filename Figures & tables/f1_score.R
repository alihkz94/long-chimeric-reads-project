# Load necessary packages
library(ggplot2)

# Define the data without the UCHIME_denovo mediocre module
data <- data.frame(
  module = factor(c("UCHIME_denovo default", "DADA2 default", "Chimeras_denovo default",
                   "UCHIME_denovo customized", "Chimeras_denovo customized"),
                 levels = c("UCHIME_denovo default", "DADA2 default", "Chimeras_denovo default",
                           "UCHIME_denovo customized", "Chimeras_denovo customized")),
  F1_score = c(0.8009, 0.14, 0.088,
               0.92, 0.48)
)

# Create the plot
p <- ggplot(data, aes(x = module, y = F1_score, fill = module)) +
  geom_point(shape = 21, size = 6) +
  labs(title = "F1 Scores of Modules",
       x = "Module",
       y = "F1 Score") +
  scale_color_manual(values = c("brown", "green", "black",
                               "blue", "purple")) +
  scale_fill_manual(values = c("brown", "green", "black",
                               "blue", "purple")) +
  geom_text(aes(label = round(F1_score, 2)), vjust = -1.2, size = 3.5) + # Adjusted text size for clarity
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"), # Rotate x-axis labels for better fitting
        plot.margin = margin(12, 12, 12, 12, "mm")) + # Adjust plot margins
  ylim(0, 1.2) # Extend the y-axis

# Optionally, save as PNG for digital uses
ggsave("F1_Scores_modules_no_uchime_mediocre.png", plot = p, width = 8, height = 6, dpi = 1200)
