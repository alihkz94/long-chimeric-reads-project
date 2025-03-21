# Load necessary libraries
library(dplyr)      # Data manipulation
library(ggplot2)    # Visualization
library(car)        # For residual diagnostics
library(ggpubr)     # For enhanced plotting (optional)

# Read data
data <- read.csv("Jamy__Procrustes_residuals_VsRaw.csv")

# Convert relevant variables to factors
data$table2 <- factor(data$table2)           # Chimera removal method
data$isolation_source <- factor(data$isolation_source)

### Community Richness Analysis ###

# Visualize the distribution of richness
hist(data$Richness, main = "Histogram of Community Richness", xlab = "Richness")

# Fit a linear model for richness
richness_model <- lm(Richness ~ table2 * isolation_source, data = data)
summary(richness_model)

# Check model residuals for normality (Shapiro-Wilk test)
richness_res <- residuals(richness_model)
shapiro_result <- shapiro.test(richness_res)
print(shapiro_result)

# If residuals are non-normal (p-value < 0.05), consider non-parametric tests:
kruskal_richness_table2 <- kruskal.test(Richness ~ table2, data = data)
kruskal_richness_isolation <- kruskal.test(Richness ~ isolation_source, data = data)
print(kruskal_richness_table2)
print(kruskal_richness_isolation)

# Boxplot for community richness
ggplot(data, aes(x = table2, y = Richness, fill = isolation_source)) +
  geom_boxplot() +
  labs(x = "Chimera Removal Method",
       y = "Community Richness",
       title = "Community Richness by Chimera Removal Method and Isolation Source") +
  theme_minimal()

### Shannon Diversity Analysis ###

# Visualize the distribution of Shannon diversity
hist(data$Shannon, main = "Histogram of Shannon Diversity", xlab = "Shannon Index")

# Fit a linear model for Shannon diversity
shannon_model <- lm(Shannon ~ table2 * isolation_source, data = data)
summary(shannon_model)

# Check model residuals for normality
shannon_res <- residuals(shannon_model)
shapiro_result_shannon <- shapiro.test(shannon_res)
print(shapiro_result_shannon)

# If residuals are non-normal, consider non-parametric tests:
kruskal_shannon_table2 <- kruskal.test(Shannon ~ table2, data = data)
kruskal_shannon_isolation <- kruskal.test(Shannon ~ isolation_source, data = data)
print(kruskal_shannon_table2)
print(kruskal_shannon_isolation)

# Boxplot for Shannon diversity
ggplot(data, aes(x = table2, y = Shannon, fill = isolation_source)) +
  geom_boxplot() +
  labs(x = "Chimera Removal Method",
       y = "Shannon Diversity Index",
       title = "Shannon Diversity by Chimera Removal Method and Isolation Source") +
  theme_minimal()
