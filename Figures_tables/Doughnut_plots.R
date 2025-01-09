######### chimeric reads #########
#Libraries
library(ggplot2)
library(svglite)
library(scales)  # For percentage formatting
library(patchwork)  # For combining plots

# Data
chimera_data <- data.frame(
  Module = c("UCHIME_denovo", "Chimeras_denovo", "DADA2"),
  False_positive_chimeras = c(4521, 121963, 176001),
  Absolute_chimeras = c(34469, 89136, 143066),
  Borderline = c(2097, 576, 965)
)

# Function to create a doughnut chart with the module name in the center
create_doughnut_chart <- function(module_data, module_name) {
  # Reshape the data for the doughnut chart
  pie_data <- data.frame(
    Category = c("False Positive Chimeras", "Absolute Chimeras", "Borderline"),
    Value = c(module_data$False_positive_chimeras, module_data$Absolute_chimeras, module_data$Borderline)
  )
  
  # Calculate percentage for labels
  pie_data$Percentage <- pie_data$Value / sum(pie_data$Value) * 100
  pie_data$Label <- paste0(round(pie_data$Percentage, 1), "%")  # Format percentage labels
  
  # Create the doughnut chart
  doughnut_chart <- ggplot(pie_data, aes(x = 2, y = Value, fill = Category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +  # Create space in the center for the doughnut effect
    geom_text(aes(label = Label), position = position_stack(vjust = 0.5), color = "white", size = 4) +  # Add percentage labels
    theme_void() +  # Remove axes and gridlines
    theme(legend.position = "none") +  # Hide the legend since categories will be labeled outside the chart
    scale_fill_manual(values = c("#F8766D", "#00BFC4", "#7CAE00")) +  # Custom colors for the slices
    annotate("text", x = 0, y = 0, label = module_name, size = 6, fontface = "bold")  # Add module name in center
  
  return(doughnut_chart)
}

# Create a list to hold the individual plots
doughnut_plots <- list()

# Loop over each module and create corresponding doughnut chart
for (i in 1:nrow(chimera_data)) {
  doughnut_plots[[i]] <- create_doughnut_chart(chimera_data[i, ], chimera_data$Module[i])
}

# Combine all three doughnut charts horizontally using patchwork
combined_chart <- doughnut_plots[[1]] + doughnut_plots[[2]] + doughnut_plots[[3]] + 
  plot_layout(ncol = 3) + 
  plot_annotation(title = "Chimera Distribution Across Modules") & 
  theme(plot.title = element_text(hjust = 0.5))

# Add the category names at the bottom
combined_chart <- combined_chart + 
  plot_spacer() + 
  ggplot() + 
  annotate("text", x = 1:3, y = 0, label = c("False Positive Chimeras", "Absolute Chimeras", "Borderline"), size = 5) + 
  theme_void()

# Save the combined chart as an SVG
svg_filename <- "combined_doughnut_chart.svg"
svglite(svg_filename, width = 13.33, height = 7.5)  # PowerPoint slide dimensions
print(combined_chart)
dev.off()



####### non-chimeric ########
# Required Libraries
library(ggplot2)
library(svglite)
library(scales)  # For percentage formatting
library(patchwork)  # For combining plots

# Data
chimera_data <- data.frame(
  Module = c("UCHIME_denovo", "Chimeras_denovo", "DADA2"),
  False_positive_chimeras = c(982444, 889302, 866957),
  Absolute_chimeras = c(980821, 915282, 834527),
  Borderline = c(47088, 48355, 49160)
)

# Function to create a doughnut chart with the module name in the center
create_doughnut_chart <- function(module_data, module_name) {
  # Reshape the data for the doughnut chart
  pie_data <- data.frame(
    Category = c("False Positive Chimeras", "Absolute Chimeras", "Borderline"),
    Value = c(module_data$False_positive_chimeras, module_data$Absolute_chimeras, module_data$Borderline)
  )
  
  # Calculate percentage for labels
  pie_data$Percentage <- pie_data$Value / sum(pie_data$Value) * 100
  pie_data$Label <- paste0(round(pie_data$Percentage, 1), "%")  # Format percentage labels
  
  # Create the doughnut chart
  doughnut_chart <- ggplot(pie_data, aes(x = 2, y = Value, fill = Category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +  # Create space in the center for the doughnut effect
    geom_text(aes(label = Label), position = position_stack(vjust = 0.5), color = "white", size = 4) +  # Add percentage labels
    theme_void() +  # Remove axes and gridlines
    theme(legend.position = "none") +  # Hide the legend since categories will be labeled outside the chart
    scale_fill_manual(values = c("#F8766D", "#00BFC4", "#7CAE00")) +  # Custom colors for the slices
    annotate("text", x = 0, y = 0, label = module_name, size = 6, fontface = "bold")  # Add module name in center
  
  return(doughnut_chart)
}

# Create a list to hold the individual plots
doughnut_plots <- list()

# Loop over each module and create corresponding doughnut chart
for (i in 1:nrow(chimera_data)) {
  doughnut_plots[[i]] <- create_doughnut_chart(chimera_data[i, ], chimera_data$Module[i])
}

# Combine all three doughnut charts horizontally using patchwork
combined_chart <- doughnut_plots[[1]] + doughnut_plots[[2]] + doughnut_plots[[3]] + 
  plot_layout(ncol = 3) + 
  plot_annotation(title = "Chimera Distribution Across Modules") & 
  theme(plot.title = element_text(hjust = 0.5))

# Add the category names at the bottom
combined_chart <- combined_chart + 
  plot_spacer() + 
  ggplot() + 
  annotate("text", x = 1:3, y = 0, label = c("False Positive Chimeras", "Absolute Chimeras", "Borderline"), size = 5) + 
  theme_void()

# Save the combined chart as an SVG
svg_filename <- "combined_doughnut_chart_nonchimeric.svg"
svglite(svg_filename, width = 13.33, height = 7.5)  # PowerPoint slide dimensions
print(combined_chart)
dev.off()




