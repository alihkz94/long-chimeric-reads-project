## Create a table with formatted data and save as a high-resolution image ##

# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(grid)

# Create the data with better formatting and preserve spaces in column names
data <- data.frame(
  Method = c("UCHIME_denovo", "Chimeras_denovo", "DADA2"),
  `Total sequences` = c("41,087", "211,675", "320,032"),
  `False positive chimeras` = c("4,521 (11.00%)", "121,963 (57.62%)", "176,001 (54.99%)"),
  `Absolute chimeras` = c("36,566 (89.00%)", "89,712 (42.38%)", "144,031 (45.01%)"),
  stringsAsFactors = FALSE,  # Prevent unwanted factor levels
  check.names = FALSE         # Preserve spaces in column names
)

# Customize the theme for better appearance
theme_table <- ttheme_default(
  core = list(
    fg_params = list(hjust = 0.5, x = 0.5),
    bg_params = list(fill = c("white", "lightgrey")),
    fontface = "plain",
    fontsize = 14
  ),
  colhead = list(
    fg_params = list(hjust = 0.5, x = 0.5, fontface = "bold"),
    bg_params = list(fill = "grey80"),
    fontsize = 14
  )
)

# Generate the table
table_plot <- tableGrob(data, rows = NULL, theme = theme_table)

# Save the table as a high-resolution image (600 dpi)
png("chimera_table_fixed.png", width = 12, height = 4, units = 'in', res = 600)
grid.draw(table_plot)
dev.off()