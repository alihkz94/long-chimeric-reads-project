library(seqinr)
library(ggplot2)

setwd("~/Documents/simulated_data/analysis/DADA2/test_settings/fasta_files")

##### TWO FILE COMPARISON ####

# Function to read and process FASTA files
process_fasta <- function(fasta_file) {
  fasta = read.fasta(fasta_file, seqtype = "DNA", as.string = TRUE, seqonly = FALSE)
  seq_names = getName(fasta)
  seqs = unlist(getSequence(fasta, as.string = TRUE))
  break_pts = gsub("_.*", "", (gsub(".*_at_", "", seq_names)))
  lens = nchar(seqs)
  
  scores = numeric(length(lens))
  
  for (i in 1:length(lens)) {
    part_a = as.numeric(lens[i]) - as.numeric(break_pts[i])
    part_b = as.numeric(lens[i]) - part_a
    
    if (part_a <= part_b) {
      scores[i] = round(part_a / as.numeric(lens[i]), digits = 4)
    } else {
      scores[i] = round(part_b / as.numeric(lens[i]), digits = 4)
    }
  }
  
  return(data.frame(seq_names, scores))
}

# File name mapping based on maxshift values for labeling
file_label_mapping <- c(
  "11" = "maxshift_100",
  "12" = "maxshift_150",
  "13" = "maxshift_200",
  "14" = "maxshift_250",
  "15" = "maxshift_300",
  "16" = "maxshift_325",
  "17" = "maxshift_350"
)

# Original file names for reading from the directory
file_name_mapping <- c(
  "11" = "removed_chimeric_sequences_11_chimeric.fasta",
  "12" = "removed_chimeric_sequences_12_chimeric.fasta",
  "13" = "removed_chimeric_sequences_13_chimeric.fasta",
  "14" = "removed_chimeric_sequences_14_chimeric.fasta",
  "15" = "removed_chimeric_sequences_15_chimeric.fasta",
  "16" = "removed_chimeric_sequences_16_chimeric.fasta",
  "17" = "removed_chimeric_sequences_17_chimeric.fasta"
)

# Process both FASTA files
df1 = process_fasta(file_name_mapping["11"])
df2 = process_fasta(file_name_mapping["17"])

# Identify unique sequences in both files
unique_to_df1 = setdiff(df1$seq_names, df2$seq_names)
unique_to_df2 = setdiff(df2$seq_names, df1$seq_names)

# Filter data frames to only include unique sequences
unique_df1 = df1[df1$seq_names %in% unique_to_df1,]
unique_df2 = df2[df2$seq_names %in% unique_to_df2,]

# Visualization
all_data = rbind(
  data.frame(Score = unique_df1$scores, Type = paste0("Unique to ", file_label_mapping["11"])),
  data.frame(Score = unique_df2$scores, Type = paste0("Unique to ", file_label_mapping["17"])),
  data.frame(Score = df1$scores, Type = paste0("All in ", file_label_mapping["11"])),
  data.frame(Score = df2$scores, Type = paste0("All in ", file_label_mapping["17"]))
)

ggplot(all_data, aes(x = Type, y = Score)) + geom_boxplot() + ggtitle("Breakpoint Scores Comparison")

# Enhanced ggplot2 visualization
enhanced_plot <- ggplot(all_data, aes(x = Type, y = Score, fill = Type)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c(
    `Unique to maxshift_100` = "red", 
    `Unique to maxshift_350` = "blue", 
    `All in maxshift_100` = "green", 
    `All in maxshift_350` = "purple"
  )) +
  labs(title = "Comparison of Breakpoint Scores Across Files",
       subtitle = "Breakpoint scores of unique and common chimeras",
       x = "File Type",
       y = "Breakpoint Score",
       fill = "Type") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12)
  )

# Display the enhanced plot
print(enhanced_plot)

##### ALL FILES ####

# Initialize an empty data frame for all chimeras and unique chimeras
all_data = data.frame()
unique_data = data.frame()

# Process the original file (number 11)
df1 = process_fasta(file_name_mapping["11"])

# Add unique sequences for the first file (maxshift_100) to the unique_data data frame
unique_to_df1_all = df1$seq_names
for(i in 12:17) {
  df_current = process_fasta(file_name_mapping[as.character(i)])
  unique_to_df1_all = setdiff(unique_to_df1_all, df_current$seq_names)
}
unique_df1_all = df1[df1$seq_names %in% unique_to_df1_all,]
unique_data = rbind(unique_data, data.frame(Score = unique_df1_all$scores, Type = "Unique to maxshift_100"))

# Loop through files 12 to 17 to find unique sequences compared to file 11
for(i in 12:17) {
  file_name = file_name_mapping[as.character(i)]
  df_current = process_fasta(file_name)
  
  all_data = rbind(all_data, data.frame(Score = df_current$scores, Type = paste0("All in ", file_label_mapping[as.character(i)])))
  
  unique_to_current = setdiff(df_current$seq_names, df1$seq_names)
  unique_df_current = df_current[df_current$seq_names %in% unique_to_current,]
  
  unique_data = rbind(unique_data, data.frame(Score = unique_df_current$scores, Type = paste0("Unique to ", file_label_mapping[as.character(i)])))
}

# Visualization for all chimeras with new x-axis labels
all_plot = ggplot(all_data, aes(x = Type, y = Score, fill = Type)) +
  geom_boxplot(alpha = 0.7) +
  scale_x_discrete(labels = file_label_mapping) +
  labs(title = "Breakpoint Scores for All Chimeras",
       x = "Maxshift Conditions",
       y = "Breakpoint Score") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold"))
print(all_plot)

# Visualization for unique chimeras with new x-axis labels
unique_plot = ggplot(unique_data, aes(x = Type, y = Score, fill = Type)) +
  geom_boxplot(alpha = 0.7) +
  scale_x_discrete(labels = file_label_mapping) +
  labs(title = "Breakpoint Scores for Unique Chimeras",
       x = "Maxshift Conditions",
       y = "Breakpoint Score") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold"))
print(unique_plot)


# Save the enhanced plot as a PDF with increased width for clarity
ggsave("enhanced_plot.pdf", plot = enhanced_plot, width = 12, height = 6)

# Save the all chimeras plot as a PDF with increased width for clarity
ggsave("all_chimeras_plot.pdf", plot = all_plot, width = 12, height = 6)

# Save the unique chimeras plot as a PDF with increased width for clarity
ggsave("unique_chimeras_plot.pdf", plot = unique_plot, width = 16, height = 6)
