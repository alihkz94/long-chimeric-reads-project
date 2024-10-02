import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Sample data for custom methods
methods_custom = ['UCHIME_denovo custom', 'Chimeras_denovo custom', 'DADA2 custom']
uchime_custom_coverage = [38.12, 56.90, 27.31, 34.54, 53.00, 41.75, 51.95, 67.88, 90.06, 53.35, 39.78, 67.97, 68.72, 45.51, 35.45, 41.26, 28.26, 59.16]
chimeras_custom_coverage = [36.65, 53.37, 22.99, 32.83, 51.45, 40.40, 51.37, 66.68, 88.79, 58.74, 38.66, 63.73, 69.17, 41.37, 31.51, 36.20, 29.20, 52.62]
dada2_custom_coverage = [35.50, 50.98, 20.95, 31.38, 49.66, 39.39, 50.79, 65.87, 87.91, 55.65, 37.56, 59.46, 67.68, 39.37, 30.12, 33.58, 29.03, 49.69]

# Data from "uchime default.txt"
uchime_default_coverage = [38.01, 56.79, 27.24, 34.61, 52.93, 41.60, 51.99, 67.86, 89.95, 53.71, 40.08, 67.91, 70.00, 45.51, 35.40, 41.02, 28.39, 58.84]

# Data from "chimeras denovo default.txt"
chimeras_default_coverage = [35.08, 52.20, 22.70, 32.60, 49.25, 39.51, 49.93, 65.52, 87.70, 55.70, 37.50, 62.31, 67.99, 40.27, 32.14, 34.71, 28.54, 51.40]

# Combine data into a DataFrame
methods = ['UCHIME_denovo custom']*len(uchime_custom_coverage) + \
          ['Chimeras_denovo custom']*len(chimeras_custom_coverage) + \
          ['DADA2']*len(dada2_custom_coverage) + \
          ['UCHIME_denovo default']*len(uchime_default_coverage) + \
          ['Chimeras_denovo default']*len(chimeras_default_coverage)

coverage = uchime_custom_coverage + chimeras_custom_coverage + dada2_custom_coverage + \
           uchime_default_coverage + chimeras_default_coverage

data = {'Method': methods, 'Coverage': coverage}
df = pd.DataFrame(data)

# Create the plot with larger fonts and increased label padding
plt.figure(figsize=(12, 8))
sns.boxplot(x='Method', y='Coverage', data=df, palette="Set3")
sns.stripplot(x='Method', y='Coverage', data=df, color='black', alpha=0.5, jitter=True)

# Add titles and labels with larger font sizes and increased padding
plt.xlabel('Chimeras Detection Methods', fontsize=16, labelpad=20)
plt.ylabel('Query Coverage Percentage', fontsize=16, labelpad=20)
plt.title('Query Coverage of Reads Across Different Chimeras Detection Methods', fontsize=20, pad=30)
plt.xticks(fontsize=14, rotation=15)
plt.yticks(fontsize=14)
plt.grid(True)

# Save the figure with high resolution
plt.savefig('Query_Coverage_Boxplot_Updated.png', dpi=300, bbox_inches='tight')
plt.close()

