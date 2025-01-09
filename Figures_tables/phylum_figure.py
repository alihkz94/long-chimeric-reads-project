###########################################################
##Top 3 phyla abundance in FASTA files using BLAST results#
###########################################################
import os
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
from collections import defaultdict
from multiprocessing import Pool, cpu_count

# Fixed color map for phyla
fixed_color_map = {
    "Basidiomycota": "#1f77b4",  # Blue
    "Ascomycota": "#ff7f0e",     # Orange
    "Glomeromycota": "#2ca02c",  # Green
    "Gastrotricha": "#d62728",   # Red
    "Chrysophyta": "#9467bd",    # Purple
    "Arthropoda": "#8c564b",     # Brown
    "Cnidaria": "#e377c2",       # Pink
    "Diplonemia": "#7f7f7f",     # Gray
    "Radiolaria": "#bcbd22",     # Olive
    "Dinoflagellata": "#17becf", # Teal
    "Chlorophyta": "#ffbb78",    # Light Orange
    "Pelagophyta": "#98df8a",    # Light Green
    "Ciliophora": "#c5b0d5",     # Light Purple
    "Rozellomycota": "#ff9896",  # Light Red
    "Hyphochytriomycota": "#f7b6d2", # Light Pink
    "Dictyochophyta": "#9edae5", # Light Teal
    "Cercozoa": "#aec7e8",       # Light Blue
    "Nematoda": "#c49c94",       # Dark Brown
    "Evosea": "#dbdb8d",         # Light Olive
    "Annelida": "#17becf",       # Dark Teal
    "Perkinsia": "#ffbb78",      # Yellow-Orange
    "Katablepharidophyta": "#c7c7c7", # Grayish
    "MAST-9": "#98df8a",         # Mint Green
    "Bryophyta": "#ff9896",      # Coral
    "Bacillariophyta": "#c5b0d5",# Lavender
    "Tracheophyta": "#d62728",   # Dark Red
    "Petalomonadia": "#ff7f0e",  # Dark Orange
    "Kinetoplastia": "#7f7f7f",  # Gray
}

# Function to extract phylum and size from the BLAST result line
def extract_phylum_and_size(blast_line):
    parts = blast_line.split('+')
    if len(parts) > 1:
        taxonomy = parts[1].split(';')
        size_part = parts[0].split(';size=')[-1]
        size = int(size_part.split('+')[0])
        for taxon in taxonomy:
            if taxon.startswith('p__'):
                phylum = taxon.split('__')[1]
                if "_phy_Incertae_sedis" not in phylum:
                    return phylum, size
    return 'Unknown', 0

# Function to process a single BLAST result file
def process_blast_file(filepath):
    print(f"Processing BLAST file: {filepath}")
    results = []
    with open(filepath) as f:
        for line in f:
            if line.strip() and '+' in line:
                seq_id = line.split(';')[0]
                phylum, size = extract_phylum_and_size(line)
                if phylum != 'Unknown':
                    results.append((seq_id, phylum, size))
    return results

# Parse BLAST result files in parallel
def parse_blast_results_parallel(blast_dir, output_file):
    print("Starting parallel processing of BLAST results...")
    blast_files = [os.path.join(blast_dir, f) for f in os.listdir(blast_dir) if f.endswith(".txt")]
    with Pool(cpu_count()) as pool:
        results = pool.map(process_blast_file, blast_files)

    flattened_results = [item for sublist in results for item in sublist]
    df = pd.DataFrame(flattened_results, columns=['seq_id', 'phylum', 'size'])
    df.to_csv(output_file, index=False)
    print(f"BLAST results written to {output_file}")

# Function to map sequences to phyla in FASTA files
def process_fasta_file(args):
    filename, blast_results_df = args
    print(f"Processing FASTA file: {filename}")
    phylum_counts = defaultdict(int)
    sample_name = os.path.basename(filename).replace('_rescued.fasta', '')
    for record in SeqIO.parse(filename, "fasta"):
        seq_id = record.id.split(';')[0]
        phylum_data = blast_results_df[blast_results_df['seq_id'] == seq_id]
        for _, row in phylum_data.iterrows():
            phylum_counts[(sample_name, row['phylum'])] += row['size']
    return phylum_counts

# Parse FASTA files and count phylum occurrences
def count_phylum_abundance_parallel(fasta_dir, blast_results_file):
    print("Starting parallel processing of FASTA files...")
    blast_results_df = pd.read_csv(blast_results_file)
    fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if f.endswith(".fasta")]
    with Pool(cpu_count()) as pool:
        results = pool.map(process_fasta_file, [(f, blast_results_df) for f in fasta_files])

    combined_counts = defaultdict(lambda: defaultdict(int))
    for result in results:
        for (sample, phylum), count in result.items():
            combined_counts[sample][phylum] += count
    print("FASTA files processed.")
    return combined_counts

# Function to filter and get top N phyla per sample
def get_top_phyla(phylum_counts, top_n=3):
    top_phyla_counts = defaultdict(dict)
    for sample, phyla in phylum_counts.items():
        sorted_phyla = sorted(phyla.items(), key=lambda item: item[1], reverse=True)
        for phylum, count in sorted_phyla[:top_n]:
            top_phyla_counts[sample][phylum] = count
    return top_phyla_counts

# Normalize phylum counts to percentages
def normalize_phylum_counts(phylum_counts):
    normalized_counts = defaultdict(dict)
    for sample, phyla in phylum_counts.items():
        total_reads = sum(phyla.values())
        for phylum, count in phyla.items():
            normalized_counts[sample][phylum] = (count / total_reads) * 100  # Convert to percentage
    return normalized_counts

# Create a stacked percentage bar plot for top N phyla
def plot_percentage_phylum_abundance(phylum_counts, output_file, top_n=3):
    print("Creating percentage-based phylum abundance plot...")
    top_phyla_counts = get_top_phyla(phylum_counts, top_n)
    df = pd.DataFrame(top_phyla_counts).fillna(0).T
    df = df.div(df.sum(axis=1), axis=0) * 100  # Convert to percentage
    df.sort_index(axis=0, inplace=True)

    # Assign colors based on the fixed color map
    colors = [fixed_color_map.get(phylum, '#333333') for phylum in df.columns]  # Default color for phyla not in the map

    # Plot the stacked percentage bar plot with larger, bold fonts
    fig, ax = plt.subplots(figsize=(15, 10))  # Adjust size for PowerPoint
    df.plot(kind='bar', stacked=True, color=colors, ax=ax)
    
    # Set larger, bold fonts
    ax.set_xlabel('Sample', fontsize=16, fontweight='bold')
    ax.set_ylabel('Relative Abundance (%)', fontsize=16, fontweight='bold')
    # Set larger, bold tick labels
    ax.tick_params(axis='x', labelsize=14, labelrotation=90, labelcolor='black', width=2)
    ax.tick_params(axis='y', labelsize=14, labelcolor='black', width=2)
    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Rotate x-axis labels and align them with the bars
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center')

    # Customize legend
    legend = ax.legend(title='Phylum', title_fontsize=18, fontsize=18, bbox_to_anchor=(1.05, 1), loc='upper left')
    legend.get_title().set_fontweight('bold')
    for text in legend.get_texts():
        text.set_fontweight('bold')

    plt.tight_layout()
    
    # Save the figure in SVG format for scalability in PowerPoint
    plt.savefig(output_file, format='svg')
    plt.show()
    print(f"Percentage-based phylum abundance plot saved to {output_file}")

# Main function to execute the workflow
def main(fasta_dir, blast_dir, intermediate_file, output_file, top_n=3):
    print("Starting script execution...")
    parse_blast_results_parallel(blast_dir, intermediate_file)
    phylum_counts = count_phylum_abundance_parallel(fasta_dir, intermediate_file)
    normalized_counts = normalize_phylum_counts(phylum_counts)
    plot_percentage_phylum_abundance(normalized_counts, output_file, top_n)
    print("Script completed.")

# Run the main function
if __name__ == '__main__':
    main('./rescued_reads', './blast_output', './blast_results_intermediate.csv', 'percentage_phylum_abundance.svg', top_n=3)
