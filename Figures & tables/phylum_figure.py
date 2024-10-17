###########################################################
##Top 3 phyla abundance in FASTA files using BLAST results#
###########################################################
import os
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
from collections import defaultdict
from multiprocessing import Pool, cpu_count

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

# Function to process a single BLAST result file and write to intermediate file
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

# Parse all BLAST result files in parallel and write to intermediate file
def parse_blast_results_parallel(blast_dir, output_file):
    print("Starting parallel processing of BLAST results...")
    blast_files = [os.path.join(blast_dir, f) for f in os.listdir(blast_dir) if f.endswith(".txt")]
    with Pool(cpu_count()) as pool:
        results = pool.map(process_blast_file, blast_files)
    
    # Flatten the results list and write to CSV
    flattened_results = [item for sublist in results for item in sublist]
    df = pd.DataFrame(flattened_results, columns=['seq_id', 'phylum', 'size'])
    df.to_csv(output_file, index=False)
    print(f"BLAST results written to {output_file}")

# Function to map sequences in FASTA files to phyla using intermediate BLAST results file
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

# Parse FASTA files and count phylum occurrences with sizes in parallel
def count_phylum_abundance_parallel(fasta_dir, blast_results_file):
    print("Starting parallel processing of FASTA files...")
    blast_results_df = pd.read_csv(blast_results_file)
    fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if f.endswith(".fasta")]
    with Pool(cpu_count()) as pool:
        results = pool.map(process_fasta_file, [(f, blast_results_df) for f in fasta_files])

    # Combine results from all files
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

# Create a stacked bar plot for top N phylum abundance
def plot_phylum_abundance(phylum_counts, output_file, top_n=3):
    print("Creating phylum abundance plot...")
    top_phyla_counts = get_top_phyla(phylum_counts, top_n)
    df = pd.DataFrame(top_phyla_counts).fillna(0).T
    df.sort_index(axis=0, inplace=True)

    # Define distinct colors for up to 30 phyla
    colors = [
        "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",  # Red, Green, Yellow, Blue, Orange
        "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",  # Purple, Cyan, Magenta, Lime, Soft Pink
        "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000",  # Teal, Lavender, Brown, Pale Yellow, Maroon
        "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080",  # Mint Green, Olive, Peach, Navy, Gray
        "#67001f", "#053061", "#fddbc7", "#c7eae5", "#ca0020",  # Burgundy, Dark Blue, Soft Coral, Light Teal, Dark Red
        "#0571b0", "#f4a582", "#92c5de", "#4d4d4d", "#d01c8b"   # Medium Blue, Salmon, Sky Blue, Dark Gray, Bright Pink
]
    
    # Plot
    df.plot(kind='bar', stacked=True, figsize=(15, 10), color=colors[:len(df.columns)])
    plt.title(f'Top {top_n} Phylum Abundance in FASTA Files')
    plt.xlabel('Sample')
    plt.ylabel('Abundance')
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='Phylum', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_file, dpi = 600)
    plt.show()
    print(f"Phylum abundance plot saved to {output_file}")

# Main function to execute the workflow
def main(fasta_dir, blast_dir, intermediate_file, output_file, top_n=3):
    print("Starting script execution...")
    parse_blast_results_parallel(blast_dir, intermediate_file)
    phylum_counts = count_phylum_abundance_parallel(fasta_dir, intermediate_file)
    plot_phylum_abundance(phylum_counts, output_file, top_n)
    print("Script completed.")

# Define input directories and output file
fasta_dir = './rescued_reads'
blast_dir = './blast_output'
output_file = 'top_phylum_abundance.png'
intermediate_file = './blast_results_intermediate.csv'

# Run the main function
if __name__ == '__main__':
    main(fasta_dir, blast_dir, intermediate_file, output_file, top_n=3)
