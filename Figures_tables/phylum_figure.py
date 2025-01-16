###############################################################
## Phyla relative abundance in FASTA files using BLAST results#
###############################################################
import os
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
from collections import defaultdict
from multiprocessing import Pool, cpu_count

# Biome mapping
biome_mapping = {
    'ERR6454461': 'freshwater', 'ERR6454462': 'freshwater',
    'ERR6454463': 'freshwater', 'ERR6454464': 'freshwater',
    'ERR6454465': 'Soil', 'ERR6454466': 'Soil',
    'ERR6454467': 'Soil', 'ERR6454468': 'Soil',
    'ERR6454469': 'Soil', 'ERR6454470': 'marine',
    'ERR6454471': 'marine', 'ERR6454472': 'marine',
    'ERR6454473': 'marine', 'ERR6454474': 'marine',
    'ERR6454475': 'marine', 'ERR6454476': 'marine',
    'ERR6454477': 'marine', 'ERR6454478': 'marine'
}

# Fixed color map for phyla, with custom colors for "Others (Fungi)" and "Others (Eukaryotes)"
fixed_color_map = {
    "Arthropoda": "#1f77b4",         # Blue 
    "Dinoflagellata": "#ff7f0e",     # Orange 
    "Chrysophyta": "#2ca02c",        # Green 
    "Gastrotricha": "#d62728",       # Red 
    "Bryophyta": "#7a44a3",          # Sharper Purple
    "Ciliophora": "#8c564b",         # Brown 
    "Bacillariophyta": "#d63a94",    # Sharper Pink
    "Radiolaria": "#a1a30d",         # Sharper Olive
    "Nematoda": "#128593",           # Sharper Teal
    "Basidiomycota": "#5e5e5e",      # Darker Gray
    "Ascomycota": "#ffbb78",         # Sharper Light Orange
    "Tracheophyta": "#6bbc47",       # Sharper Light Green
    "Cnidaria": "#9b71c5",           # Sharper Light Purple
    "Chlorophyta": "#ff7b7b",        # Sharper Light Coral
    "Annelida": "#47d2d8",           # Sharper Light Teal
    "Cercozoa": "#8b5c4a",           # Sharper Dark Brown
    "Katablepharidophyta": "#b3b34d",# Sharper Light Olive
    "Dictyochophyta": "#9c9c9c",     # Sharper Light Gray
    "Rozellomycota": "#e469a3",      # Sharper Light Pink
    "Diplonemia": "#10676f",         # Unique Deep Teal
    "Pelagophyta": "#4cc4cc",        # Unique Bright Cyan
    "Glomeromycota": "#e75c5c"       # Unique Coral
}

# Define fungal phyla known in our dataset
fungal_phyla = {"Ascomycota", "Basidiomycota", "Glomeromycota", "Rozellomycota"}

# Assign distinct colors for the new "Others" categories
others_fungi_color = "#FFD700"       # Gold
others_eukaryotes_color = "#333333"  # Black
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

def parse_blast_results_parallel(blast_dir, output_file):
    print("Starting parallel processing of BLAST results...")
    blast_files = [os.path.join(blast_dir, f) for f in os.listdir(blast_dir) if f.endswith(".txt")]
    with Pool(cpu_count()) as pool:
        results = pool.map(process_blast_file, blast_files)

    flattened_results = [item for sublist in results for item in sublist]
    df = pd.DataFrame(flattened_results, columns=['seq_id', 'phylum', 'size'])
    df.to_csv(output_file, index=False)
    print(f"BLAST results written to {output_file}")

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

def normalize_and_categorize_phylum_counts(phylum_counts):
    normalized_counts = defaultdict(dict)
    for sample, phyla in phylum_counts.items():
        total_reads = sum(phyla.values())
        others_fungi_total = 0.0
        others_euk_total = 0.0
        for phylum, count in phyla.items():
            percentage = (count / total_reads) * 100
            # Check if phylum is not in the fixed_color_map or under 10%
            if percentage < 10 or phylum not in fixed_color_map or phylum == "MAST-9":
                if phylum in fungal_phyla:
                    others_fungi_total += percentage
                else:
                    others_euk_total += percentage
            else:
                normalized_counts[sample][phylum] = percentage

        # Add "Others (Fungi)" if any fungal low-abundance phyla exist
        if others_fungi_total > 0:
            normalized_counts[sample]["Others (Fungi)"] = others_fungi_total

        # Add "Others (Eukaryotes)" if any non-fungal low-abundance phyla exist
        if others_euk_total > 0:
            normalized_counts[sample]["Others (Eukaryotes)"] = others_euk_total

    return normalized_counts

def plot_percentage_phylum_abundance_by_biome(phylum_counts, biome_mapping, output_file):
    print("Creating percentage-based phylum abundance plot separated by biomes...")

    # Convert dictionary to DataFrame
    df = pd.DataFrame(phylum_counts).fillna(0).T
    df['Sample'] = df.index
    df['Biome'] = df['Sample'].map(biome_mapping)

    # Define custom order for biomes
    biome_order = ['freshwater', 'Soil', 'marine']
    df['Biome'] = pd.Categorical(df['Biome'], categories=biome_order, ordered=True)

    # Sort samples by biome and sample name
    df.sort_values(by=['Biome', 'Sample'], inplace=True)

    # Check if DataFrame is empty after mapping and sorting
    if df.empty:
        print("No data available for plotting. Ensure the input data is correct.")
        return

    # Get all phyla columns
    phyla_columns = [col for col in df.columns if col not in ['Sample', 'Biome']]

    # Assign colors
    color_map = []
    for phylum in phyla_columns:
        if phylum in fixed_color_map:
            color_map.append(fixed_color_map[phylum])
        elif phylum == "Others (Fungi)":
            color_map.append(others_fungi_color)
        elif phylum == "Others (Eukaryotes)":
            color_map.append(others_eukaryotes_color)
        else:
            color_map.append("#000000")

    # Unique biomes and their relative sample counts
    unique_biomes = df['Biome'].cat.categories
    sample_counts = [len(df[df['Biome'] == biome]) for biome in unique_biomes]

    # Safeguard against empty sample counts or zero total samples
    total_samples = sum(sample_counts)
    if total_samples == 0:
        print("Total samples count is zero. Check the biome mapping or input data for inconsistencies.")
        return

    # Adjust width ratios based on sample counts
    width_ratios = [count / total_samples * len(sample_counts) for count in sample_counts]

    # Create subplots with adjusted widths
    fig, axes = plt.subplots(
        1, len(unique_biomes),
        figsize=(5 * len(unique_biomes) + 5, 10),
        gridspec_kw={'width_ratios': width_ratios},
        sharey=True
    )

    # Ensure axes is iterable
    if len(unique_biomes) == 1:
        axes = [axes]

    # Plot each biome
    for ax, biome in zip(axes, unique_biomes):
        biome_df = df[df['Biome'] == biome]
        if biome_df.empty:
            print(f"No data for biome: {biome}. Skipping plot for this biome.")
            continue

        biome_df.set_index('Sample', inplace=False)[phyla_columns].plot(
            kind='bar', stacked=True, color=color_map, ax=ax, legend=False
        )
        ax.set_title(f"{biome.capitalize()} Biome", fontsize=18, fontweight='bold')
        ax.set_xlabel('Sample', fontsize=14, fontweight='bold')
        ax.set_ylabel('Relative Abundance (%)', fontsize=18, fontweight='bold')
        ax.tick_params(axis='x', labelsize=14, labelrotation=90, labelcolor='black', width=2)
        ax.tick_params(axis='y', labelsize=14, labelcolor='black', width=2)
        for label in ax.get_xticklabels():
            label.set_fontweight('bold')
        for label in ax.get_yticklabels():
            label.set_fontweight('bold')

    # Unified legend on the right
    handles, labels = axes[0].get_legend_handles_labels()
    fig.subplots_adjust(right=0.85)  # Space for legend
    legend = fig.legend(handles, labels, title='Phylum', title_fontsize=18, fontsize=14,
                        loc='center right', bbox_to_anchor=(1.0, 0.5))
    legend.get_title().set_fontweight('bold')
    for text in legend.get_texts():
        text.set_fontweight('bold')

    # Adjust layout
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(output_file, format='svg')
    plt.show()
    print(f"Percentage-based phylum abundance plot by biome saved to {output_file}")

def main(fasta_dir, blast_dir, intermediate_file, output_file):
    print("Starting script execution...")
    parse_blast_results_parallel(blast_dir, intermediate_file)
    phylum_counts = count_phylum_abundance_parallel(fasta_dir, intermediate_file)
    normalized_counts = normalize_and_categorize_phylum_counts(phylum_counts)
    # Plot by biome
    plot_percentage_phylum_abundance_by_biome(normalized_counts, biome_mapping, output_file)
    print("Script completed.")

if __name__ == '__main__':
    main('./rescued_reads', './blast_results', './blast_results_intermediate.csv', 'percentage_phylum_abundance_by_biome.svg')
