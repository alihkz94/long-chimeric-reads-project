#The script designed for extracting the Phyla total abundance from the fasta files and the BLAST output of these seuqnces. (NOTE: for uchime_denovo some phyla needs to be removed from valid phyla)
import os
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from multiprocessing import Pool, cpu_count

# Define fungal phyla known in our dataset
fungal_phyla = {"Ascomycota", "Basidiomycota", "Glomeromycota", "Rozellomycota", "Chytridiomycota", "Mucoromycota"}

# Define a mapping for standardizing phylum names to match the plot
phylum_standardization = {
    "MAST-12": "Others (Eukaryotes)",
    "MAST-9": "Others (Eukaryotes)",
    # Add any other mappings as needed
}

# Adjust valid_phyla to only include major groups that should remain separate.
valid_phyla = {
    "Chrysophyta",
    "Cercozoa",
    "Gastrotricha",
    "Ascomycota",
    "Basidiomycota",
    "Bryophyta",
    "Ciliophora",
    "Chlorophyta",
    "Nematoda",
    "Arthropoda",
    "Dinoflagellata",
    "Radiolaria",
    "Cnidaria",
    "Katablepharidophyta",
    "Tracheophyta",
    "Annelida",
    "Bacillariophyta",
    "Others (Eukaryotes)",
    "Others (Fungi)"
}

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
                    if phylum in phylum_standardization:
                        return phylum_standardization[phylum], size
                    return phylum, size
    return 'Unknown', 0

def process_blast_file(filepath):
    print("Processing BLAST file:", filepath)
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
    
    if not blast_files:
        print("WARNING: No .txt files found in", blast_dir, ". Please check the path.")
        return defaultdict(list)
        
    with Pool(cpu_count()) as pool:
        results = pool.map(process_blast_file, blast_files)

    flattened_results = [item for sublist in results for item in sublist]
    
    df = pd.DataFrame(flattened_results, columns=["seq_id", "phylum", "size"])
    df.to_csv(output_file, index=False)
    print("BLAST results written to", output_file)
    
    blast_dict = defaultdict(list)
    for seq_id, phylum, size in flattened_results:
        blast_dict[seq_id].append((phylum, size))
    
    return blast_dict

def process_fasta_file(args):
    filename, blast_dict = args
    print("Processing FASTA file:", filename)
    phylum_counts = defaultdict(int)
    sample_name = os.path.basename(filename).replace("_rescued.fasta", "")
    
    try:
        for record in SeqIO.parse(filename, "fasta"):
            seq_id = record.id.split(";")[0]
            if seq_id in blast_dict:
                for phylum, size in blast_dict[seq_id]:
                    phylum_counts[(sample_name, phylum)] += size
    except Exception as e:
        print("Error processing", filename, ":", e)
    
    return phylum_counts

def count_phylum_abundance_parallel(fasta_dir, blast_dict):
    print("Starting parallel processing of FASTA files...")
    fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if f.endswith(".fasta")]
    
    if not fasta_files:
        print("WARNING: No .fasta files found in", fasta_dir, ". Please check the path.")
        return defaultdict(lambda: defaultdict(int))
    
    chunk_size = min(len(fasta_files), cpu_count())
    results = []
    
    for i in range(0, len(fasta_files), chunk_size):
        chunk = fasta_files[i:i+chunk_size]
        with Pool(len(chunk)) as pool:
            chunk_results = pool.map(process_fasta_file, [(f, blast_dict) for f in chunk])
            results.extend(chunk_results)
            print("Processed {}/{} FASTA files".format(i + len(chunk), len(fasta_files)))
    
    combined_counts = defaultdict(lambda: defaultdict(int))
    for result in results:
        for (sample, phylum), count in result.items():
            combined_counts[sample][phylum] += count
    print("FASTA files processed.")
    return combined_counts

def categorize_phylum_counts(phylum_counts, threshold_percentage=1):
    total_counts_by_phylum = defaultdict(int)
    for sample, phyla in phylum_counts.items():
        for phylum, count in phyla.items():
            total_counts_by_phylum[phylum] += count

    total_reads = sum(total_counts_by_phylum.values())
    phylum_percentages = {phylum: (count / total_reads) * 100 for phylum, count in total_counts_by_phylum.items()}
    
    total_phylum_counts = defaultdict(int)
    others_fungi_count = 0
    others_euk_count = 0
    
    categorized_as_others_fungi = set()
    categorized_as_others_euk = set()
    
    print("\nDetected phyla and their percentages:")
    for phylum, percentage in sorted(phylum_percentages.items(), key=lambda x: x[1], reverse=True):
        print("{}: {:.2f}%".format(phylum, percentage))
    
    for phylum, count in total_counts_by_phylum.items():
        percentage = phylum_percentages[phylum]
        # If the phylum is in valid_phyla (major groups) or its percentage exceeds the threshold, keep it separate.
        if phylum in valid_phyla or percentage >= threshold_percentage:
            total_phylum_counts[phylum] = count
        elif phylum in fungal_phyla:
            others_fungi_count += count
            categorized_as_others_fungi.add(phylum)
        else:
            others_euk_count += count
            categorized_as_others_euk.add(phylum)
    
    if others_fungi_count > 0:
        total_phylum_counts["Others (Fungi)"] = others_fungi_count

    if others_euk_count > 0:
        total_phylum_counts["Others (Eukaryotes)"] = others_euk_count

    print("\nPhyla categorized as 'Others (Fungi)':", ", ".join(categorized_as_others_fungi))
    print("Phyla categorized as 'Others (Eukaryotes)':", ", ".join(categorized_as_others_euk))
    
    return total_phylum_counts

def generate_total_abundance_table(phylum_counts, output_file):
    print("Generating total phylum abundance table...")
    sorted_phyla = sorted(phylum_counts.items(), key=lambda x: x[1], reverse=True)
    
    df = pd.DataFrame({
        "Phylum": [phylum for phylum, _ in sorted_phyla],
        "Total_Abundance": [count for _, count in sorted_phyla]
    })
    
    total_reads = df["Total_Abundance"].sum()
    df["Percentage"] = (df["Total_Abundance"] / total_reads * 100).round(2)
    
    df.to_csv(output_file, index=False)
    print("Total phylum abundance table saved to", output_file)
    
    print("\nPhylum Total Abundance Table:")
    print(df)
    
    return df

def main(fasta_dir, blast_dir, intermediate_file, output_table):
    print("Starting script execution...")
    print("FASTA directory:", fasta_dir)
    print("BLAST directory:", blast_dir)
    
    if not os.path.exists(fasta_dir):
        print("ERROR: FASTA directory", fasta_dir, "does not exist!")
        return
    
    if not os.path.exists(blast_dir):
        print("ERROR: BLAST directory", blast_dir, "does not exist!")
        return
    
    if os.path.exists(intermediate_file):
        print("Loading existing BLAST results from", intermediate_file)
        blast_df = pd.read_csv(intermediate_file)
        blast_dict = defaultdict(list)
        for _, row in blast_df.iterrows():
            phylum = row["phylum"]
            if phylum in phylum_standardization:
                phylum = phylum_standardization[phylum]
            blast_dict[row["seq_id"]].append((phylum, row["size"]))
    else:
        blast_dict = parse_blast_results_parallel(blast_dir, intermediate_file)
    
    phylum_counts = count_phylum_abundance_parallel(fasta_dir, blast_dict)
    
    # Use a threshold of 1% to group low-abundance phyla
    total_phylum_counts = categorize_phylum_counts(phylum_counts, threshold_percentage=1)
    
    generate_total_abundance_table(total_phylum_counts, output_table)
    
    print("Script completed successfully.")

if __name__ == '__main__':
    main("./rescued_reads", 
         "./blast_results",  
         "./blast_results_intermediate.csv", 
         "./total_phylum_abundance.csv")
