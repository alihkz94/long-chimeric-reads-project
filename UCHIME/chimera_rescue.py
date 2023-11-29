import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import Counter

def rescue_chimeric_sequences(input_dir, min_occurrence, output_dir):
    aggregated_sequences = []
    sequence_occurrences = Counter()

    # Read and aggregate sequences
    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            file_path = os.path.join(input_dir, filename)
            with open(file_path, 'r') as fasta_file:
                seqs_in_file = list(SeqIO.parse(fasta_file, "fasta"))
                print(f"Found {len(seqs_in_file)} sequences in {filename}")

                for record in seqs_in_file:
                    sequence = str(record.seq)
                    aggregated_sequences.append(sequence)
                    sequence_occurrences[sequence] += 1

    print(f"Total unique sequences aggregated: {len(set(aggregated_sequences))}")

    # Filter sequences
    rescued_sequences = {seq for seq, count in sequence_occurrences.items() if count >= min_occurrence}
    print(f"Number of sequences above occurrence threshold: {len(rescued_sequences)}")

    # Write output
    if rescued_sequences:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        output_file = os.path.join(output_dir, "RescuedChimera.fa")
        with open(output_file, "w") as output_handle:
            for seq in rescued_sequences:
                seq_record = SeqRecord(SeqIO.Seq(seq), id="RescuedChimera", description="")
                SeqIO.write(seq_record, output_handle, "fasta")

        print(f"Rescued {len(rescued_sequences)} unique chimeric sequences to {output_file}")
    else:
        print("No sequences met the rescue criteria.")

if __name__ == "__main__":
    input_directory = "/home/ali/Documents/simulated_data/analysis/Blast/treat_chimera/chimeras/combined_chimeras" # Path to input directory containing FASTA files
    min_occurrence = 2  # Minimum occurrence threshold
    output_directory = "/home/ali/Documents/simulated_data/analysis/Blast/treat_chimera/chimeras/combined_chimeras/rescued_chimeras" # Path to output directory

    try:
        rescue_chimeric_sequences(input_directory, min_occurrence, output_directory)
    except Exception as e:
        print(f"An error occurred: {e}")
