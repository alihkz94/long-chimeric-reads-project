import os
import random
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def generate_chimeras(selected_input_file, input_files, output_file, chimera_info_file, chimera_id_prefix="chimera"):
    records = []
    for input_file in input_files:
        records.extend(list(SeqIO.parse(input_file, "fasta")))

    main_records = list(SeqIO.parse(selected_input_file, "fasta"))
    mixed_records = records

    # Calculate the number of chimeras based on the input file's total reads
    total_reads = len(main_records)
    num_chimeras = int(total_reads * random.uniform(0.01, 0.03))  # 1-3% chimeric reads

    chimeras = []

    # Calculate the abundance ratios of the original sequences in the fasta file
    original_ratios = calculate_abundance_ratio(main_records)

    with open(chimera_info_file, "w") as chimera_info_handle:
        chimera_info_handle.write("chimera_id\tseq1_id\tseq2_id\tbreakpoint\treversed\tratio\tlength\n")

        i = 0
        while i < num_chimeras:
            if i < int(num_chimeras * 0.1):
                seq1, seq2 = random.sample(main_records, 2)
            else:
                seq1, seq2 = random.sample(mixed_records, 2)

            # If seq1.id is not present in the original_ratios dictionary, skip this iteration
            if seq1.id not in original_ratios:
                continue

            breakpoint = random.randint(1, len(seq1) - 1)
            chimera_seq = seq1.seq[:breakpoint] + seq2.seq[breakpoint:]

            if i < int(num_chimeras * 0.04):
                chimera_seq = chimera_seq.reverse_complement()

            chimera_id = f"{chimera_id_prefix}_{i}"
            chimera_record = SeqRecord(Seq(str(chimera_seq)), id=chimera_id, description="")

            chimeras.append(chimera_record)

            reversed_status = "yes" if i < int(num_chimeras * 0.04) else "no"

            # Calculate the chimera ratio such that the parent ratio is 1.5 to 10 times greater
            min_ratio = 0.1  # minimum value for the chimera ratio
            max_ratio = min(original_ratios[seq1.id] * 10, 1)  # maximum value for the chimera ratio
            ratio = random.uniform(min_ratio, max_ratio)  # Select a ratio within the allowed range

            chimera_info_handle.write(f"{chimera_id}\t{seq1.id}\t{seq2.id}\t{breakpoint}\t{reversed_status}\t{ratio}\t{len(chimera_seq)}\n")

            i += 1

    # Insert chimeric reads at specified positions
    for chimera in chimeras:
        position = random.randint(0, len(main_records) - 1)
        main_records.insert(position, chimera)

    # Write the output file with chimeric reads
    output_directory = os.path.dirname(output_file)
    os.makedirs(output_directory, exist_ok=True)
    with open(output_file, "w") as output_handle:
        SeqIO.write(main_records, output_handle, "fasta")

def calculate_abundance_ratio(records):
    abundance_ratios = {}
    total_reads = len(records)

    for record in records:
        abundance_ratios[record.id] = abundance_ratios.get(record.id, 0) + 1

    for seq_id, count in abundance_ratios.items():
        abundance_ratios[seq_id] = count / total_reads

    return abundance_ratios

# Main function to execute the script
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} selected_input_file input_directory output_fasta chimera_info_file")
        sys.exit(1)
    selected_input_file = sys.argv[1]
    input_directory = sys.argv[2]
    output_fasta = sys.argv[3]
    chimera_info_file = sys.argv[4]

    # Get input_files in the input_directory with a .fasta extension
    input_files = [os.path.join(input_directory, file) for file in os.listdir(input_directory) if file.endswith(".fasta")]

    # Generate chimeric reads and write them to the output file
    generate_chimeras(selected_input_file, input_files, output_fasta, chimera_info_file)
