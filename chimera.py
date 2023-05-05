import os
import random
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# Function to mutate a sequence based on the mutation_rate
def mutate_sequence(seq_record, mutation_rate):
    sequence = str(seq_record.seq)
    mutated_sequence = []

    # Iterate through each base and mutate with the specified mutation_rate
    for base in sequence:
        if random.random() < mutation_rate:
            mutated_sequence.append(random.choice('ACTG'.replace(base, '')))
        else:
            mutated_sequence.append(base)

    # Create a new SeqRecord with the mutated sequence
    mutated_seq_string = "".join(mutated_sequence)
    mutated_seq_record = SeqRecord(Seq(mutated_seq_string, IUPAC.ambiguous_dna), id=seq_record.id, description=seq_record.description)

    return mutated_seq_record

# Function to generate chimeric reads
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

    with open(chimera_info_file, "w") as chimera_info_handle:
        chimera_info_handle.write("chimera_id\tseq1_id\tseq2_id\tbreakpoint\treversed\n")

        for i in range(num_chimeras):
            seq1, seq2 = random.sample(mixed_records, 2)
            breakpoint = random.randint(1, len(seq1) - 1)
            chimera_seq = seq1.seq[:breakpoint] + seq2.seq[breakpoint:]

            if i < int(num_chimeras * 0.04):
                chimera_seq = chimera_seq.reverse_complement()

            chimera_id = f"{chimera_id_prefix}_{i}"
            chimera_record = SeqRecord(Seq(str(chimera_seq), IUPAC.ambiguous_dna), id=chimera_id, description="")

            # Introduce mutations
            mutation_rate = 0.005  # Adjust this value based on the desired mutation rate
            chimera_record = mutate_sequence(chimera_record, mutation_rate)

            chimeras.append(chimera_record)

            reversed_status = "yes" if i < int(num_chimeras * 0.04) else "no"
            chimera_info_handle.write(f"{chimera_id}\t{seq1.id}\t{seq2.id}\t{breakpoint}\t{reversed_status}\n")

    # Insert chimeric reads at specified positions
    for chimera in chimeras:
        position = random.randint(0, len(main_records) - 1)
        main_records.insert(position, chimera)

    # Write the output file with chimeric reads
    with open(output_file, "w") as output_handle:
        SeqIO.write(main_records, output_handle, "fasta")

# Main function to execute the script
if __name__ == "__main__":
    if len(sys.argv) != 6:
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
