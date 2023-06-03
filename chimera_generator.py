import os  # module for interacting with the operating system
import random  # module for generating random numbers
from Bio import SeqIO  # module for working with biological sequence data
from Bio.Seq import Seq  # class for representing biological sequences
from Bio.SeqRecord import SeqRecord  # class for representing a sequence record


def generate_chimeras(chimera_id_prefix="chimera"):
    # Get the current working directory
    input_directory = os.getcwd()

    # Create the output directory path
    output_directory = os.path.join(input_directory, "chimeric_reads")

    # Create the chimera information file path
    chimera_info_file = os.path.join(output_directory, "chimera_info.tsv")

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Get a list of input files in the input directory that end with ".fasta"
    input_files = [os.path.join(input_directory, file) for file in os.listdir(input_directory) if file.endswith(".fasta")]

    # Iterate over each selected input file
    for selected_input_file in input_files:
        records = []
        # Read all records from each input file and store them in a list
        for input_file in input_files:
            records.extend(list(SeqIO.parse(input_file, "fasta")))

        # Read the records from the selected input file
        main_records = list(SeqIO.parse(selected_input_file, "fasta"))

        # Set the mixed records to be the combined records from all input files
        mixed_records = records

        # Calculate the total number of reads in the selected input file
        total_reads = len(main_records)

        # Calculate the number of chimeric reads to generate (1-3% of total reads)
        num_chimeras = int(total_reads * random.uniform(0.01, 0.03))

        # Create an empty list to store the generated chimeras
        chimeras = []

        # Calculate the original abundance ratios of the main records
        original_ratios = calculate_abundance_ratio(main_records)

        # Open the chimera info file in append mode
        with open(chimera_info_file, "a") as chimera_info_handle:
            # Write the header to the chimera info file if it's empty
            if os.path.getsize(chimera_info_file) == 0:
                chimera_info_handle.write("chimera_id\tseq1_id\tseq2_id\tbreakpoints\treversed\tratio\tlength\n")

            i = 0
            # Generate the specified number of chimeras
            while i < num_chimeras:
                # Select two random sequences, either from the main records or mixed records
                if i < int(num_chimeras * 0.1):
                    seq1, seq2 = random.sample(main_records, 2)
                else:
                    seq1, seq2 = random.sample(mixed_records, 2)

                # Skip the iteration if seq1 is not present in the original ratios
                if seq1.id not in original_ratios:
                    continue

                # Generate a random breakpoint position between seq1 and seq2
                breakpoint = random.randint(1, len(seq1) - 1)

                # Create the chimera sequence by combining the sequences at the breakpoint
                chimera_seq = seq1.seq[:breakpoint] + seq2.seq[breakpoint:]

                # Reverse complement the chimera sequence for a certain percentage of chimeras
                if i < int(num_chimeras * 0.04):
                    chimera_seq = chimera_seq.reverse_complement()

                # Check the length of the chimera sequence to ensure it falls within a valid range
                if len(seq1) < 1.5 * len(chimera_seq) or len(seq1) > 10 * len(chimera_seq):
                    continue

                # Generate the chimera ID using the specified prefix and a unique number
                chimera_id = f"{chimera_id_prefix}_{i}"

                # Create a SeqRecord object for the chimera sequence
                chimera_record = SeqRecord(Seq(str(chimera_seq)), id=chimera_id, description="")

                # Add the chimera record to the list of chimeras
                chimeras.append(chimera_record)

                # Determine if the chimera is reversed or not based on the percentage
                reversed_status = "yes" if i < int(num_chimeras * 0.04) else "no"

                # Calculate the ratio of the chimera compared to the original sequence
                min_ratio = 0.1  # minimum value for the chimera ratio
                max_ratio = original_ratios[seq1.id] / 1.5  # maximum value for the chimera ratio
                ratio = random.uniform(min_ratio, min(max_ratio, 10))  # Select a ratio within the allowed range

                # Write the chimera information to the chimera info file
                chimera_info_handle.write(f"{chimera_id}\t{seq1.id}\t{seq2.id}\t{breakpoint}\t{reversed_status}\t{ratio}\t{len(chimera_seq)}\n")

                # Increment the counter
                i += 1

        # Create the output file path by joining the output directory and the selected input file name
        output_file = os.path.join(output_directory, os.path.basename(selected_input_file))

        # Insert the generated chimeras randomly into the main records
        for chimera in chimeras:
            position = random.randint(0, len(main_records) - 1)
            main_records.insert(position, chimera)

        # Write the modified main records (including chimeras) to the output file
        with open(output_file, "w") as output_handle:
            SeqIO.write(main_records, output_handle, "fasta")


def calculate_abundance_ratio(records):
    # Dictionary to store the abundance ratios
    abundance_ratios = {}

    # Calculate the total number of reads
    total_reads = len(records)

    # Count the occurrence of each sequence ID in the records
    for record in records:
        abundance_ratios[record.id] = abundance_ratios.get(record.id, 0) + 1

    # Calculate the abundance ratio for each sequence ID
    for seq_id, count in abundance_ratios.items():
        abundance_ratios[seq_id] = count / total_reads

    return abundance_ratios


if __name__ == "__main__":
    # Call the generate_chimeras function
    generate_chimeras()
