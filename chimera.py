# Import necessary libraries
import random
import sys
from Bio import SeqIO

# Define the generate_chimeras function with appropriate parameters
def generate_chimeras(fastq_file, output_file, chimera_info_file, num_chimeras, chimera_id_prefix="chimera"):
    # Read the input FASTQ file and store the records in a list
    records = list(SeqIO.parse(fastq_file, "fastq"))
    # Initialize an empty list to store chimeras
    chimeras = []

    # Open the chimera_info_file for writing
    with open(chimera_info_file, "w") as chimera_info_handle:
        # Write the header for the chimera info file
        chimera_info_handle.write("chimera_id\tseq1_id\tseq2_id\tbreakpoint\n")

        # Iterate for the desired number of chimeras
        for i in range(num_chimeras):
            # Randomly sample two sequences from the input records
            seq1, seq2 = random.sample(records, 2)
            # Generate a random breakpoint within the range of sequence length
            breakpoint = random.randint(1, len(seq1) - 1)
            # Create the chimeric sequence by joining seq1 and seq2 at the breakpoint
            chimera_seq = seq1[:breakpoint] + seq2[breakpoint:]
            # Create the chimera ID with a prefix and index
            chimera_id = f"{chimera_id_prefix}_{i}"
            # Assign the chimera ID to the sequence
            chimera_seq.id = chimera_id
            # Clear the sequence description
            chimera_seq.description = ""
            # Append the chimeric sequence to the list of chimeras
            chimeras.append(chimera_seq)

            # Write the chimera info to the chimera_info_file
            chimera_info_handle.write(f"{chimera_id}\t{seq1.id}\t{seq2.id}\t{breakpoint}\n")

    # Open the output FASTQ file for writing
    with open(output_file, "w") as output_handle:
        # Write the original records and chimeras to the output file
        SeqIO.write(records + chimeras, output_handle, "fastq")

# Check if the script is being executed as the main program
if __name__ == "__main__":
    # Check if the number of command-line arguments is correct
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} input_fastq output_fastq chimera_info_file num_chimeras")
        sys.exit(1)

    # Parse the command-line arguments
    input_fastq = sys.argv[1]
    output_fastq = sys.argv[2]
    chimera_info_file = sys.argv[3]
    num_chimeras = int(sys.argv[4])

    # Call the generate_chimeras function with the parsed arguments
    generate_chimeras(input_fastq, output_fastq, chimera_info_file, num_chimeras)
