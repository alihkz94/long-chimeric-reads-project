import os
import random
import sys
from Bio import SeqIO
from Bio.Seq import Seq

# Define the generate_chimeras function with appropriate parameters
def generate_chimeras(input_files, output_file, chimera_info_file, num_chimeras, chimera_id_prefix="chimera"):
    # Initialize an empty list to store all records from input files
    records = []
    # Read all the input FASTA files and store the records in the list
    for input_file in input_files:
        records.extend(list(SeqIO.parse(input_file, "fasta")))

    # Select 10% of records from the first file
    main_records = records[:len(records)//10]
    # Use all records for 90% chimeric reads
    mixed_records = records

    # Initialize an empty list to store chimeras
    chimeras = []

    # Open the chimera_info_file for writing
    with open(chimera_info_file, "w") as chimera_info_handle:
        # Write the header for the chimera info file
        chimera_info_handle.write("chimera_id\tseq1_id\tseq2_id\tbreakpoint\treversed\n")

        # Iterate for the desired number of chimeras
        for i in range(num_chimeras):
            # Randomly sample two sequences from the mixed_records
            seq1, seq2 = random.sample(mixed_records, 2)
            # Generate a random breakpoint within the range of sequence length
            breakpoint = random.randint(1, len(seq1) - 1)
            # Create the chimeric sequence by joining seq1 and seq2 at the breakpoint
            chimera_seq = seq1[:breakpoint] + seq2[breakpoint:]

            # Reverse and complement 4% of the 10% portion of chimeric sequences
            if i < int(num_chimeras * 0.04):
                chimera_seq = chimera_seq.reverse_complement()

            # Create the chimera ID with a prefix and index
            chimera_id = f"{chimera_id_prefix}_{i}"
            # Assign the chimera ID to the sequence
            chimera_seq.id = chimera_id
            # Clear the sequence description
            chimera_seq.description = ""
            # Append the chimeric sequence to the list of chimeras
            chimeras.append(chimera_seq)

            # Write the chimera info to the chimera_info_file
            reversed_status = "yes" if i < int(num_chimeras * 0.04) else "no"
            chimera_info_handle.write(f"{chimera_id}\t{seq1.id}\t{seq2.id}\t{breakpoint}\t{reversed_status}\n")

    # Open the output FASTA file for writing
    with open(output_file, "w") as output_handle:
        # Write the main records and chimeras to the output file
        SeqIO.write(main_records + chimeras, output_handle, "fasta")

# Check if the script is being executed as the main program
if __name__ == "__main__":
    # Check if the number of command-line arguments is correct
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} input_directory output_fasta chimera_info_file num_chimeras")
        sys.exit(1)

    # Parse the command-line arguments
    input_directory = sys.argv[1]
    output_fasta = sys.argv[2]
    chimera_info_file = sys.argv[3]
    num




#run the command with the command below: 
#python generate_chimeras.py input_fasta1.fasta input_fasta2.fasta output_chimeras.fasta chimera_info.txt num_chimeras

