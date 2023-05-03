import random
import sys
from Bio import SeqIO

# Define the generate_chimeras function with appropriate parameters
def generate_chimeras(fasta_file1, fasta_file2, output_file, chimera_info_file, num_chimeras, chimera_id_prefix="chimera"):
    # Read the input FASTA files and store the records in lists
    records1 = list(SeqIO.parse(fasta_file1, "fasta"))
    records2 = list(SeqIO.parse(fasta_file2, "fasta"))
    
    # Combine the records from both files
    records_combined = records1 + records2
    
    # Initialize an empty list to store chimeras
    chimeras = []

    # Open the chimera_info_file for writing
    with open(chimera_info_file, "w") as chimera_info_handle:
        # Write the header for the chimera info file
        chimera_info_handle.write("chimera_id\tseq1_id\tseq2_id\tbreakpoint\n")

        # Iterate for the desired number of chimeras
        for i in range(num_chimeras):
            # Select sequences for chimera generation based on the desired proportions
            if i < num_chimeras * 0.1:
                seq1, seq2 = random.sample(records1, 2)
            else:
                seq1, seq2 = random.sample(records_combined, 2)
            
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

    # Open the output FASTA file for writing
    with open(output_file, "w") as output_handle:
        # Write the chimeras to the output file
        SeqIO.write(chimeras, output_handle, "fasta")

# Check if the script is being executed as the main program
if __name__ == "__main__":
    # Check if the number of command-line arguments is correct
    if len(sys.argv) != 6:
        print(f"Usage: {sys.argv[0]} input_fasta1 input_fasta2 output_fasta chimera_info_file num_chimeras")
        sys.exit(1)

    # Parse the command-line arguments
    input_fasta1 = sys.argv[1]
    input_fasta2 = sys.argv[2]
    output_fasta = sys.argv[3]
    chimera_info_file = sys.argv[4]
    num_chimeras = int(sys.argv[5])

    # Call the generate_chimeras function with the parsed arguments
    generate_chimeras(input_fasta1, input_fasta2, output_fasta, chimera_info_file, num_chimeras)



#run the command with the command below: 
#python generate_chimeras.py input_fasta1.fasta input_fasta2.fasta output_chimeras.fasta chimera_info.txt num_chimeras

