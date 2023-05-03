import random
import sys
from Bio import SeqIO

def generate_chimeras(fasta_file1, fasta_file2, output_file, chimera_info_file, num_chimeras, chimera_id_prefix="chimera"):
    records1 = list(SeqIO.parse(fasta_file1, "fasta"))
    records2 = list(SeqIO.parse(fasta_file2, "fasta"))
    records_combined = records1 + records2
    
    chimeras = []

    with open(chimera_info_file, "w") as chimera_info_handle:
        chimera_info_handle.write("chimera_id\tseq1_id\tseq2_id\tbreakpoint\n")

        for i in range(num_chimeras):
            if i < num_chimeras * 0.1:
                seq1, seq2 = random.sample(records1, 2)
            else:
                seq1, seq2 = random.sample(records_combined, 2)
            breakpoint = random.randint(1, len(seq1) - 1)
            chimera_seq = seq1[:breakpoint] + seq2[breakpoint:]
            chimera_id = f"{chimera_id_prefix}_{i}"
            chimera_seq.id = chimera_id
            chimera_seq.description = ""
            chimeras.append(chimera_seq)

            chimera_info_handle.write(f"{chimera_id}\t{seq1.id}\t{seq2.id}\t{breakpoint}\n")

    with open(output_file, "w") as output_handle:
        SeqIO.write(chimeras, output_handle, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print(f"Usage: {sys.argv[0]} input_fasta1 input_fasta2 output_fasta chimera_info_file num_chimeras")
        sys.exit(1)

    input_fasta1 = sys.argv[1]
    input_fasta2 = sys.argv[2]
    output_fasta = sys.argv[3]
    chimera_info_file = sys.argv[4]
    num_chimeras = int(sys.argv[5])

    generate_chimeras(input_fasta1, input_fasta2, output_fasta, chimera_info_file, num_chimeras)


#run the command with the command below: 
#python generate_chimeras.py input_fasta1.fasta input_fasta2.fasta output_chimeras.fasta chimera_info.txt num_chimeras

