import random
import sys
from Bio import SeqIO

def generate_chimeras(fastq_file, output_file, chimera_info_file, num_chimeras, chimera_id_prefix="chimera"):
    records = list(SeqIO.parse(fastq_file, "fastq"))
    chimeras = []

    with open(chimera_info_file, "w") as chimera_info_handle:
        chimera_info_handle.write("chimera_id\tseq1_id\tseq2_id\tbreakpoint\n")

        for i in range(num_chimeras):
            seq1, seq2 = random.sample(records, 2)
            breakpoint = random.randint(1, len(seq1) - 1)
            chimera_seq = seq1[:breakpoint] + seq2[breakpoint:]
            chimera_id = f"{chimera_id_prefix}_{i}"
            chimera_seq.id = chimera_id
            chimera_seq.description = ""
            chimeras.append(chimera_seq)

            chimera_info_handle.write(f"{chimera_id}\t{seq1.id}\t{seq2.id}\t{breakpoint}\n")

    with open(output_file, "w") as output_handle:
        SeqIO.write(records + chimeras, output_handle, "fastq")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} input_fastq output_fastq chimera_info_file num_chimeras")
        sys.exit(1)

    input_fastq = sys.argv[1]
    output_fastq = sys.argv[2]
    chimera_info_file = sys.argv[3]
    num_chimeras = int(sys.argv[4])

    generate_chimeras(input_fastq, output_fastq, chimera_info_file, num_chimeras)
