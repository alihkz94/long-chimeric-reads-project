import os
import random
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def generate_chimeras(input_files, output_file, chimera_info_file, num_chimeras, chimera_id_prefix="chimera"):
    records = []
    for input_file in input_files:
        records.extend(list(SeqIO.parse(input_file, "fasta")))

    main_records = records[:len(records)//10]  # 10% of records from the first file
    mixed_records = records  # 90% of chimeric reads from a mix of sequences

    chimeras = []

    with open(chimera_info_file, "w") as chimera_info_handle:
        chimera_info_handle.write("chimera_id\tseq1_id\tseq2_id\tbreakpoint\treversed\n")

        for i in range(num_chimeras):
            seq1, seq2 = random.sample(mixed_records, 2)
            breakpoint = random.randint(1, len(seq1) - 1)
            chimera_seq = seq1[:breakpoint] + seq2[breakpoint:]

            if i < int(num_chimeras * 0.04):  # 4% of 10% portion with reversed sequences
                chimera_seq = chimera_seq.reverse_complement()

            chimera_id = f"{chimera_id_prefix}_{i}"
            chimera_seq.id = chimera_id
            chimera_seq.description = ""
            chimeras.append(chimera_seq)

            reversed_status = "yes" if i < int(num_chimeras * 0.04) else "no"
            chimera_info_handle.write(f"{chimera_id}\t{seq1.id}\t{seq2.id}\t{breakpoint}\t{reversed_status}\n")

    with open(output_file, "w") as output_handle:
        SeqIO.write(main_records + chimeras, output_handle, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} input_directory output_fasta chimera_info_file num_chimeras")
        sys.exit(1)

    input_directory = sys.argv[1]
    output_fasta = sys.argv[2]
    chimera_info_file = sys.argv[3]
    num_chimeras = int(sys.argv[4])

    input_files = [os.path.join(input_directory, file) for file in os.listdir(input_directory) if file.endswith(".fasta")]

    generate_chimeras(input_files, output_fasta, chimera_info_file, num_chimeras)
