import os
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def generate_chimeras(chimera_id_prefix="chimera"):
    input_directory = os.getcwd()
    output_directory = os.path.join(input_directory, "chimeric_reads")
    chimera_info_file = os.path.join(output_directory, "chimera_info.tsv")

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    input_files = [file for file in os.listdir(input_directory) if file.endswith(".fasta")]

    for selected_input_file in input_files:
        records = []
        for input_file in input_files:
            records.extend([(rec, input_file) for rec in SeqIO.parse(os.path.join(input_directory, input_file), "fasta")])

        main_records = [(rec, selected_input_file) for rec in SeqIO.parse(os.path.join(input_directory, selected_input_file), "fasta")]
        mixed_records = records

        total_reads = len(main_records)
        num_chimeras = int(total_reads * random.uniform(0.01, 0.03))

        chimeras = []
        original_ratios = calculate_abundance_ratio(main_records)

        with open(chimera_info_file, "a") as chimera_info_handle:
            if os.path.getsize(chimera_info_file) == 0:
                chimera_info_handle.write("chimera_id\tseq1_id\tseq1_file\tseq2_id\tseq2_file\tbreakpoints\treversed\tratio\tlength\n")

            i = 0
            last_reverse_status = False
            while i < num_chimeras:
                if i < int(num_chimeras * 0.1):
                    seq1, seq2 = random.sample(main_records, 2)
                else:
                    seq1, seq2 = random.sample(mixed_records, 2)

                seq1_rec, seq1_file = seq1
                seq2_rec, seq2_file = seq2

                if seq1_rec.id not in original_ratios:
                    continue

                breakpoint = random.randint(1, len(seq1_rec) - 1)
                chimera_seq = seq1_rec.seq[:breakpoint] + seq2_rec.seq[breakpoint:]

                if (i % 25 == 0) and (not last_reverse_status):
                    chimera_seq = chimera_seq.reverse_complement()
                    last_reverse_status = True
                else:
                    last_reverse_status = False

                if len(seq1_rec) < 1.5 * len(chimera_seq) or len(seq1_rec) > 10 * len(chimera_seq):
                    continue

                reversed_status = "yes" if last_reverse_status else "no"
                chimera_id = f"{chimera_id_prefix}_{seq1_rec.id}_and_{seq2_rec.id}_at_{breakpoint}_reversed_{reversed_status}_{i}"
                chimera_record = SeqRecord(Seq(str(chimera_seq)), id=chimera_id, description="")
                chimeras.append((chimera_record, seq1_file))

                min_ratio = 0.1
                max_ratio = original_ratios[seq1_rec.id] / 1.5
                ratio = random.uniform(min_ratio, min(max_ratio, 10))

                chimera_info_handle.write(f"{chimera_id}\t{seq1_rec.id}\t{seq1_file}\t{seq2_rec.id}\t{seq2_file}\t{breakpoint}\t{reversed_status}\t{ratio}\t{len(chimera_seq)}\n")

                i += 1

        output_file = os.path.join(output_directory, os.path.basename(selected_input_file))
        for chimera, file in chimeras:
            position = random.randint(0, len(main_records) - 1)
            main_records.insert(position, (chimera, file))

        with open(output_file, "w") as output_handle:
            SeqIO.write([rec for rec, file in main_records], output_handle, "fasta")

def calculate_abundance_ratio(records):
    abundance_ratios = {}
    total_reads = len(records)

    for record, file in records:
        abundance_ratios[record.id] = abundance_ratios.get(record.id, 0) + 1

    for seq_id, count in abundance_ratios.items():
        abundance_ratios[seq_id] = count / total_reads

    return abundance_ratios

if __name__ == "__main__":
    generate_chimeras()
