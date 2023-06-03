import os  # module for interacting with the operating system
import random  # module for generating random numbers
from Bio import SeqIO  # module for working with biological sequence data
from Bio.Seq import Seq  # class for representing biological sequences
from Bio.SeqRecord import SeqRecord  # class for representing a sequence record


def generate_chimeras(chimera_id_prefix="chimera"):
    input_directory = os.getcwd()  # get the current working directory
    output_directory = os.path.join(input_directory, "chimeric_reads")  # create the output directory path
    chimera_info_file = os.path.join(output_directory, "chimera_info.tsv")  # create the chimera information file path

    if not os.path.exists(output_directory):  # create the output directory if it doesn't exist
        os.makedirs(output_directory)

    input_files = [os.path.join(input_directory, file) for file in os.listdir(input_directory) if
                   file.endswith(".fasta")]  # get a list of input files in the input directory that end with ".fasta"

    for selected_input_file in input_files:
        records = []
        for input_file in input_files:
            records.extend([(rec, input_file) for rec in SeqIO.parse(input_file, "fasta")])

        main_records = [(rec, selected_input_file) for rec in SeqIO.parse(selected_input_file, "fasta")]
        mixed_records = records

        total_reads = len(main_records)
        num_chimeras = int(total_reads * random.uniform(0.01, 0.03))  # calculate the number of chimeric reads (1-3% of total reads)

        chimeras = []  # create an empty list to store the generated chimeras

        original_ratios = calculate_abundance_ratio(main_records)  # calculate the original abundance ratios of the main records

        with open(chimera_info_file, "a") as chimera_info_handle:  # open the chimera info file in append mode
            if os.path.getsize(chimera_info_file) == 0:  # write the header to the chimera info file if it's empty
                chimera_info_handle.write(
                    "chimera_id\tseq1_id\tseq1_file\tseq2_id\tseq2_file\tbreakpoints\treversed\tratio\tlength\n")

            i = 0
            while i < num_chimeras:
                if i < int(num_chimeras * 0.1):
                    seq1, seq2 = random.sample(main_records, 2)  # select two random sequences from the main records
                else:
                    seq1, seq2 = random.sample(mixed_records, 2)  # select two random sequences from the mixed records

                seq1_rec, seq1_file = seq1
                seq2_rec, seq2_file = seq2

                if seq1_rec.id not in original_ratios:  # skip the iteration if seq1 is not present in the original ratios
                    continue

                breakpoint = random.randint(1, len(seq1_rec) - 1)  # generate a random breakpoint position between seq1 and seq2
                chimera_seq = seq1_rec.seq[:breakpoint] + seq2_rec.seq[breakpoint:]  # create the chimera sequence by combining the sequences at the breakpoint

                if i < int(num_chimeras * 0.04):
                    chimera_seq = chimera_seq.reverse_complement()  # reverse complement the chimera sequence for a certain percentage of chimeras

                if len(seq1_rec) < 1.5 * len(chimera_seq) or len(seq1_rec) > 10 * len(chimera_seq):  # check the length of the chimera sequence to ensure it falls within a valid range
                    continue

                chimera_id = f"{chimera_id_prefix}_{i}"  # generate the chimera ID using the specified prefix and a unique number
                chimera_record = SeqRecord(Seq(str(chimera_seq)), id=chimera_id, description="")  # create a SeqRecord object for the chimera sequence
                chimeras.append((chimera_record, seq1_file))  # add the chimera record and the origin file of seq1 to the list of chimeras

                reversed_status = "yes" if i < int(num_chimeras * 0.04) else "no"  # determine if the chimera is reversed or not based on the percentage

                min_ratio = 0.1  # minimum value for the chimera ratio
                max_ratio = original_ratios[seq1_rec.id] / 1.5  # maximum value for the chimera ratio
                ratio = random.uniform(min_ratio, min(max_ratio, 10))  # select a ratio within the allowed range

                chimera_info_handle.write(
                    f"{chimera_id}\t{seq1_rec.id}\t{os.path.basename(seq1_file)}\t{seq2_rec.id}\t{os.path.basename(seq2_file)}\t{breakpoint}\t{reversed_status}\t{ratio}\t{len(chimera_seq)}\n")  # write the chimera information to the chimera info file

                i += 1

        output_file = os.path.join(output_directory, os.path.basename(selected_input_file))  # create the output file path
        for chimera, file in chimeras:
            position = random.randint(0, len(main_records) - 1)
            main_records.insert(position, (chimera, file))  # insert the generated chimeras randomly into the main records

        with open(output_file, "w") as output_handle:
            SeqIO.write([rec for rec, file in main_records], output_handle, "fasta")  # write the modified main records (including chimeras) to the output file


def calculate_abundance_ratio(records):
    abundance_ratios = {}  # dictionary to store the abundance ratios
    total_reads = len(records)  # calculate the total number of reads

    for record, file in records:
        abundance_ratios[record.id] = abundance_ratios.get(record.id, 0) + 1  # count the occurrence of each sequence ID in the records

    for seq_id, count in abundance_ratios.items():
        abundance_ratios[seq_id] = count / total_reads  # calculate the abundance ratio for each sequence ID

    return abundance_ratios


if __name__ == "__main__":
    generate_chimeras()  # call the generate_chimeras function
