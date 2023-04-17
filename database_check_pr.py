from Bio import SeqIO

def is_nucleotide_sequence(seq):
    nucleotide_symbols = "ACGTURYKMSWBDHVN"
    return all(char.upper() in nucleotide_symbols for char in seq)

input_file = "/home/ali/Music/BOLD_Public.14-Apr-2023.fasta"
output_file = "/home/ali/Music/databse.fasta"

with open(output_file, "w") as output_handle:
    for record in SeqIO.parse(input_file, "fasta"):
        if is_nucleotide_sequence(record.seq):
            SeqIO.write(record, output_handle, "fasta")
