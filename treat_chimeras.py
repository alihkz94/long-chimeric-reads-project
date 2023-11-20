import pandas as pd
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.pyplot as plt

# Function to parse arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Chimeric Read Analysis Module')
    parser.add_argument('--blast_output', help='BLAST output file', required=True)
    parser.add_argument('--input_fasta', help='Input FASTA file', required=True)
    parser.add_argument('--output', help='Output CSV file', required=True)
    parser.add_argument('--min_occurrence', type=int, default=2, help='Minimum sequence occurrence')
    return parser.parse_args()

# Function to parse BLAST output
def parse_blast_output(file_path):
    columns = ['qseqid', 'stitle', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'length', 'nident', 'mismatch', 'gapopen', 'gaps', 'sstrand', 'qcovs', 'pident']
    return pd.read_csv(file_path, sep='\t', names=columns, header=0)

# Function to process each BLAST record
def process_blast_record(record, min_occurrence, sequence_counts, sequences):
    adjCov = record['qcovs']
    if record['slen'] - record['qlen'] < 0:
        adjCov = ((record['send'] - record['sstart'] + 1) / record['slen']) * 100

    occurrence = sequence_counts.get(record['qseqid'], 0)
    is_false_positive = adjCov > 91
    is_rescued = occurrence > min_occurrence

    return {
        'SeqID': record['qseqid'],
        'IsChimeric': is_false_positive or is_rescued,
        'Sequence': sequences.get(record['qseqid'], '')
    }

# Function to generate FASTA file from processed data
def generate_fasta(processed_data, output_path):
    records = [SeqRecord(Seq(row['Sequence']), id=row['SeqID'], description="") for index, row in processed_data.iterrows() if not row['IsChimeric']]
    SeqIO.write(records, output_path, "fasta")

# Function to plot data
def plot_data(processed_data):
    lengths = processed_data['Sequence'].apply(len)
    plt.hist(lengths, bins=30)
    plt.title('Distribution of Sequence Lengths')
    plt.xlabel('Length')
    plt.ylabel('Frequency')
    plt.show()

# Main function
def main():
    args = parse_args()
    blast_data = parse_blast_output(args.blast_output)

    sequences = {record.id: str(record.seq) for record in SeqIO.parse(args.input_fasta, "fasta")}
    sequence_counts = blast_data['qseqid'].value_counts().to_dict()

    processed_data = pd.DataFrame([process_blast_record(record, args.min_occurrence, sequence_counts, sequences) for _, record in blast_data.iterrows()])
    processed_data.to_csv(args.output, index=False)

    generate_fasta(processed_data, "processed_sequences.fasta")
    plot_data(processed_data)

if __name__ == "__main__":
    main()

#The script adjusts the coverage value (adjCov) if the subject sequence is shorter than the query sequence. 
#It also checks if the sequence occurs more frequently than the minimum occurrence threshold. Sequences are flagged as chimeric based on these criteria.
