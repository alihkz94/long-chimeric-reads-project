import pandas as pd
import argparse
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob
import logging

# Basic logging configuration
logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    """
    Parses command line arguments.
    Requires three directories: chimeras_dir, blast_output_dir, and nonchimeric_dir.
    """
    parser = argparse.ArgumentParser(description='Chimeric Read Analysis Module')
    parser.add_argument('--chimeras_dir', help='Directory containing chimeras files', required=True)
    parser.add_argument('--blast_output_dir', help='Directory containing BLAST output files', required=True)
    parser.add_argument('--nonchimeric_dir', help='Directory containing nonchimeric files', required=True)
    return parser.parse_args()

def parse_blast_output(file_path):
    """
    Reads and processes a BLAST output file.
    Extracts relevant information into a DataFrame for further processing.
    """
    columns = ['qseqid', 'stitle', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'length', 'nident', 'mismatch', 'gapopen', 'gaps', 'sstrand', 'qcovs', 'pident']
    try:
        return pd.read_csv(file_path, sep='\+', names=columns, header=None, engine='python')
    except Exception as e:
        logging.error(f"Error reading file {file_path}: {e}")
        return pd.DataFrame(columns=columns)

def process_files(chimeras_dir, blast_output_dir, nonchimeric_dir):
    """
    Main processing function.
    Iterates through chimeras, compares with BLAST output, and recovers nonchimeric sequences.
    """
    # Directory setup for recovered reads
    recovered_dir = "rescued_reads"
    os.makedirs(recovered_dir, exist_ok=True)

    report_data = []
    total_rescued = 0

    # Processing each chimeras file
    for chimeras_file in glob.glob(f"{chimeras_dir}/*.fasta"):
        # File path management
        base_name = os.path.basename(chimeras_file).split('.')[0]
        blast_file = f"{blast_output_dir}/{base_name}.txt"
        recovered_file = f"{recovered_dir}/{base_name}_rescued.fasta"

        if os.path.exists(blast_file):
            blast_data = parse_blast_output(blast_file)
            blast_data['qcovs'] = pd.to_numeric(blast_data['qcovs'], errors='coerce')

            # Sequence processing and recovery
            sequences = {record.id: str(record.seq) for record in SeqIO.parse(chimeras_file, "fasta")}
            rescued_seqs = []

            # Recovering nonchimeric sequences
            for _, record in blast_data.iterrows():
                if record['qcovs'] > 91:
                    seq_id = record['qseqid']
                    if seq_id in sequences:
                        rescued_seqs.append(SeqRecord(Seq(sequences[seq_id]), id=seq_id, description=""))

            rescued_count = len(rescued_seqs)
            total_rescued += rescued_count

            # Write rescued sequences without line breaks
            if rescued_seqs:
                with open(recovered_file, 'w') as output_handle:
                    for seq in rescued_seqs:
                        output_handle.write(f">{seq.id}\n{seq.seq}\n")

            report_data.append([base_name, rescued_count])

    # Writing a report
    with open('report.txt', 'w') as report_file:
        report_file.write('Filename,Rescued\n')
        for line in report_data:
            report_file.write(','.join(map(str, line)) + '\n')
        report_file.write(f'Total Rescued Sequences: {total_rescued}')

def main():
    """
    Main function to initialize and execute script processes.
    """
    args = parse_args()
    process_files(args.chimeras_dir, args.blast_output_dir, args.nonchimeric_dir)

if __name__ == "__main__":
    main()

#usage: python blast_recovery.py --chimeras_dir /media/ali/data_store/test_blast/chimeras/ --blast_output_dir /media/ali/data_store/test_blast/blast_output/ --nonchimeric_dir /media/ali/data_store/test_blast/