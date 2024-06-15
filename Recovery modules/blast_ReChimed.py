import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob
import logging

# Basic logging configuration
logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_custom_blast_output(file_path):
    """
    Reads and processes a BLAST output file with non-standard delimiters.
    Extracts relevant information into a DataFrame for further processing.
    """
    try:
        data = []
        with open(file_path, 'r') as file:
            for line in file:
                parts = line.strip().split('+')
                if len(parts) == 17:
                    data.append(parts)
        
        columns = ['qseqid', 'stitle', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'length', 'nident', 'mismatch', 'gapopen', 'gaps', 'sstrand', 'qcovs', 'pident']
        return pd.DataFrame(data, columns=columns)
    except Exception as e:
        logging.error(f"Error reading file {file_path}: {e}")
        return pd.DataFrame(columns=columns)

def calculate_adjusted_coverage(row):
    """
    Calculates adjusted coverage based on the given row data.
    """
    slen = int(row['slen'])
    qlen = int(row['qlen'])
    sstart = int(row['sstart'])
    send = int(row['send'])
    qcovs = float(row['qcovs'])

    if (slen - qlen) < 0:
        adjCov = ((send - sstart + 1) / slen) * 100
    else:
        adjCov = qcovs

    return adjCov

def process_files():
    """
    Main processing function.
    Iterates through fasta files, compares with BLAST output, and classifies sequences.
    """
    report_data = []

    # Processing each fasta file
    for fasta_file in glob.glob("*.fasta"):
        # File path management
        base_name = os.path.basename(fasta_file).split('.')[0]
        blast_file = f"{base_name}.txt"

        if os.path.exists(blast_file):
            blast_data = parse_custom_blast_output(blast_file)
            blast_data['qcovs'] = pd.to_numeric(blast_data['qcovs'], errors='coerce')

            # Sequence processing and classification
            sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
            nonchimeric_seqs = []
            chimeric_seqs = []

            # Classifying sequences
            for _, record in blast_data.iterrows():
                adjCov = calculate_adjusted_coverage(record)

                seq_id = record['qseqid']
                if seq_id in sequences:
                    if adjCov > 91:
                        nonchimeric_seqs.append(SeqRecord(Seq(sequences[seq_id]), id=seq_id, description=""))
                    else:
                        chimeric_seqs.append(SeqRecord(Seq(sequences[seq_id]), id=seq_id, description=""))

            nonchimeric_count = len(nonchimeric_seqs)
            chimeric_count = len(chimeric_seqs)

            # Write nonchimeric sequences without line breaks
            if nonchimeric_seqs:
                with open(f"{base_name}_nonchimeric.fasta", 'w') as output_handle:
                    for seq in nonchimeric_seqs:
                        output_handle.write(f">{seq.id}\n{seq.seq}\n")

            # Write chimeric sequences without line breaks
            if chimeric_seqs:
                with open(f"{base_name}_chimeric.fasta", 'w') as output_handle:
                    for seq in chimeric_seqs:
                        output_handle.write(f">{seq.id}\n{seq.seq}\n")

            report_data.append([base_name, nonchimeric_count, chimeric_count])

    # Writing a report
    with open('report.txt', 'w') as report_file:
        report_file.write('Filename,Nonchimeric,Chimeric\n')
        for line in report_data:
            report_file.write(','.join(map(str, line)) + '\n')

def main():
    """
    Main function to initialize and execute script processes.
    """
    process_files()

if __name__ == "__main__":
    main()
