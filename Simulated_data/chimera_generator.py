"""
=====================================================================
Chimera Generator: Creating Chimeric Reads from Multiple FASTA Files
=====================================================================

**Description:**
This script generates chimeric reads by combining sequences from multiple FASTA files.
It randomly selects pairs of sequences, determines breakpoints, and creates new sequences by
joining parts of each sequence. The generated chimeras are saved in a designated output
directory, and detailed information about each chimera is logged for further analysis.

**Features:**
- Combines sequences from multiple FASTA files to create chimeric reads.
- Configurable prefix for chimera IDs.
- Ensures chimeric sequences meet specified length and ratio criteria.
- Optionally reverse complements some chimeras for diversity.
- Logs detailed information about each generated chimera, including parent sequences and breakpoints.

**Usage:**
1. **Prepare Input:**
   - Place all input FASTA files in the working directory where the script will be executed.

2. **Run the Script:**
   ```bash
   python generate_chimeras.py
   ```
"""   
#load libraries
import os
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def generate_chimeras(chimera_id_prefix="chimera"):

    # Set up input and output directories
    input_directory = os.getcwd()
    output_directory = os.path.join(input_directory, "chimeric_reads")
    chimera_info_file = os.path.join(output_directory, "chimera_info.tsv")

    # Create output directory if needed
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Get list of input FASTA files
    input_files = [file for file in os.listdir(input_directory) if file.endswith(".fasta")]

    # Track short chimeras to avoid too many
    short_chimeras_count = 0

    for selected_input_file in input_files:

        # Parse all input files
        records = []
        for input_file in input_files:
            records.extend([(rec, input_file) for rec in SeqIO.parse(os.path.join(input_directory, input_file), "fasta")])

        # Parse seqs from selected "main" file
        main_records = [(rec, selected_input_file) for rec in SeqIO.parse(os.path.join(input_directory, selected_input_file), "fasta")]
        mixed_records = records

        # Determine number of chimeras to generate
        total_reads = len(main_records)
        num_chimeras = int(total_reads * random.uniform(0.01, 0.03))

        chimeras = []
        original_ratios = calculate_abundance_ratio(main_records)

        with open(chimera_info_file, "a") as chimera_info_handle:
        
            # Write header if new file
            if os.path.getsize(chimera_info_file) == 0:
                chimera_info_handle.write("chimera_id\tseq1_id\tseq1_file\tseq2_id\tseq2_file\tbreakpoints\treversed\tratio\tlength\n")

            # Generate chimeras
            i = 0
            last_reverse_status = False
            while i < num_chimeras:

                # Select parent seqs
                if i < int(num_chimeras * 0.1):
                    seq1, seq2 = random.sample(main_records, 2)  
                else:
                    seq1, seq2 = random.sample(mixed_records, 2)

                # Get seq records
                seq1_rec, seq1_file = seq1
                seq2_rec, seq2_file = seq2

                # Skip if no original ratio data
                if seq1_rec.id not in original_ratios:
                    continue

                # Generate breakpoint
                breakpoint = random.randint(1, len(seq1_rec) - 1)
                chimera_seq = seq1_rec.seq[:breakpoint] + seq2_rec.seq[breakpoint:]

                # Check if chimera components are too short
                if len(seq1_rec.seq[:breakpoint]) < 50 or len(seq2_rec.seq[breakpoint:]) < 50:
                    short_chimeras_count += 1
                    if short_chimeras_count / num_chimeras > 0.1:
                        continue
                else:
                    short_chimeras_count = max(0, short_chimeras_count - 1) 

                # Reverse complement some chimeras
                if (i % 25 == 0) and (not last_reverse_status):
                    chimera_seq = chimera_seq.reverse_complement()
                    last_reverse_status = True
                else:
                    last_reverse_status = False

                # Skip if chimera size is too skewed
                if len(seq1_rec) < 1.5 * len(chimera_seq) or len(seq1_rec) > 10 * len(chimera_seq):
                    continue

                # Generate chimera ID and record
                reversed_status = "yes" if last_reverse_status else "no"
                chimera_id = f"{chimera_id_prefix}_{seq1_rec.id}_and_{seq2_rec.id}_at_{breakpoint}_reversed_{reversed_status}_{i}"
                chimera_record = SeqRecord(Seq(str(chimera_seq)), id=chimera_id, description="")
                chimeras.append((chimera_record, seq1_file))

                # Calculate abundance
                min_ratio = 0.1
                max_ratio = original_ratios[seq1_rec.id] / 1.5
                ratio = random.uniform(min_ratio, min(max_ratio, 10))

                # Log chimera info
                chimera_info_handle.write(f"{chimera_id}\t{seq1_rec.id}\t{seq1_file}\t{seq2_rec.id}\t{seq2_file}\t{breakpoint}\t{reversed_status}\t{ratio}\t{len(chimera_seq)}\n")

                i += 1

        # Write chimeras to output FASTA
        output_file = os.path.join(output_directory, os.path.basename(selected_input_file))
        for chimera, _ in chimeras:
            position = random.randint(0, len(main_records) - 1)
            main_records.insert(position, (chimera, selected_input_file))

        with open(output_file, "w") as output_handle:
            SeqIO.write([rec for rec, file in main_records], output_handle, "fasta")
            
# Calculate original abundance ratios
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
