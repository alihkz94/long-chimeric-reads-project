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
import os
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import numpy as np

class ChimeraGenerator:
    def __init__(self, chimera_id_prefix="chimera"):
        self.chimera_id_prefix = chimera_id_prefix
        # Dataset-specific parameters
        self.min_seq_len = 260
        self.max_seq_len = 1393
        self.avg_seq_len = 547.4
        self.expected_seqs = 28254
        
        # ITS region boundaries adjusted for average sequence length
        self.its_regions = {
            'ITS1': (0, 0.35),    # Slightly expanded ITS1 region
            '5.8S': (0.35, 0.48), # Adjusted 5.8S region
            'ITS2': (0.48, 1.0)   # Adjusted ITS2 region
        }
        
        # PacBio-specific parameters
        self.chimera_rate_bounds = (0.08, 0.12)  # 8-12% chimera rate for PacBio
        self.length_ratio_bounds = (0.6, 1.8)    # More stringent length ratio constraints
        
    def validate_chimera_length(self, chimera_seq):
        """
        Validate chimera length against dataset parameters
        """
        length = len(chimera_seq)
        # Allow slightly shorter/longer sequences than original bounds
        min_allowed = max(self.min_seq_len * 0.9, 50)  # Never below 50
        max_allowed = min(self.max_seq_len * 1.1, 1500) # Never above 1500
        
        return min_allowed <= length <= max_allowed

    def get_weighted_breakpoint(self, seq_length):
        """
        Generate breakpoints with dataset-specific weighting
        """
        weights = []
        its1_peak = int(seq_length * 0.25)  # Peak probability in ITS1
        its2_peak = int(seq_length * 0.65)  # Peak probability in ITS2
        
        for i in range(seq_length):
            pos_ratio = i / seq_length
            
            # Create gaussian-like distributions around typical breakpoints
            its1_weight = np.exp(-0.5 * ((i - its1_peak) / (seq_length * 0.1)) ** 2)
            its2_weight = np.exp(-0.5 * ((i - its2_peak) / (seq_length * 0.1)) ** 2)
            
            if self.its_regions['ITS1'][0] <= pos_ratio <= self.its_regions['ITS1'][1]:
                weights.append(its1_weight)
            elif self.its_regions['ITS2'][0] <= pos_ratio <= self.its_regions['ITS2'][1]:
                weights.append(its2_weight)
            elif self.its_regions['5.8S'][0] <= pos_ratio <= self.its_regions['5.8S'][1]:
                weights.append(0.2)  # Low probability in 5.8S
            else:
                weights.append(0.3)  # Moderate probability in other regions
                
        weights = np.array(weights) / sum(weights)
        return np.random.choice(range(seq_length), p=weights)

    def generate_chimeras(self):
        input_directory = os.getcwd()
        output_directory = os.path.join(input_directory, "chimeric_reads")
        chimera_info_file = os.path.join(output_directory, "chimera_info.tsv")
        
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        input_files = [f for f in os.listdir(input_directory) if f.endswith(".fasta")]
        chimera_stats = defaultdict(int)
        
        for selected_input_file in input_files:
            # Load and validate input sequences
            records = []
            for input_file in input_files:
                file_records = list(SeqIO.parse(os.path.join(input_directory, input_file), "fasta"))
                # Filter sequences based on length constraints
                valid_records = [(rec, input_file) for rec in file_records 
                               if self.min_seq_len <= len(rec.seq) <= self.max_seq_len]
                records.extend(valid_records)

            main_records = [(rec, selected_input_file) for rec in 
                           SeqIO.parse(os.path.join(input_directory, selected_input_file), "fasta")
                           if self.min_seq_len <= len(rec.seq) <= self.max_seq_len]
            
            mixed_records = records
            total_reads = len(main_records)
            
            # Calculate chimera numbers based on PacBio-specific rates
            chimera_rate = random.uniform(*self.chimera_rate_bounds)
            num_chimeras = int(total_reads * chimera_rate)
            
            chimeras = []
            original_ratios = self._calculate_abundance_ratio(main_records)

            with open(chimera_info_file, "a") as chimera_info_handle:
                if os.path.getsize(chimera_info_file) == 0:
                    chimera_info_handle.write("chimera_id\tseq1_id\tseq1_file\tseq2_id\tseq2_file\t"
                                           "breakpoint\tbreakpoint_region\treversed\tratio\tlength\t"
                                           "parent1_length\tparent2_length\tits_region_break\n")

                i = 0
                attempts = 0
                max_attempts = num_chimeras * 3  # Limit attempts to avoid infinite loops
                
                while i < num_chimeras and attempts < max_attempts:
                    attempts += 1
                    
                    # Select parent sequences with length-based weighting
                    if i < int(num_chimeras * 0.15):
                        seq1, seq2 = random.sample(main_records, 2)
                    else:
                        seq1, seq2 = random.sample(mixed_records, 2)

                    seq1_rec, seq1_file = seq1
                    seq2_rec, seq2_file = seq2

                    if seq1_rec.id not in original_ratios:
                        continue

                    # Generate biased breakpoint
                    breakpoint = self.get_weighted_breakpoint(len(seq1_rec))
                    pos_ratio = breakpoint / len(seq1_rec)
                    
                    # Determine breakpoint region
                    breakpoint_region = next(region for region, (start, end) in 
                                          self.its_regions.items() if start <= pos_ratio <= end)

                    # Create chimeric sequence
                    chimera_seq = seq1_rec.seq[:breakpoint] + seq2_rec.seq[breakpoint:]

                    # Validate chimera
                    if not self.validate_chimera_length(chimera_seq):
                        chimera_stats['length_rejected'] += 1
                        continue

                    # Length ratio validation
                    parent_ratio = len(seq1_rec) / len(seq2_rec)
                    if not (self.length_ratio_bounds[0] <= parent_ratio <= self.length_ratio_bounds[1]):
                        chimera_stats['ratio_rejected'] += 1
                        continue

                    # PacBio-specific reverse complement probability
                    should_reverse = random.random() < 0.08  # 8% chance based on PacBio characteristics
                    if should_reverse:
                        chimera_seq = chimera_seq.reverse_complement()

                    # Generate chimera record
                    reversed_status = "yes" if should_reverse else "no"
                    chimera_id = (f"{self.chimera_id_prefix}_{seq1_rec.id}_and_{seq2_rec.id}_"
                                f"at_{breakpoint}_{breakpoint_region}_reversed_{reversed_status}_{i}")
                    
                    chimera_record = SeqRecord(Seq(str(chimera_seq)), id=chimera_id, description="")
                    chimeras.append((chimera_record, seq1_file))

                    # Calculate abundance with realistic distribution
                    parent_abundance = original_ratios[seq1_rec.id]
                    # Log-normal distribution for chimera abundance
                    ratio = np.random.lognormal(mean=-1.5, sigma=0.5) * parent_abundance
                    ratio = min(ratio, 0.5 * parent_abundance)  # Cap at 50% of parent abundance

                    # Log chimera information
                    chimera_info_handle.write(
                        f"{chimera_id}\t{seq1_rec.id}\t{seq1_file}\t{seq2_rec.id}\t{seq2_file}\t"
                        f"{breakpoint}\t{breakpoint_region}\t{reversed_status}\t{ratio:.4f}\t"
                        f"{len(chimera_seq)}\t{len(seq1_rec)}\t{len(seq2_rec)}\t{pos_ratio:.3f}\n"
                    )

                    i += 1
                    chimera_stats['accepted'] += 1

            # Write output with chimeras
            output_file = os.path.join(output_directory, os.path.basename(selected_input_file))
            self._write_output_with_chimeras(output_file, main_records, chimeras)
            
        # Write detailed statistics
        stats_file = os.path.join(output_directory, "chimera_generation_stats.txt")
        self._write_statistics(stats_file, chimera_stats)

    def _calculate_abundance_ratio(self, records):
        """Calculate abundance ratios with log-normal distribution"""
        abundance_ratios = defaultdict(int)
        total_reads = len(records)

        for record, _ in records:
            abundance_ratios[record.id] += 1

        return {seq_id: count/total_reads for seq_id, count in abundance_ratios.items()}

    def _write_output_with_chimeras(self, output_file, main_records, chimeras):
        """Write output with improved chimera distribution"""
        # Insert chimeras following a more natural distribution
        positions = np.random.choice(
            len(main_records), 
            size=len(chimeras), 
            replace=True
        )
        
        for (chimera, _), pos in zip(chimeras, sorted(positions)):
            main_records.insert(pos, (chimera, os.path.basename(output_file)))

        with open(output_file, "w") as output_handle:
            SeqIO.write([rec for rec, _ in main_records], output_handle, "fasta-2line")

    def _write_statistics(self, stats_file, stats):
        """Write comprehensive statistics about chimera generation"""
        with open(stats_file, "w") as f:
            f.write("Chimera Generation Statistics\n")
            f.write("===========================\n")
            f.write(f"Input Dataset Properties:\n")
            f.write(f"  Total sequences: {self.expected_seqs:,}\n")
            f.write(f"  Average length: {self.avg_seq_len:.1f}\n")
            f.write(f"  Length range: {self.min_seq_len}-{self.max_seq_len}\n\n")
            
            f.write("Chimera Generation Results:\n")
            f.write(f"  Total chimeras accepted: {stats['accepted']}\n")
            f.write(f"  Rejected due to length: {stats['length_rejected']}\n")
            f.write(f"  Rejected due to invalid ratio: {stats['ratio_rejected']}\n")
            
            total_attempted = sum(stats.values())
            if total_attempted > 0:
                success_rate = stats['accepted']/total_attempted
                f.write(f"  Overall success rate: {success_rate:.2%}\n")

if __name__ == "__main__":
    generator = ChimeraGenerator()
    generator.generate_chimeras()
