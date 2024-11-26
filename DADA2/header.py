import os
from Bio import SeqIO

# Input folder paths
original_folder = "./original"  # Replace with your folder path containing original FASTA files
header_folder = "."  # Replace with your folder path containing files with new headers
output_folder = "./output"  # Replace with your folder path to save updated FASTA files

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Function to replace headers and save files
def replace_headers(original_folder, header_folder, output_folder):
    for original_file in os.listdir(original_folder):
        # Check for matching file in the header folder
        if original_file in os.listdir(header_folder):
            original_path = os.path.join(original_folder, original_file)
            header_path = os.path.join(header_folder, original_file)
            output_path = os.path.join(output_folder, original_file)

            # Read sequences from original and header files
            original_records = {record.id: record for record in SeqIO.parse(original_path, "fasta")}
            header_records = {record.id: record for record in SeqIO.parse(header_path, "fasta")}

            # Update headers if matches exist
            updated_records = []
            for seq_id, original_record in original_records.items():
                if seq_id in header_records:
                    updated_record = header_records[seq_id]
                    updated_record.seq = original_record.seq  # Replace sequence
                    updated_records.append(updated_record)
                else:
                    # Retain original if no matching header is found
                    updated_records.append(original_record)

            # Write updated sequences to output file with continuous sequence lines
            with open(output_path, "w") as output_handle:
                SeqIO.write(updated_records, output_handle, "fasta-2line")
            print(f"Processed and saved: {original_file}")

# Run the function
replace_headers(original_folder, header_folder, output_folder)
