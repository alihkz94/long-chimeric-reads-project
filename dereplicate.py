#### To check the duplicates and print them out : 
lines = []

line_numbers = {}


with open("/home/ali/Videos/dro_unique.fasta", "r") as f:

	for i, line in enumerate(f):

		if line in lines:

			if line not in line_numbers:

				line_numbers[line] = [i+1]

			else:

				line_numbers[line].append(i+1)

		else:

			lines.append(line)



duplicate_lines = {k: v for k, v in line_numbers.items() if len(v) > 1}


if duplicate_lines:

	print("Duplicate lines found:")

	for line, line_nums in duplicate_lines.items():

		print(f"Line {line_nums}: {line}")

else:

	print("No duplicate lines found.")



### to check and remove the duplicates in the FASTA file: 
lines = []
sequences = {}

with open("/home/ali/Videos/dro.fasta", "r") as f:
    current_sequence = ""
    current_header = ""
    for line in f:
        if line.startswith(">"): # Header line
            if current_sequence:
                sequences[current_header] = current_sequence
                current_sequence = ""
            current_header = line.strip()
        else: # Sequence line
            current_sequence += line.strip()
    # Add the last sequence
    sequences[current_header] = current_sequence

# Remove duplicates
unique_sequences = {}
for header, sequence in sequences.items():
    if sequence not in unique_sequences.values():
        unique_sequences[header] = sequence

# Output the results to a new file
with open("/home/ali/Videos/dro_unique.fasta", "w") as f:
    for header, sequence in unique_sequences.items():
        f.write(header + "\n")
        f.write(sequence + "\n")
