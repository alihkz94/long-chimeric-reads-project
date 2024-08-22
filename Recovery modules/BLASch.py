###### Script to compare the genus information from two BLAST output files and generate reports for all categories ######
import os
import pandas as pd

def parse_blast_output(file_path):
    genus_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('ERR'):
                parts = line.split('+')
                sequence_id = parts[0]
                taxonomy = parts[1].split(';')
                
                # Search specifically for genus level (g__)
                genus = 'g__unclassified'
                for taxon in taxonomy:
                    if taxon.startswith('g__'):
                        genus = taxon.strip().lower()  # Get the genus and normalize to lowercase
                        break
                
                # Add the parsed genus to the dictionary
                genus_dict[sequence_id] = genus

    return genus_dict

def compare_genus(genus_dict1, genus_dict2):
    matching_genus = {}
    different_genus = {}
    unclassified_genus = {}
    unclassified_count = 0

    for seq_id, genus1 in genus_dict1.items():
        if seq_id in genus_dict2:
            genus2 = genus_dict2[seq_id]

            if genus1 == 'g__unclassified' and genus2 == 'g__unclassified':
                unclassified_count += 1
                unclassified_genus[seq_id] = (genus1, genus2)
            elif genus1 == genus2:
                matching_genus[seq_id] = genus1
            else:
                different_genus[seq_id] = (genus1, genus2)

    return matching_genus, different_genus, unclassified_genus, unclassified_count

# Define the directories
begin_dir = '/home/ali/Documents/pipe/long/begin'
end_dir = '/home/ali/Documents/pipe/long/end'

# Initialize a list to store the results
results = []
difference_details = []
similar_details = []
unclassified_details = []

# Loop over each file in the begin directory
for file_name in os.listdir(begin_dir):
    if file_name.endswith('_begin.txt'):
        base_name = file_name.replace('_begin.txt', '')

        # Construct the corresponding file path in the end directory
        begin_file = os.path.join(begin_dir, f"{base_name}_begin.txt")
        end_file = os.path.join(end_dir, f"{base_name}_end.txt")

        print(f"Processing file: {base_name}")

        # Ensure the corresponding file exists in the end directory
        if os.path.exists(end_file):
            print(f"Found matching file in end folder for: {base_name}")
            
            # Parse the BLAST output files
            genus_begin = parse_blast_output(begin_file)
            genus_end = parse_blast_output(end_file)

            # Compare the genus information
            matching_genus, different_genus, unclassified_genus, unclassified_count = compare_genus(genus_begin, genus_end)

            # Store the results
            results.append({
                'Base Name': base_name,
                'Similar': len(matching_genus),
                'Different': len(different_genus),
                'Unclassified (Both)': unclassified_count
            })
            
            # Prepare details for the text reports
            for seq_id, (genus1, genus2) in different_genus.items():
                difference_details.append(f"{seq_id}\t{genus1}\t{genus2}\n")
            
            for seq_id, genus in matching_genus.items():
                similar_details.append(f"{seq_id}\t{genus}\t{genus}\n")
            
            for seq_id, (genus1, genus2) in unclassified_genus.items():
                unclassified_details.append(f"{seq_id}\t{genus1}\t{genus2}\n")

        else:
            print(f"Matching file not found in end folder for: {base_name}")

# Convert the results into a DataFrame
df = pd.DataFrame(results)

# Save the results to an Excel file
output_file = '/home/ali/Documents/pipe/long/genus_comparison_report.xlsx'
df.to_excel(output_file, index=False)

# Generate the text report for differences
difference_report_file = '/home/ali/Documents/pipe/long/genus_differences_report.txt'
with open(difference_report_file, 'w') as report_file:
    report_file.write("qseqid\tbegin\tend\n")
    report_file.writelines(difference_details)

# Generate the text report for similar genera
similar_report_file = '/home/ali/Documents/pipe/long/genus_similar_report.txt'
with open(similar_report_file, 'w') as report_file:
    report_file.write("qseqid\tbegin\tend\n")
    report_file.writelines(similar_details)

# Generate the text report for unclassified genera
unclassified_report_file = '/home/ali/Documents/pipe/long/genus_unclassified_report.txt'
with open(unclassified_report_file, 'w') as report_file:
    report_file.write("qseqid\tbegin\tend\n")
    report_file.writelines(unclassified_details)

if df.empty:
    print("No matching data found. Please check the input files and ensure they contain valid genus information.")
else:
    print(f"Report successfully generated: {output_file}")
    print(f"Differences report generated: {difference_report_file}")
    print(f"Similar report generated: {similar_report_file}")
    print(f"Unclassified report generated: {unclassified_report_file}")
