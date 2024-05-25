import sys
import re
from collections import Counter

def find_common_words(species_names, min_count=2):
    """
    Function to find common words in species names.

    Parameters:
    species_names (list): List of species names.
    min_count (int): Minimum count to consider a word as common. Default is 2.

    Returns:
    list: List of common words.
    """
    word_counter = Counter()
    for name in species_names:
        words = re.findall(r'\b\w+\b', name)
        word_counter.update(words)

    common_words = [word for word, count in word_counter.items() if count >= min_count]
    return common_words

def extract_species_name(header):
    """
    Function to extract species name from header.

    Parameters:
    header (str): Header string.

    Returns:
    str: Species name.
    """
    return header.split('|', 1)[-1].strip()

def extract_headers_with_common_words(fasta_file, output_tsv):
    """
    Function to extract headers with common words from a fasta file and write them to a tsv file.

    Parameters:
    fasta_file (str): Path to the input fasta file.
    output_tsv (str): Path to the output tsv file.
    """
    headers = []
    species_names = []

    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip()
                headers.append(header)
                species_names.append(extract_species_name(header))

    common_words = find_common_words(species_names)

    with open(output_tsv, 'w') as out:
        for header, species_name in zip(headers, species_names):
            if any(common_word in species_name for common_word in common_words):
                out.write(f"{header}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_fasta_headers.py <input_fasta> <output_tsv>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_tsv = sys.argv[2]

    extract_headers_with_common_words(fasta_file, output_tsv)
