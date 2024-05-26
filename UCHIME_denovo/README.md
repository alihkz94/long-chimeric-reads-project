UCHIME-related scripts and files
## Chimera_rescue explanation 
The Python code is a function named `rescue_chimeric_sequences` that is designed to process a set of FASTA files, aggregate sequences, filter them based on a minimum occurrence threshold, and write the filtered sequences to an output file. 

The function takes three parameters: `input_dir` (the directory containing the FASTA files), `min_occurrence` (the minimum occurrence threshold for a sequence to be considered), and `output_dir` (the directory where the output file will be written).

The function starts by initializing two empty lists: `aggregated_sequences` and `sequence_occurrences`. The former will hold all the sequences found in the FASTA files, while the latter, a Counter object, will keep track of how many times each sequence occurs.

The function then iterates over all files in the `input_dir` directory. If a file ends with ".fasta" or ".fa", it is considered a FASTA file and is processed. The file is opened and parsed using the `SeqIO.parse` function from the BioPython library, which returns a list of sequence records. Each sequence record is then converted to a string and added to the `aggregated_sequences` list, and its occurrence count is incremented in the `sequence_occurrences` Counter.

After all sequences have been aggregated, the function creates a set of 'rescued' sequences. These are the sequences that occur at least `min_occurrence` times in the dataset. This is done using a set comprehension that iterates over the items in the `sequence_occurrences` Counter.

Finally, if there are any rescued sequences, the function checks if the `output_dir` directory exists and creates it if it doesn't. It then opens a new file in this directory, "RescuedChimera.fa", and writes the rescued sequences to this file in FASTA format. Each sequence is written as a `SeqRecord` object, with the sequence ID set to "RescuedChimera" and no description. If no sequences met the rescue criteria, a message is printed to the console.
