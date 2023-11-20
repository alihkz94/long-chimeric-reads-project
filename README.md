# chimera
Scripts related to generating simulated reads and chimera tools benchmarking 

1. started by cutadapt with trimming ITS regions.
2. Later on, the "ITS.fasta" file was inserted for the simulation with SimLoRD.
3. Later on, quality filtering was applied to reads with the PipeCarft-VSEARCH module, and the script was: 
``` bash
vsearch --fastq_filter input_file --fastq_maxee 1 --fastq_maxns 0  --fastq_minlen 50 --threads 8 --fastq_qmax 93 --fastq_qmin 0   --fastqout /input/qualFiltered_out/output_file.fastq
``` 
4. The fastq reads converted to the Fasta ones.
5.  Later on ITSx was applied to the data.
6.  After the ITSx, The full and partial one was selected from ITSx outputs.
7.  Chimera simulator will be applied to get chimeric read out of these Fasta reads.
8.  After this chimeric reads were generated they were concatenated into one file.
9.  UCHIME reference-based chimera filtering ran over the files. 


## Best settings 
The best settings that Uchime achieved were with settings --mindiv 0.4 --dn 1.6 --minh 0.08. 

## treat chimeras

we can run it by this script: 
``` bash
python chimera_recovery.py --blast_output blast_1st_tophit.txt --output recovered.csv --min_occurrence 2 --input_fasta combined_chimeras.fasta
``` 
