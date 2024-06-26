# Data and scripts for "Enhancing long-read amplicon sequencing: Overcoming chimeric sequence challenges in biodiversity studies with VSEARCH and DADA2" (Hakimzadeh et al. 2024)

### Structure
This repository contains the data and part of the analysis stack for the abovementioned paper. It is structured as follows:

[Chimeras_denovo](https://github.com/alihkz94/long-chimeric-reads-project/tree/main/Chimeras_denovo) holds scripts related to the long read module of the VSEARCH applied for the simulated dataset and statistical part for calculating the F1 score. 

[DADA2](https://github.com/alihkz94/long-chimeric-reads-project/tree/main/DADA2) holds scripts related to DADA2 quality filtering and chimera filtering related to the real dataset. Moreover the scripts for the simulated dataset and statical analysis for calculating F1 score. 

[Recovery modules](https://github.com/alihkz94/long-chimeric-reads-project/tree/main/Recovery%20modules) contains the BLAST scripts for alignment and recovery module for it. Besides it contains ReChime (v1) module. 

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
