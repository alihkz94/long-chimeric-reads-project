#!/bin/bash

#Cutadapt_ITS : 

cutadapt -a KCNGTWGGWGAACCWGC \
-g "CATATHANTAAGSSSAGG;min_overlap=16" \
-e 2 -j 8 \
--discard-untrimmed \
--action=retain \
-o ITS.fasta \
fungal_species.fasta > result_ITS.txt


#Cutadapt_SSU: 

cutadapt -a TACCTGGTTGATYCTGCCAGT \
-g "KCNGTWGGWGAACCWGC;min_overlap=15" \
-e 2 -j 8 \
--discard-untrimmed \
--action=retain \
-o SSU.fasta \
fungal_species.fasta > result_SSU.txt

