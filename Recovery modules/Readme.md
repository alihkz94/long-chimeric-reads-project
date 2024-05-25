## RECOVERY MODULES FOR THE CHIMERIC READS

### BLAST
#### Dependencies: Python packages: Biopython, pandas

To apply the BLAST module to the chimeric queries, the *BLASTn.sh* module must run first on them; later, the *BLAST_recovery.py* module can rescue the nonchimeric reads. The module can be run by the command below specifying the directories related to each output:

```bash
python BLAST_recovery.py --chimeras_dir ./chimeras --blast_output_dir ./blast_output --nonchimeric_dir .
```

### ReChime
Please read the ReChime_RUNME.txt file for the instructions in the **ReChime_v1.zip** file
