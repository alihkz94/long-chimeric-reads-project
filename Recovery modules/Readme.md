## RECOVERY MODULES FOR THE CHIMERIC READS

### BLAST
#### Dependencies:
#### Python packages: Biopython, pandas

To apply the BLAST module to the chimeric queries, the *BLASTn.sh* module needs to run first on them; later on, the *BLAST_recovery.py* module can rescue the nonchimeric reads. The module can be run by the command below specifying the directories related to each output:

```bash
python BLAST_recovery.py --chimeras_dir ./chimeras --blast_output_dir ./blast_output --nonchimeric_dir .
```

### ReChime
