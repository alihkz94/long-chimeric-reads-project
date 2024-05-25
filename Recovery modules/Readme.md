## RECOVERY MODULES FOR THE CHIMERIC READS

### BLAST
 For applying the BLAST module on the chimeric, read the *BLASTn.sh* module needs to run first on the chimeric reads, and later on, the *BLAST_recovery.py* module can be used to rescue the nonchimeric reads. The module can be run by the command below with specifying the directories related to each output:

```bash
python blast_recovery.py --chimeras_dir ./chimeras/ --blast_output_dir ./blast_output/ --nonchimeric_dir .
```

### ReChime
