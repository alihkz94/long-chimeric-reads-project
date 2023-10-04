in our tests we got this results so far: 

 For the consensus:
Total original sequences: 45263

Total filtered chimeras: 6039

Total chimera info records: 793

Debugging Info: TP+FP+TN+FN = 45263

Calculated Metrics: {'True Positives': 427, 'False Positives': 0, 'True Negatives': 44470, 'False Negatives': 366, 'Sensitivity': 0.5384615384615384, 'Specificity': 1.0}

 

for the pooled method: 
Total original sequences: 45263

Total filtered chimeras: 925

Total chimera info records: 793

Debugging Info: TP+FP+TN+FN = 45263

Calculated Metrics: {'True Positives': 293, 'False Positives': 0, 'True Negatives': 44470, 'False Negatives': 500, 'Sensitivity': 0.3694829760403531, 'Specificity': 1.0}

Moreover, in the case of Vsearch is also in this way: 

Calculated Metrics: {'True Positives': 142, 'False Positives': 0, 'True Negatives': 44470, 'False Negatives': 629, 'Sensitivity': 0.18417639429312582, 'Specificity': 1.0}

## Number of short reads

Number of chimeric reads with one part shorter than 60 bps: 53

Number of chimeric reads with one part shorter than 70 bps: 103

Number of chimeric reads with one part shorter than 100 bps: 283
