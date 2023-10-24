## OPTIMUM SETTINGS

The best settings for the DADA2 for simulated dataset was with : 

`list(minSampleFraction = 0.9, ignoreNNegatives = 1, minFoldParentOverAbundance = 1.5, 
minParentAbundance = 2, allowOneOff = TRUE, minOneOffParentDistance = 2, maxShift = 300),
`
## The results obtained based on these options:  
`True Positives: 623
False Positives: 6013
True Negatives: 38457
False Negatives: 148
Sensitivity: 0.808041504539559
Specificity: 0.8647852484821228
Precision: 0.0938818565400844
Recall: 0.808041504539559
F1 Score: 0.16821925205886326
Confusion Matrix:
[[38457  6013]
 [  148   623]]`
