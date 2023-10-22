## OPTIMUM SETTINGS

The best settings for the DADA2 for simulated dataset was with : 

`list(minSampleFraction = 0.9, ignoreNNegatives = 1, minFoldParentOverAbundance = 1.5, 
minParentAbundance = 2, allowOneOff = TRUE, minOneOffParentDistance = 2, maxShift = 300),
`

If we decrease the `minOneOffParentDistance`  the number of the FALSE POSITIVE increase. it can be vital if we want to push more for detection of seqeunces. 
