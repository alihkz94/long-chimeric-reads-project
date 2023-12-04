## OPTIMUM SETTINGS

The best settings for the DADA2 for the simulated dataset were with : 

`list(minSampleFraction = 0.9, ignoreNNegatives = 1, minFoldParentOverAbundance = 1.5, 
minParentAbundance = 2, allowOneOff = TRUE, minOneOffParentDistance = 2, maxShift = 350),
`
## The results obtained based on these options:  

True Positives: 622

False Positives: 5986

True Negatives: 38484

False Negatives: 149

Sensitivity: 0.80674448767834

Specificity: 0.865392399370362 

Precision: 0.0941283292978208

Recall: 0.80674448767834 

F1 Score: 0.168586529340019



## get chimeric reads out:
asv_seqs = colnames(ASV_tab.nochim)
asv_headers = openssl::sha1(asv_seqs)


new_seqtab = so that headers are sha1
new_nonchimeric_seqtab = so that headers are sha1 

then compare those two tables to write out chimeras with sha1 headers

sha1 is HEADER:
SHA-1 (Secure Hash Algorithm 1) is a hash function which takes an input and produces a 160-bit (20-byte)
hash value known as a message digest – typically rendered as 40 hexadecimal digits.

OpenSSL library you can generate these asv headers. also, check the pipecraft headers to see how  it work.

asv_headers=openssl::sha1(asv_seqs)
