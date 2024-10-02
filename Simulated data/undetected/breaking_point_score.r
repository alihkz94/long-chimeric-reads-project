
### R ###
# Calculate the "breaking point scores"  - represents the percentage of the shorter chimeric fragment from the full chimering sequence, and the length of the shorter chimeric fragment

# OUTPUT is "breaking_point_scores.txt"
    #sequence name
    #percentage of the shorter chimeric fragment
    #shorter chimeric fragment length

library(seqinr)

fasta = read.fasta("chimeric.fasta", seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE, seqonly = FALSE)
seq_names = getName(fasta)
seqs = unlist(getSequence(fasta, as.string = TRUE))

#get breaking point value per seq
break_pts = gsub("_.*", "", (gsub(".*_at_", "", seq_names)))

#get sequence lengths
lens = nchar(seqs)

#check if output file existst and delete if it does. 
if (file.exists("breaking_point_scores.txt")) {
  file.remove("breaking_point_scores.txt")
}

# calculate breaking point scores and output a tab delimited file
for (i in 1:length(lens)) {
    part_a = as.numeric(lens[i])-as.numeric(break_pts[i])
    part_b = as.numeric(lens[i])-(as.numeric(lens[i])-as.numeric(break_pts[i]))

    if (part_a <= part_b) {
        score = round(as.numeric(part_a)/as.numeric(lens[i]), digits = 4)
        cat(seq_names[i], score, part_a, "\n", sep = "\t", file = "breaking_point_scores.txt", append = TRUE)

    } else if (part_b < part_a) {
        score = round(as.numeric(part_b)/as.numeric(lens[i]), digits = 4)
        cat(seq_names[i], score, part_b, "\n", sep = "\t", file = "breaking_point_scores.txt", append = TRUE)

    }
}

