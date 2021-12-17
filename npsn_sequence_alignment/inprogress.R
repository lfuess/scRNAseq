
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("msa")
##this file creates our supplementary figure which is a multiple sequence alignment of zebrafish and stickleback npsn##
##we use a fasta of the protein sequences of each of our npsn transcripts to align##
library(msa)

##load in sequence data##
  mySequenceFile = "npsn_proteins.fa"
  mySequences <- readAAStringSet(mySequenceFile)
  mySequences  

##run an alignment with default parameters
  myFirstAlignment <- msa(mySequences)
  print(myFirstAlignment, show="complete")

##make a pretty figure##
  msaPrettyPrint(myFirstAlignment, output = c("pdf"), showNames="none",
                 showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
  ?msaPrettyPrint
  
  