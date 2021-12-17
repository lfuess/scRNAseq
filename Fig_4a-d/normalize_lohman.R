library("DESeq2")
colData = read.csv("Lohman_Traits.csv")
names(colData)
colData = colData[,c(1,24)]

countData <- read.csv("lohman_counts.csv", row.names="", check.names= FALSE)

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Exposed) ##doesn't matter##
dds <- estimateSizeFactors(dds) 
dds <- estimateDispersions(dds)

vst <- getVarianceStabilizedData(dds)

write.csv(vst, file = "normalizedreads_lohman.csv")
