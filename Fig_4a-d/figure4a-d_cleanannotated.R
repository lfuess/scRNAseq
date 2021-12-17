##this code generates figures 4a-4d in our single cell RNAseq data##
##it looks back at old traditional transcriptomic data sets with accompanying flow cytometry data
##specifically creates correlations between flow cytometry data and expression of two key marker genes 

library(data.table)
library(tidyverse)
library(ggpubr)

##we start looking at correlations in Fuess et al. 2021
  ##first lets isolate our flow cytometry data##
  flowc = read.csv("Masterdata_TagSeq_SDH_23July2018.csv")
  names(flowc)
  flowc=flowc[,c(2,43:44)]
  flowc
  ##next bring in normalized reads and format##
  exp = read.csv("normalizedreads_tagseq.csv", check.names=FALSE)
  exp = as.data.frame(t(exp), stringsAsFactors = FALSE)
  exp[1:10,1:10]
  my.names <- exp[1,]
  
  colnames(exp) <- my.names
  exp = exp[-c(1), ] 
  dim(exp)
  
  
  exp=setDT(exp, keep.rownames = TRUE)[]
  
  ##merge them together##
  merged = merge(exp, flowc, by.x = "rn", by.y = "sample_ID")
  merged[1:10,1:10]
  merged = merged %>% remove_rownames %>% column_to_rownames(var="rn")
  
  merged2 = sapply(merged, as.numeric)
  merged2[1:10,1:10]
  rownames(merged2) = rownames (merged)
  final = as.data.frame(merged2)
  final[1:10,1:10]
  ##run some correlations##
  ##first normalize our frequencies##
  final$A_Freq.Gran_Norm = asin(sqrt(final$A_Freq.Gran/100))
  final$A_Freq.Lymph_Norm = asin(sqrt(final$A_Freq.Lymph/100))
  
  ##cd79a##
  p1=ggscatter(final, x="A_Freq.Lymph_Norm", y="ENSGACT00000004238.1", add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Frequency of Lymphocytes", ylab = "cd79a Expression")
  
  ##npsn b##
  p2 = ggscatter(final, x="A_Freq.Gran_Norm", y="ENSGACT00000019050.1", add = "reg.line", conf.int = TRUE,
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "Frequency of Granulocytes", ylab = "npsn.b Expression")


##Do this for lohman's data as well##

  ##first lets isolate our flow cytometry data##
  flowc = read.csv("Lohman_Traits.csv")
  names(flowc)
  flowc=flowc[,c(1,37:38)]
  flowc
  ##next bring in normalized reads and format##
  exp = read.csv("normalizedreads_lohman.csv", check.names=FALSE)
  exp = as.data.frame(t(exp), stringsAsFactors = FALSE)
  exp[1:10,1:10]
  my.names <- exp[1,]
  
  colnames(exp) <- my.names
  exp = exp[-c(1), ] 
  dim(exp)
  
  exp=setDT(exp, keep.rownames = TRUE)[]
  
  ##merge them together##
  merged = merge(exp, flowc, by.x = "rn", by.y = "SampleID")
  merged[1:10,1:10]
  merged = merged %>% remove_rownames %>% column_to_rownames(var="rn")
  
  merged2 = sapply(merged, as.numeric)
  merged2[1:10,1:10]
  rownames(merged2) = rownames (merged)
  final = as.data.frame(merged2)
  ##run some correlations##
  ##normalize##
  final$A_Freq.Gran_Norm = asin(sqrt(final$A_Freq.Gran/100))
  final$A_Freq.Lymph_Norm = asin(sqrt(final$A_Freq.Lymph/100))
  
  ##cd79a##
  cor.test(final$A_Freq.Lymph_Norm, final$ENSGACG00000003230)
  p3=ggscatter(final, x="A_Freq.Lymph_Norm", y="ENSGACG00000003230", add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Frequency of Lymphocytes", ylab = "cd79a Expression")
  
  ##test the npsn.B##
  cor.test(final$A_Freq.Gran_Norm, final$ENSGACG00000014415)
  p4=ggscatter(final, x="A_Freq.Gran_Norm", y="ENSGACG00000014415", add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Frequency of Granuloyctes", ylab = "npsn.b Expression")


##arrange the 4 plots for figure making##
multiplot(p1,p2,p3,p4,cols = 2)


