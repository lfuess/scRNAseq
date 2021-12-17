##this is the code we will use to create the violin plots featured in figure 3b of our manuscript##
##the genes selected were hand currated, so here we'll parse the raw data matrix, then make individual graphs
##for each gene of interest
library(Matrix)
library(reshape2)
library(data.table)
library(dplyr)
library(tidyverse)
library(vioplot)
library(RColorBrewer)
library(viridis)

##process the raw data##
  ##create the sparse matrix##
  matrix_dir = "./filtered_feature_bc_matrix/"
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  summ <- summary(mat)
  
  ##convert to dense matrix##
  test = data.frame(Gene      = rownames(mat)[summ$i],
                    Cell = colnames(mat)[summ$j],
                    Count      = summ$x)
  
  ##free up some space##
  rm(mat)
  rm(summ)
  
  ##pull out just those that match our clustered cells##
  clusters=read.csv("Bcells_popinfo.csv",stringsAsFactors = FALSE)
  barcodes=as.vector(clusters$Barcode)
  test1= filter(test, Cell %in% barcodes)
  test2=as.data.table(test1)
  ##pivot the matrix##
  test3 = data.table::dcast(test2, Cell ~ Gene,  value.var = "Count")
  
  test4=test3 %>% remove_rownames %>% column_to_rownames(var="Cell")
  ##log transform it##
  test_log=log2(test4)
  
##manipulate the data to format for violin plotting##
  ##data frame manipulation##
  violin.data = as.data.frame(test_log)
  violin.data = setDT(violin.data, keep.rownames = TRUE)[]
  names(violin.data)
  
  ##and the population info##
  populations = read.csv("Bcells_popinfo.csv",stringsAsFactors = FALSE)
  colnames(populations) <- c("Barcode","Population")
  ##make the final combined data frame for a violin plot##
  vcombined = merge(violin.data, populations, by.x = "rn", by.y = "Barcode")
  vcombined$Population = factor(vcombined$Population)
  vcombined[is.na(vcombined)] <- 0  
  ##set population levels
  vcombined$Population = factor(vcombined$Population, levels = c("SAY", "GOS", "ROB"))

##And actually plot!##
  
  ##ENSGACG00000012769 (Ig Mu heavy chain)
  par(cex.lab=1.1, cex.axis=.8) 
  vioplot(vcombined$ENSGACG00000012769[vcombined$Pop=="SAY"],
          vcombined$ENSGACG00000012769[vcombined$Pop=="GOS"],
          vcombined$ENSGACG00000012769[vcombined$Pop=="ROB"],
          names=c("SAY", "GOS", "ROB"),
          col=c("#46337E","#4AC16D","#FDE725"),
          rectCol="grey", colMed = 'black',
          xlab = "Cell Type",
          ylab = "Expression", areaEqual = T)
  
  ##ENSGACG00000012783 (Ig Mu heavy chain)
  par(cex.lab=1.1, cex.axis=.8) 
  vioplot(vcombined$ENSGACG00000012783[vcombined$Pop=="SAY"],
          vcombined$ENSGACG00000012783[vcombined$Pop=="GOS"],
          vcombined$ENSGACG00000012783[vcombined$Pop=="ROB"],
          names=c("SAY", "GOS", "ROB"),
          col=c("#46337E","#4AC16D","#FDE725"),
          rectCol="grey", colMed = 'black',
          xlab = "Cell Type",
          ylab = "Expression", areaEqual = T)
  
  ##ENSGACG00000009315 (Ig kappa chain C region)
  par(cex.lab=1.1, cex.axis=.8) 
  vioplot(vcombined$ENSGACG00000009315[vcombined$Pop=="SAY"],
          vcombined$ENSGACG00000009315[vcombined$Pop=="GOS"],
          vcombined$ENSGACG00000009315[vcombined$Pop=="ROB"],
          names=c("SAY", "GOS", "ROB"),
          col=c("#46337E","#4AC16D","#FDE725"),
          rectCol="grey", colMed = 'black',
          xlab = "Cell Type",
          ylab = "Expression", areaEqual = T)
  
  ##ENSGACG00000009519 (Ig kappa chain C region)
  par(cex.lab=1.1, cex.axis=.8) 
  vioplot(vcombined$ENSGACG00000009519[vcombined$Pop=="SAY"],
          vcombined$ENSGACG00000009519[vcombined$Pop=="GOS"],
          vcombined$ENSGACG00000009519[vcombined$Pop=="ROB"],
          names=c("SAY", "GOS", "ROB"),
          col=c("#46337E","#4AC16D","#FDE725"),
          rectCol="grey", colMed = 'black',
          xlab = "Cell Type",
          ylab = "Expression", areaEqual = T)
  
  ##ENSGACG00000010720 (Ig lambda-3 chain C region)
  par(cex.lab=1.1, cex.axis=.8) 
  vioplot(vcombined$ENSGACG00000010720[vcombined$Pop=="SAY"],
          vcombined$ENSGACG00000010720[vcombined$Pop=="GOS"],
          vcombined$ENSGACG00000010720[vcombined$Pop=="ROB"],
          names=c("SAY", "GOS", "ROB"),
          col=c("#46337E","#4AC16D","#FDE725"),
          rectCol="grey", colMed = 'black',
          xlab = "Cell Type",
          ylab = "Expression", areaEqual = T)
  