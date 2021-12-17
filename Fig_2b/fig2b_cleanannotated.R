##this is the code we will use to create the violin plots featured in figure 2b of our manuscript##
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
  clusters=read.csv("RBC_clusts.csv",stringsAsFactors = FALSE)
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
  
  ##and the cluster info##
  clusters = read.csv("RBC_clusts.csv",stringsAsFactors = FALSE)
  clusters$Custom.Cluster[clusters$Custom.Cluster=="Cluster 21"] <- "Neutrophil-like"
  clusters$Custom.Cluster[clusters$Custom.Cluster=="Cluster 22"] <- "Bcell-like"
  colnames(clusters) <- c("Barcode","Cluster")
  ##make the final combined data frame for a violin plot##
  vcombined = merge(violin.data, clusters, by.x = "rn", by.y = "Barcode")
  vcombined$Cluster = factor(vcombined$Cluster)
  vcombined[is.na(vcombined)] <- 0

##And actually plot!##
  ##b-cell markers##
    ##swap70a##
    par(cex.lab=1.1, cex.axis=.8) 
    vioplot(vcombined$ENSGACG00000015680[vcombined$Cluster=="Neutrophil-like"],
            vcombined$ENSGACG00000015680[vcombined$Cluster=="Bcell-like"],
            names=c("Neutrophil like", "B-cell like"),
            col=c("#35aa93","#8fd0c3"),
            rectCol="grey", colMed = 'black',
            xlab = "RBC Subdivision",
            ylab = "Expression", areaEqual = T)
    ##dap1b
    par(cex.lab=1.1, cex.axis=.8) 
    vioplot(vcombined$ENSGACG00000005778[vcombined$Cluster=="Neutrophil-like"],
            vcombined$ENSGACG00000005778[vcombined$Cluster=="Bcell-like"],
            names=c("Neutrophil like", "B-cell like"),
            col=c("#35aa93","#8fd0c3"),
            rectCol="grey", colMed = 'black',
            xlab = "RBC Subdivision",
            ylab = "Expression", areaEqual = T)
    ##cd74a
    par(cex.lab=1.1, cex.axis=.8) 
    vioplot(vcombined$ENSGACG00000018016[vcombined$Cluster=="Neutrophil-like"],
            vcombined$ENSGACG00000018016[vcombined$Cluster=="Bcell-like"],
            names=c("Neutrophil like", "B-cell like"),
            col=c("#35aa93","#8fd0c3"),
            rectCol="grey", colMed = 'black',
            xlab = "RBC Subdivision",
            ylab = "Expression", areaEqual = T)
  ##neutrophil markers##
    ##npsn(A)
    par(cex.lab=1.1, cex.axis=.8) 
    vioplot(vcombined$ENSGACG00000000834[vcombined$Cluster=="Neutrophil-like"],
            vcombined$ENSGACG00000000834[vcombined$Cluster=="Bcell-like"],
            names=c("Neutrophil like", "B-cell like"),
            col=c("#35aa93","#8fd0c3"),
            rectCol="grey", colMed = 'black',
            xlab = "RBC Subdivision",
            ylab = "Expression", areaEqual = T)
    ##alox5ap
    par(cex.lab=1.1, cex.axis=.8) 
    vioplot(vcombined$ENSGACG00000020547[vcombined$Cluster=="Neutrophil-like"],
            vcombined$ENSGACG00000020547[vcombined$Cluster=="Bcell-like"],
            names=c("Neutrophil like", "B-cell like"),
            col=c("#35aa93","#8fd0c3"),
            rectCol="grey", colMed = 'black',
            xlab = "RBC Subdivision",
            ylab = "Expression", areaEqual = T)
    ##fbp2
    par(cex.lab=1.1, cex.axis=.8) 
    vioplot(vcombined$ENSGACG00000014527[vcombined$Cluster=="Neutrophil-like"],
            vcombined$ENSGACG00000014527[vcombined$Cluster=="Bcell-like"],
            names=c("Neutrophil like", "B-cell like"),
            col=c("#35aa93","#8fd0c3"),
            rectCol="grey", colMed = 'black',
            xlab = "RBC Subdivision",
            ylab = "Expression", areaEqual = T)
    ##pygl
    par(cex.lab=1.1, cex.axis=.8) 
    vioplot(vcombined$ENSGACG00000011880[vcombined$Cluster=="Neutrophil-like"],
            vcombined$ENSGACG00000011880[vcombined$Cluster=="Bcell-like"],
            names=c("Neutrophil like", "B-cell like"),
            col=c("#35aa93","#8fd0c3"),
            rectCol="grey", colMed = 'black',
            xlab = "RBC Subdivision",
            ylab = "Expression", areaEqual = T)
    ##npsn(B)
    par(cex.lab=1.1, cex.axis=.8) 
    vioplot(vcombined$ENSGACG00000014415[vcombined$Cluster=="Neutrophil-like"],
            vcombined$ENSGACG00000014415[vcombined$Cluster=="Bcell-like"],
            names=c("Neutrophil like", "B-cell like"),
            col=c("#35aa93","#8fd0c3"),
            rectCol="grey", colMed = 'black',
            xlab = "RBC Subdivision",
            ylab = "Expression", areaEqual = T)
    ##mpx
    par(cex.lab=1.1, cex.axis=.8) 
    vioplot(vcombined$ENSGACG00000017166[vcombined$Cluster=="Neutrophil-like"],
            vcombined$ENSGACG00000017166[vcombined$Cluster=="Bcell-like"],
            names=c("Neutrophil like", "B-cell like"),
            col=c("#35aa93","#8fd0c3"),
            rectCol="grey", colMed = 'black',
            xlab = "RBC Subdivision",
            ylab = "Expression", areaEqual = T)
    
  