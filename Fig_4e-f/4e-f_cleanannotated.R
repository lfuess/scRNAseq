##this code creates figures 4e-f in our manuscript
##we plot expression trends following parasite infection for apc (e) and b-cell (f) markers
## for the infection expression we use data from Fuess et al. 2021 (mol ecol)

##package info##
library(janitor)
library(data.table)
library(tidyverse)
library(plyr)
library(tidyr)
library(Rmisc)
library(biomaRt)
library(dplyr)

##we also need to set our biomart info so we can translate transcript ids to gene ids##
bm <- useMart("ensembl")
bm <- useDataset("gaculeatus_gene_ensembl", mart=bm)

##and let's create a list of our significant infection DEGs from our data set##
  genes = read.csv("DESeq_Infection_Full_StringParam_new.csv")
  genes[1:10,1:7]
  ##select sig genes
  genes = genes[genes$padj < .1,]
  ##switch transcript IDs to gene names##
  gene_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id_version',"ensembl_gene_id"),filters = "ensembl_transcript_id_version", values = genes$X, bmHeader = T, )
  genes_final = merge(genes, gene_names, by.x = "X", by.y = "Transcript stable ID version")
  ##create final matrix
  genes_final[1:10,1:8]
  genes_final=genes_final[,c(1,8,7)]


##we'll start with our APC graph##
  ##create a list of marker genes that match significant DEGs from our data set##
    ##first generate our list of singificant markers##
    APC = read.csv("APC_MarkerData.csv")
    APC_markers = APC[APC$Antigen.Presenting.Myeloid.Cells.P.Value < .1,]
    APC_markers=APC_markers[,c(1:2)]
    ##and merge our sig genes with this##
    sig_APC=merge(genes_final, APC_markers, by.x = "Gene stable ID", by.y = "FeatureID")
  ##now we have that list we can merge it with our normalized read data and create our graph##
    ##read in read data##
    reads = read.csv("normalizedreads_TagSeq.csv", check.names = FALSE)
    sig_APC[1:10,1:4]
    sig_APC = sig_APC[,c(2,4)]
    reads[1:10,1:10]
    names(reads)[names(reads) == ''] <- 'FeatureID'
    names(sig_APC)[names(sig_APC) == 'X'] <- 'FeatureID'
    
    ##create a matrix of read data for the genes of interest (neutrophil markers) ##
      matrix=join(sig_APC, reads, by = "FeatureID")
      matrix[1:10,1:10]
      matrix=matrix[,c(1, 3:392)]
      matrix=matrix %>% remove_rownames %>% column_to_rownames(var="FeatureID")
      ##flip matrix##
      matrix2 = as.data.frame(t(matrix))
      matrix2[1:10,1:10]
      ##convert row names to first column##
      matrix2= setDT(matrix2, keep.rownames = TRUE)[]
      matrix2[1:10,1:10]
    
    ##manipulate the data to get it in the format me we need it for our graphs (pivot the matrix)##
      longmatrix <- gather(matrix2,gene,value,-rn)
      longmatrix
    ##merge in infection info##
      infect = read.csv("ExpDesign.csv")
      infect[1:10,1:8]
      infect=infect[,c(1,6)]
      combined_matrix = merge(longmatrix, infect, by.x = "rn", by.y = "Sample")
      combined_matrix[1:10,1:4]  
    ##now we can create summary stats for each gene and condition (infected/uninfected##
      sum = summarySE(combined_matrix, measurevar="value", groupvars=c("gene", "worm_present"))
    ##and then we plot our graphs##
      ##relevel infection data#
      sum$worm_present <- factor(sum$worm_present, levels = c("FALSE", "TRUE"))
      p_APC=ggplot(sum, aes(x=worm_present, y=value, group=gene)) + geom_line() + geom_point() +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black")) + 
        xlab("Infection Status") + ylab("Normalized Expression")
      
      
##next is B cell markers##
      ##create a list of marker genes that match significant DEGs from our data set##
      ##first generate our list of singificant markers##
      Bcell = read.csv("Bcell_MarkerData.csv")
      Bcell_markers = Bcell[Bcell$B.Cells.P.Value < .1,]
      Bcell_markers=Bcell_markers[,c(1:2)]
      ##and merge our sig genes with this##
      sig_Bcell=merge(genes_final, Bcell_markers, by.x = "Gene stable ID", by.y = "FeatureID")
      ##now we have that list we can merge it with our normalized read data and create our graph##
      ##read in read data##
      reads = read.csv("normalizedreads_TagSeq.csv", check.names = FALSE)
      sig_Bcell[1:10,1:4]
      sig_Bcell = sig_Bcell[,c(2,4)]
      reads[1:10,1:10]
      names(reads)[names(reads) == ''] <- 'FeatureID'
      names(sig_Bcell)[names(sig_Bcell) == 'X'] <- 'FeatureID'
      
      ##create a matrix of read data for the genes of interest (neutrophil markers) ##
      matrix_Bcell=join(sig_Bcell, reads, by = "FeatureID")
      matrix_Bcell[1:10,1:10]
      matrix_Bcell=matrix_Bcell[,c(1, 3:392)]
      matrix_Bcell=matrix_Bcell %>% remove_rownames %>% column_to_rownames(var="FeatureID")
      ##flip matrix##
      matrix_Bcell2 = as.data.frame(t(matrix_Bcell))
      matrix_Bcell2[1:10,1:10]
      ##convert row names to first column##
      matrix_Bcell2= setDT(matrix_Bcell2, keep.rownames = TRUE)[]
      matrix_Bcell2[1:10,1:10]
      
      ##manipulate the data to get it in the format me we need it for our graphs (pivot the matrix)##
      longmatrix_Bcell <- gather(matrix_Bcell2,gene,value,-rn)
      longmatrix_Bcell
      ##merge in infection info##
      pop = read.csv("ExpDesign.csv")
      pop[1:10,1:8]
      pop=pop[,c(1,3)]
      combined_matrix_Bcell = merge(longmatrix_Bcell, pop, by.x = "rn", by.y = "Sample")
      combined_matrix_Bcell[1:10,1:4]  
      ##now we can create summary stats for each gene and condition (infected/uninfected##
      sum_Bcell = summarySE(combined_matrix_Bcell, measurevar="value", groupvars=c("gene", "CrossDir"))
      ##remove our F2 stats##
      sum_Bcell=sum_Bcell[sum_Bcell$CrossDir != "F2", ]
      ##and then we plot our graphs##
      ##relevel infection data#
      p_Bcell=ggplot(sum_Bcell, aes(x=CrossDir, y=value, group=gene)) + geom_line() + geom_point() +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black")) + 
        xlab("Population") + ylab("Normalized Expression")
      
##make a little multiplot##
    multiplot(p_APC, p_Bcell, cols = 2)
  
    
  