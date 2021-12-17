##we're going to make a supplemental file that displays patterns of significant differential expression between our two B-cell subsets##
library(Matrix)
library(reshape2)
library(data.table)
library(dplyr)
library(tidyverse)
library(vioplot)
library(RColorBrewer)
library(viridis)
library(pheatmap)
##first we need to get our list of significantly differentially expressed genes##
  ##read in each of the datasets##
  genes=read.csv("Bcell_cluster_DEGs.csv")
  ##filter out non-sig genes##
  siggenes=genes[genes$Cluster.13.P.Value < .1, ]
  ##sort data frame to help remove unannotated genes##
  sort_siggenes=siggenes[order(siggenes$FeatureName),] 
  ##reduce##
  finalgenes=sort_siggenes[,c(1:2,4)]
  colnames(finalgenes) <- c("gene","gene_name", "lfc")
  
##Next we need to get and filter our raw data##
  ##get our barcodes of interest and order them##
  cluster_info = read.csv("BcellSubclusters.csv")
  colnames(cluster_info) <- c("Barcode","Cluster")
  cluster_info=cluster_info %>% arrange(factor(Cluster))
  
  ##match with dataset##
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
  barcodes=as.vector(cluster_info$Barcode)
  test1= filter(test, Cell %in% barcodes)  
  test2=as.data.table(test1)
  ##pivot the matrix##
  test3 = data.table::dcast(test2, Cell ~ Gene,  value.var = "Count")
  
  test4=test3 %>% remove_rownames %>% column_to_rownames(var="Cell")
  test_log=log2(test4)
  heat.data = as.data.frame(test_log)
  heat.data = setDT(heat.data, keep.rownames = TRUE)[]
  ##merege in population info (for ordering)##
  combined = merge(cluster_info, heat.data, by.x = "Barcode", by.y = "rn")
  combined[1:10,1:10]  
  ##reorder##
  ordered=combined %>% arrange(factor(Cluster))
  ##set rownames##
  ordered_final=ordered %>% remove_rownames %>% column_to_rownames(var="Barcode")
  ordered_final = ordered_final[,1:13993]

##make a data frame with just the info for the genes were interested in##
  ##reformat the data by transposing##
  heatmap = as.data.frame(t(ordered_final))
  sapply(heatmap, class)
  heatmap=setDT(heatmap, keep.rownames = TRUE)[]
  ##pull out just the info for genes of interest##
  topcombined = merge(finalgenes, heatmap, by.x = "gene", by.y = "rn")
  colnames(topcombined)
  dim(topcombined)
  topcombined[1:10,1:10]
  heatmapfinal=topcombined[,c(2,4:12920)]
  realheatmap = heatmapfinal %>% remove_rownames %>% column_to_rownames(var="gene_name")
  realheatmap[is.na(realheatmap)] <- 0
  ##make sure it's numeric##
  realheatmap2 <- mutate_all(realheatmap, function(x) as.numeric(as.character(x)))
  
##add in cluster annotations##
  annos=cluster_info %>% remove_rownames %>% column_to_rownames(var="Barcode")
  annos$Cluster <- factor(annos$Cluster, levels = c("Cluster 12", "Cluster 13"))
  annotation_colors = list(Cluster = c("Cluster 12" = "#7db2bb", "Cluster 13" = "#d4e5e8"))
  
#and finally make the heatmap
  pheatmap(realheatmap2, color = viridis(500), cluster_cols = FALSE,
           show_colnames = FALSE, scale = "none", annotation_col = annos, annotation_colors = annotation_colors,
           annotation_legend = TRUE)  
  