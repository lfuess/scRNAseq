##this code takes data from the RBC contrast in Loupe Browser and converts it into the heatmap in Figure 2a##
##packages and stuff##
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(Matrix)
library(reshape2)
library(data.table)
library(tidyverse)
library(vioplot)
library(RColorBrewer)
library(pheatmap)

##first we need to get our list of significantly differentially expressed genes##
  ##read in the data##
  genes=read.csv("RBC_contrast.csv")
  ##filter out non-sig genes##
  siggenes=genes[genes$Cluster.22.P.Value < .1, ]
  ##sort data frame to help remove unannotated genes##
  sort_siggenes=siggenes[order(siggenes$FeatureName),] 
  ##remove rows with no gene names, mitochondrial genes, ribosomal genes##
  genes=sort_siggenes[c(1:36,77:91,99:113,145:146, 149:159),]
  ##reduce##
  genes=genes[,c(1:2,4)]
  colnames(genes) <- c("gene","gene_name", "lfc")
  ##and sort by p-val again##
  genes_final=genes[order(genes$lfc),]
  
  ##then we need to filter them to just look at genes that are markers of Bcells or Neutrophils##
    bcells = read.csv("Bcell_markers.csv")
    sig_bs=bcells[bcells$B.Cells.P.Value < .1, ]
    sig_bs=sig_bs[,c(1:2)]
    colnames(sig_bs) <- c("gene","gene_name")
    neuts=read.csv("Neut_Markers.csv")
    sig_neuts=neuts[neuts$Neutrophils.P.Value < .1, ]
    sig_neuts=sig_neuts[,c(1:2)]
    colnames(sig_neuts) <- c("gene","gene_name")
    marks=rbind(sig_bs,sig_neuts)
    
  ##filter out genes##
    features = merge(marks, genes_final, by="gene")
    ##and sort by p-val again##
    features=features[order(features$lfc),]
    ##fix our duplicate npsns##
    features_Seurat=features[,c(2)] 
   
##read in the data##   
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
    test_log=log2(test4)
    
    ##stuff for heatmaps##
    heat.data = as.data.frame(test_log)
    heat.data = setDT(heat.data, keep.rownames = TRUE)[]
    ##merege in cluster info (for ordering)##
    clusters = read.csv("RBC_Clusts.csv",stringsAsFactors = FALSE)
    clusters$Custom.Cluster[clusters$Custom.Cluster=="Cluster 21"] <- "Neutrophil-like"
    clusters$Custom.Cluster[clusters$Custom.Cluster=="Cluster 22"] <- "Bcell-like"
    combined = merge(heat.data, clusters, by.x = "rn", by.y = "Barcode")
    
    
    ##order by Cluster##
    combined = combined[order(combined$Custom.Cluster),]
    combined[1:10,13301]
    ##set rownames##
    combined_final=combined %>% remove_rownames %>% column_to_rownames(var="rn")
    combined_final = combined_final[,1:13299]
    ##reformat the data by transposing##
    heatmap = as.data.frame(t(combined_final))
    heatmap=setDT(heatmap, keep.rownames = TRUE)[]
    ##pull out info for just HC Genes##
    ##rename duplicate npsns in features list##
    ##rename duplicate npsns##
    features[63,2]="npsn(A)"
    features[59,2]="npsn(B)"
    
    topcombined = merge(features, heatmap, by.x = "gene", by.y = "rn")
    colnames(topcombined)
    dim(topcombined)
    heatmap[1:10,1:10]
    topcombined[1:10,1:10]
    heatmapfinal=topcombined[,c(2, 5:5433)]
    realheatmap = heatmapfinal %>% remove_rownames %>% column_to_rownames(var="gene_name.x")
    realheatmap[is.na(realheatmap)] <- 0
    
##make heatmap
    ##and put in cluster annotations##
    colnames(clusters) <- c("Barcode","Cluster")
    clusters$Cluster <- factor(clusters$Cluster, levels = c("Neutrophil-like", "Bcell-like"))
    annos=clusters %>% remove_rownames %>% column_to_rownames(var="Barcode")
    annotation_colors = list(
      Cluster = c("Neutrophil-like" = "#35aa93", "Bcell-like" = "#8fd0c3"))
    
    #clusters$Cluster = factor(clusters$Cluster)
    pheatmap(realheatmap, color = viridis(500), cluster_cols = FALSE,
             show_colnames = FALSE, scale = "none", annotation_col = annos, annotation_colors = annotation_colors,
             annotation_legend = TRUE)
    
  


  



  
  
  
##If we want to make it with Seurat##
  ##We need to process cluster IDs@@
    clust=read.csv("RBC_clusts.csv")
    clust$Custom.Cluster[clust$Custom.Cluster=="Cluster 21"] <- "Neutrophil-like"
    clust$Custom.Cluster[clust$Custom.Cluster=="Cluster 22"] <- "Bcell-like"
  ##input raw data##
  hk.data <- Read10X(data.dir = "./filtered_feature_bc_matrix/")
  hk <- CreateSeuratObject(counts = hk.data, project = "hk", min.cells = 0, min.features = 0)
  
  ##subset the data to remove the cells that were filtered or belong to "unknown" clusters##
  cells.use <- clust$Barcode
  good_hk <- subset(hk, cells = cells.use)
  good_hk
  
  ##add metadata (from Loupe Browser)##
  good_hk[["Clusters"]] = clust$Custom.Cluster

  ##set clusters based on Loupe Browser analysis##
  Idents(object = good_hk) <- good_hk@meta.data$'Clusters'  

  ##scaling for figure making##
  all.genes <- rownames(good_hk)
  good_hk <- ScaleData(good_hk, features = all.genes)  
  
  ##And finally we make our heatmap using the filtered, scaled, and categorized raw data + our list of interesting features##
  my_levels <- c("Neutrophil-like", "Bcell-like")
  good_hk@active.ident  = factor(x = good_hk@active.ident, levels = my_levels)
  colors = c("#35aa93","#8fd0c3")
  DoHeatmap(subset(good_hk, downsample = 100), features = features, size = 3, group.colors = colors, slot = "scale.data") + scale_fill_viridis()  
  
  