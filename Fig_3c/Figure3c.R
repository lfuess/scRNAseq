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
  genes=read.csv("RBC_GvR.csv")
  ##filter out non-sig genes##
  siggenes=genes[genes$GOS.Log2.Fold.Change < .1, ]
  ##sort data frame to help remove unannotated genes##
  sort_siggenes=siggenes[order(siggenes$FeatureName),] 
  ##remove rows with no gene names##
  genes_RBC_GvR=sort_siggenes[c(1:3,7:15),]
  ##reduce##
  genes_RBC_GvR=genes_RBC_GvR[,c(1:2,4)]
  colnames(genes_RBC_GvR) <- c("gene","gene_name", "lfc")
  ##and sort by p-val again##
  genes_RBC_GvR=genes_RBC_GvR[order(genes_RBC_GvR$lfc),]
  
  genes=read.csv("RBC_GvS.csv")
  ##filter out non-sig genes##
  siggenes=genes[genes$GOS.Log2.Fold.Change < .1, ]
  ##sort data frame to help remove unannotated genes##
  sort_siggenes=siggenes[order(siggenes$FeatureName),] 
  ##remove rows with no gene names##
  genes_RBC_GvS=sort_siggenes[c(1:4,9:14),]
  ##reduce##
  genes_RBC_GvS=genes_RBC_GvS[,c(1:2,4)]
  colnames(genes_RBC_GvS) <- c("gene","gene_name", "lfc")
  ##and sort by p-val again##
  genes_RBC_GvS=genes_RBC_GvS[order(genes_RBC_GvS$lfc),]
  
  genes=read.csv("RBC_RvS.csv")
  ##filter out non-sig genes##
  siggenes=genes[genes$ROB.Log2.Fold.Change < .1, ]
  ##sort data frame to help remove unannotated genes##
  sort_siggenes=siggenes[order(siggenes$FeatureName),] 
  ##remove rows with no gene names##
  genes_RBC_RvS=sort_siggenes[c(4:5),]
  ##reduce##
  genes_RBC_RvS=genes_RBC_RvS[,c(1:2,4)]
  colnames(genes_RBC_RvS) <- c("gene","gene_name", "lfc")
  ##and sort by p-val again##
  genes_RBC_RvS=genes_RBC_RvS[order(genes_RBC_RvS$lfc),]
  
  genes=read.csv("HCs_GvR.csv")
  ##filter out non-sig genes##
  siggenes=genes[genes$GOS.Log2.Fold.Change < .1, ]
  ##sort data frame to help remove unannotated genes##
  sort_siggenes=siggenes[order(siggenes$FeatureName),] 
  ##remove rows with no gene names##
  genes_HC_GvR=sort_siggenes[c(1:23, 33:73,75:78),]
  ##reduce##
  genes_HC_GvR=genes_HC_GvR[,c(1:2,4)]
  colnames(genes_HC_GvR) <- c("gene","gene_name", "lfc")
  ##and sort by p-val again##
  genes_HC_GvR=genes_HC_GvR[order(genes_HC_GvR$lfc),]
  
  genes=read.csv("HCs_GvS.csv")
  ##filter out non-sig genes##
  siggenes=genes[genes$GOS.Log2.Fold.Change < .1, ]
  ##sort data frame to help remove unannotated genes##
  sort_siggenes=siggenes[order(siggenes$FeatureName),] 
  ##remove rows with no gene names##
  genes_HC_GvS=sort_siggenes[c(1:13,23:49),]
  ##reduce##
  genes_HC_GvS=genes_HC_GvS[,c(1:2,4)]
  colnames(genes_HC_GvS) <- c("gene","gene_name", "lfc")
  ##and sort by p-val again##
  genes_HC_GvS=genes_HC_GvS[order(genes_HC_GvS$lfc),]
  
  genes=read.csv("HCs_RvS.csv")
  ##filter out non-sig genes##
  siggenes=genes[genes$ROB.Log2.Fold.Change < .1, ]
  ##sort data frame to help remove unannotated genes##
  sort_siggenes=siggenes[order(siggenes$FeatureName),] 
  ##remove rows with no gene names##
  genes_HC_RvS=sort_siggenes[c(1:10,18:31),]
  ##reduce##
  genes_HC_RvS=genes_HC_RvS[,c(1:2,4)]
  colnames(genes_HC_RvS) <- c("gene","gene_name", "lfc")
  ##and sort by p-val again##
  genes_HC_RvS=genes_HC_RvS[order(genes_HC_RvS$lfc),]
  
  ##combine our list of filters together##
  all_sig_genes=rbind(genes_HC_GvR, genes_HC_GvS, genes_HC_RvS, genes_RBC_GvR, genes_RBC_GvS, genes_RBC_RvS)
  ##remove duplicates
  all_sig_genes = all_sig_genes[!duplicated(all_sig_genes$gene),]

  ##filter against neutrophil markers##
  neuts=read.csv("Neut_Markers.csv")
  sig_neuts=neuts[neuts$Neutrophils.P.Value < .1, ]
  sig_neuts=sig_neuts[,c(1:2)]
  colnames(sig_neuts) <- c("gene","gene_name")
  
  ##filter out genes##
  features = merge(sig_neuts, all_sig_genes, by="gene")

##Next we need to get and filter our raw data##
  ##get our barcodes of interest and order them##
  cluster_info = read.csv("Clust_IDs.csv")
  pop_info = read.csv("pops.csv")
  pop_info_goods = merge(cluster_info, pop_info, by = "Barcode")
  
  ##order it in the least elegant way possible (sorry)##
  HCs=pop_info_goods[pop_info_goods$Simplified == 'HCs' ,]
  target <- c("SAY", "GOS", "ROB")
  HCs=HCs %>% arrange(factor(sample_name, levels = target))
  RBCs=pop_info_goods[pop_info_goods$Simplified == 'RBCs' ,]
  target <- c("SAY", "GOS", "ROB")
  RBCs=RBCs %>% arrange(factor(sample_name, levels = target))
  pop_info_goods=rbind(HCs,RBCs)

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
    combined = merge(pop_info_goods, heat.data, by.x = "Barcode", by.y = "rn")
    combined[1:10,1:10]
    
    ##order it in the least elegant way possible (sorry)##
    HCs=combined[combined$Simplified == 'HCs' ,]
    target <- c("SAY", "GOS", "ROB")
    HCs=HCs %>% arrange(factor(sample_name, levels = target))
    RBCs=combined[combined$Simplified == 'RBCs' ,]
    target <- c("SAY", "GOS", "ROB")
    RBCs=RBCs %>% arrange(factor(sample_name, levels = target))
    ordered=rbind(HCs,RBCs)
    ordered[1:10,1:10]
    ##set rownames##
    combined_final=combined %>% remove_rownames %>% column_to_rownames(var="Barcode")
    combined_final = combined_final[,1:15433]
    
    ##make a top annotated gene heatmap##
    ##reformat the data by transposing##
    heatmap = as.data.frame(t(combined_final))
    sapply(heatmap, class)
    heatmap=setDT(heatmap, keep.rownames = TRUE)[]
    ##pull out just the info for genes of interest##
    topcombined = merge(features, heatmap, by.x = "gene", by.y = "rn")
    colnames(topcombined)
    dim(topcombined)
    realheatmap[1:10,1:10]
    heatmapfinal=topcombined[,c(2,5:12853)]
    realheatmap = heatmapfinal %>% remove_rownames %>% column_to_rownames(var="gene_name.x")
    realheatmap[is.na(realheatmap)] <- 0
    ##make sure it's numeric##
    realheatmap2 <- mutate_all(realheatmap, function(x) as.numeric(as.character(x)))
    
    ##and put in cluster annotations##
    colnames(pop_info_goods) <- c("Barcode","Type", "Pop")
    annos <- pop_info_goods[,-1]
    rownames(annos) <- pop_info_goods[,1]
    annos$Pop <- factor(annos$Pop, levels = c("SAY", "GOS", "ROB"))
    annos$Type<- factor(annos$Type, levels = c("HCs", "RBCs"))
    annotation_colors = list(Pop = c("SAY" = "#46337E", "GOS" = "#4AC16D", "ROB" = "#FDE725"), Type = c("HCs"="#440154", "RBCs" = "#1FA187"))
    #make the heatmap
    pheatmap(realheatmap2, color = viridis(500), cluster_cols = FALSE,
             show_colnames = FALSE, scale = "none", annotation_col = annos, annotation_colors = annotation_colors,
             annotation_legend = TRUE)
  


  
  
  
  
  
  
  
  
  
  
  

  ##Next thing we need to do is import our metadata which we'll use to filter our cells##
  cluster_info = read.csv("Clust_IDs.csv")  
  pop_info = read.csv("pops.csv")  
  pop_info_goods = merge(cluster_info, pop_info, by = "Barcode")
  pop_info_goods = pop_info_goods[,c(1,3)] 
  
##then input and filter raw data##
  hk.data <- Read10X(data.dir = "./filtered_feature_bc_matrix/")
  hk <- CreateSeuratObject(counts = hk.data, project = "hk", min.cells = 0, min.features = 0)
  
  ##subset the data to remove the cells that were filtered or belong to "unknown" clusters##
  cells.use <- cluster_info$Barcode
  good_hk <- subset(hk, cells = cells.use)
  good_hk
  
  ##add metadata (from Loupe Browser)##
  good_hk[["Clusters"]] = cluster_info$Simplified
  good_hk[["Population"]] = pop_info_goods$sample_name
  
  ##set clusters based on Loupe Browser analysis##
  Idents(object = good_hk) <- good_hk@meta.data$'Clusters'

  ##scaling for figure making##
  all.genes <- rownames(good_hk)
  good_hk <- ScaleData(good_hk, features = all.genes)  

##and graph it##
  good_hk@meta.data$Population <- factor(x = good_hk@meta.data$Population, levels=c('SAY', 'GOS', 'ROB'))
  
  DoMultiBarHeatmap(good_hk, features=features, group.by="Clusters", additional.group.by = "Population", slot = "scale.data") + scale_fill_viridis()
  
  