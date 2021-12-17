##this code makes figure 1b: a heatmap of the top marker genes for each of our 8 major cell types##
##load necessary packages##
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)

##first we need to create a list of our top 5 marker genes per cell type##
  ##isolate top sig genes (5 or all) for each cell group##
    ##hscs##
      hsc=read.csv("HSCs.csv")
      hsc=hsc[,c(2,5)]
      ##pick out the top five annotated genes##
      hsc=hsc[c(3:7),]
      colnames(hsc) <- c("gene","pval")
    ##neuts##
      neuts=read.csv("Neuts.csv")
      neuts=neuts[,c(2,5)]
      ##pick out the top five annotated genes##
      neuts=neuts[c(2,4,6:8),]
      colnames(neuts) <- c("gene","pval")
    ##apcs##
      apcs=read.csv("APCs.csv")
      apcs=apcs[,c(2,5)]
      ##pick out the top five annotated genes##
      apcs=apcs[c(9:12,15),]
      colnames(apcs) <- c("gene","pval")
    ##bcells##
      bcells=read.csv("Bcells.csv")
      bcells=bcells[,c(2,5)]
      ##pick out the top five annotated genes##
      bcells=bcells[c(4:5,7,9,14),]
      colnames(bcells) <- c("gene","pval")
    ##rbcs##
      rbcs=read.csv("RBCs.csv")
      rbcs=rbcs[,c(2,5)]
      ##pick out the top five annotated genes##
      rbcs=rbcs[c(3,5,7,10,12),]
      colnames(rbcs) <- c("gene","pval")
    ##platlets##
      plats=read.csv("Plat.csv")
      plats=plats[,c(2,5)]
      ##pick out the top five annotated genes##
      plats=plats[c(2,4:7),]
      colnames(plats) <- c("gene","pval")
    ##fibroblasts##
      fib=read.csv("fib.csv")
      fib=fib[,c(2,5)]
      ##pick out the top five annotated genes##
      fib=fib[c(1:4,6),]
      colnames(fib) <- c("gene","pval")
    ##nkcs##
      nkcs=read.csv("NKCs.csv")
      nkcs=nkcs[,c(2,5)]
      ##pick out the top sig annotated genes##
      nkcs=nkcs[c(3:4,7),]
      colnames(nkcs) <- c("gene","pval")

  ##now we bind all those together to create a list of interesting features##
    genes=rbind(hsc,neuts,apcs,bcells,rbcs,plats,fib,nkcs)
    features=genes[,c(1)]    
    
##Next thing we need to do is import our metadata which we'll use to filter our cells##
    cluster_info = read.csv("ClusterInfo_noUnk.csv")
    cluster_info[cluster_info == "HSCs"] <- "HCs"
    
##Now we'll actually start inputing and formatting our data##
    ##input raw data##
    hk.data <- Read10X(data.dir = "./filtered_feature_bc_matrix/")
    hk <- CreateSeuratObject(counts = hk.data, project = "hk", min.cells = 0, min.features = 0)
    
    ##subset the data to remove the cells that were filtered or belong to "unknown" clusters##
    cells.use <- cluster_info$Barcode
    good_hk <- subset(hk, cells = cells.use)
    good_hk
    
    ##add metadata (from Loupe Browser)##
    good_hk[["Clusters"]] = cluster_info$Simplified
    
    ##set clusters based on Loupe Browser analysis##
    Idents(object = good_hk) <- good_hk@meta.data$'Clusters'
    
    ##scaling for figure making##
    all.genes <- rownames(good_hk)
    good_hk <- ScaleData(good_hk, features = all.genes)

##And finally we make our heatmap using the filtered, scaled, and categorized raw data + our list of interesting features##
    my_levels <- c("HCs", "Neutrophils", "APCs",
                   "B-Cells","Erythrocytes","Platelets","Fibroblasts","NKCs")
    good_hk@active.ident  = factor(x = good_hk@active.ident, levels = my_levels)
    colors = c("#440154","#46337E","#365C8D","#277F8E","#1FA187","#4AC16D","#9FDA3A",
               "#FDE725")
    DoHeatmap(subset(good_hk, downsample = 100), features = features, size = 3, group.colors = colors) + scale_fill_viridis()
    