##code to make figure 1b for the scRNAseq data; we use Seurat's visualization components only##
##clusters IDs are set manually using output from Loupe Browser##

##load Seurat##
##remove.packages(grep("spatstat", installed.packages(), value = T))
##.rs.restartR()
##devtools::install_version("spatstat", version = "1.64-1")
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)

##first lets import metadata##
cluster_info = read.csv("ClusterInfo_noUnk.csv")
pop_info = read.csv("PopID.csv")
pop_info_goods = merge(cluster_info, pop_info, by = "Barcode")
pop_info_goods = pop_info_goods[,c(1,3)]

##input raw data##
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

##set features- the top 5 marker genes for each of our 8 main clusters##
features = c("stmn1b", "h2afva", "pcna", "hmgb2a","dkc1","alox5ap","alox5a","fbp2", "npsn", "mpx",
             "krt18b","krt8", "tmem176","cxcl19","sprb","swap70a","cd79a", "irf8", "dap1b","cd74a",
             "aqp1a.1","alas2","hbae5","HBE1","tmem205","tagln","boka","myh11a","rbpms", "TPM4",
             "col1a2","dcn","pcolcea","col1a1a","rbp4","sh2d1ab","nkl.4","id3")

##Heatmap- what we actually use for the figure##
my_levels <- c("HSCs", "Neutrophils", "APCs",
               "B-Cells","Erythrocytes","Platelets","Fibroblasts","NKCs")
good_hk@active.ident  = factor(x = good_hk@active.ident, levels = my_levels)
colors = c("#440154","#46337E","#365C8D","#277F8E","#1FA187","#4AC16D","#9FDA3A",
           "#FDE725")
DoHeatmap(subset(good_hk, downsample = 100), features = features, size = 3, group.colors = colors) + scale_fill_viridis()


## other graphs!!##
RidgePlot(good_hk, features = features, ncol = 2)
VlnPlot(good_hk, features = features)
DotPlot(good_hk, features = features, col.max = 10, cols = viridis(1000)) + RotatedAxis()


