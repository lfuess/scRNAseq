##this code lays the framework for the section of our MS where we reinterpret old data set##
##here we generate a list of markers of each of our identified cell types
##and then we pull out expression data/data of interest for the expression of each of these markers in three data sets:
##Lohman 2017, frontiers in immuno
##Fuess 2021, molec ecol
##Fuess 2021, mBio

##package data
library(biomaRt)
library(dplyr)
# select mart and data set for biomart
bm <- useMart("ensembl")
bm <- useDataset("gaculeatus_gene_ensembl", mart=bm)

##first load in our marker lists##
HC=read.csv("HC_MarkerData.csv")
HC_markers = HC[HC$HCs.P.Value < .1,]
HC_markers$type = "HC"
HC_markers=HC_markers[,c(1:2,6)]

Neut = read.csv("Neutrophil_MarkerData.csv")
Neut_markers = Neut[Neut$Neutrophils.P.Value < .1,]
Neut_markers$type = "Neutrophil"
Neut_markers=Neut_markers[,c(1:2,6)]

APC = read.csv("APC_MarkerData.csv")
APC_markers = APC[APC$Antigen.Presenting.Myeloid.Cells.P.Value < .1,]
APC_markers$type = "APC"
APC_markers=APC_markers[,c(1:2,6)]

Bcell = read.csv("Bcell_MarkerData.csv")
Bcell_markers = Bcell[Bcell$B.Cells.P.Value < 0.1, ]
Bcell_markers$type = "Bcell"
Bcell_markers=Bcell_markers[,c(1:2,6)]

RBC = read.csv("RBC_MarkerData.csv")
RBC_markers = RBC[RBC$RBCs.P.Value < 0.1, ]
RBC_markers$type = "RBC"
RBC_markers=RBC_markers[,c(1:2,6)]

Plat = read.csv("Platelet_MarkerData.csv")
Plat_markers = Plat[Plat$Platelets.P.Value < 0.1, ]
Plat_markers$type = "Platelet"
Plat_markers=Plat_markers[,c(1:2,6)]

Fib = read.csv("Fibroblast_MarkerData.csv")
Fib_markers = Fib[Fib$Fibroblasts.P.Value < 0.1, ]
Fib_markers$type = "Fibroblast"
Fib_markers=Fib_markers[,c(1:2,6)]

NKC = read.csv("NKC_MarkerData.csv")
NKC_markers = NKC[NKC$Natural.Killer.Cells.P.Value < 0.1, ]
NKC_markers$type = "NKC"
NKC_markers=NKC_markers[,c(1:2,6)]

##merge them together##
markers = rbind(HC_markers,Neut_markers, APC_markers, Bcell_markers, RBC_markers, Plat_markers, Fib_markers, NKC_markers)

##now match these to old data sets##
##start with the Lohman data##

  ##first we need to isolate our three lists of differentially expressed genes##
  ##population
  ##Infectection
  ##popxInfectection
  Infect=read.csv("Control_vs_Infected_Expression.csv")
  Infect_sig <- Infect[Infect$padj < .1, ]
  Infect_sig=Infect_sig[rowSums(is.na(Infect_sig)) != ncol(Infect_sig), ]
  Infect_sig = Infect_sig[,c(1,3)]
  
  Int=read.csv("Population_InifectionStatus_Interaction_Expression.csv")
  Int_sig <- Int[Int$padj < .1,]
  Int_sig=Int_sig[rowSums(is.na(Int_sig)) != ncol(Int_sig), ]
  Int_sig = Int_sig[,c(1,3)]
  
  Pop=read.csv("Roberts_vs_Gosling.csv")
  Pop_sig <- Pop[Pop$padj < .1,]
  Pop_sig=Pop_sig[rowSums(is.na(Pop_sig)) != ncol(Pop_sig), ]
  Pop_sig = Pop_sig[,c(1,3)]
  
  ##Then we merge with our marker info##
  Infect_mark=merge(markers,Infect_sig, by.x = "FeatureID", by.y = "X")
  write.csv(Infect_mark, "Infection_markers_Lohman.csv", row.names = FALSE)
  Int_mark=merge(markers,Int_sig, by.x = "FeatureID", by.y = "X")
  write.csv(Int_mark, "Interaction_markers_Lohman.csv", row.names = FALSE)
  Pop_mark = merge(markers,Pop_sig, by.x = "FeatureID", by.y = "X")  
  write.csv(Pop_mark, "Population_markers_Lohman.csv", row.names = FALSE)

##next up is the results from my TagSeq2021 paper##
##we follow the same approach, but also convert transcript IDs to gene IDs
  ##RBCv.GBC##
  RvG_ts = read.csv("DESeq_RBCvGBC_Full_StringParam_new.csv")
  RvG_ts_sig = RvG_ts[RvG_ts$padj < .1,]  
  RvG_ts_sig = RvG_ts_sig[rowSums(is.na(RvG_ts_sig)) != ncol(RvG_ts_sig), ]
  RvG_ts_sig = RvG_ts_sig[,c(1,3)]
  RvG_ts_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id_version',"ensembl_gene_id"),filters = "ensembl_transcript_id_version", values = RvG_ts_sig$X, bmHeader = T, )
  RvG_ts_final = merge(RvG_ts_sig, RvG_ts_names, by.x = "X", by.y = "Transcript stable ID version")
  
  ##PopxFib##
  PxF_ts = read.csv("DESeq_CrossbyFibrosis_Full_StringParam_new.csv")
  PxF_ts_sig = PxF_ts[PxF_ts$padj < .1,]  
  PxF_ts_sig = PxF_ts_sig[rowSums(is.na(PxF_ts_sig)) != ncol(PxF_ts_sig), ]
  PxF_ts_sig = PxF_ts_sig[,c(1,3)]
  PxF_ts_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id_version',"ensembl_gene_id"),filters = "ensembl_transcript_id_version", values = PxF_ts_sig$X, bmHeader = T, )
  PxF_ts_final = merge(PxF_ts_sig, PxF_ts_names, by.x = "X", by.y = "Transcript stable ID version")

  ##F2vRBC##
  FvR_ts = read.csv("DESeq_F2vsRBC_Full_StringParam_new.csv")
  FvR_ts_sig = FvR_ts[FvR_ts$padj < .1,]  
  FvR_ts_sig = FvR_ts_sig[rowSums(is.na(FvR_ts_sig)) != ncol(FvR_ts_sig), ]
  FvR_ts_sig = FvR_ts_sig[,c(1,3)]
  FvR_ts_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id_version',"ensembl_gene_id"),filters = "ensembl_transcript_id_version", values = FvR_ts_sig$X, bmHeader = T, )
  FvR_ts_final = merge(FvR_ts_sig, FvR_ts_names, by.x = "X", by.y = "Transcript stable ID version")
  
  ##Fib##
  Fib_ts = read.csv("DESeq_Fibrosis_Full_StringParam_new.csv")
  Fib_ts_sig = Fib_ts[Fib_ts$padj < .1,]  
  Fib_ts_sig = Fib_ts_sig[rowSums(is.na(Fib_ts_sig)) != ncol(Fib_ts_sig), ]
  Fib_ts_sig = Fib_ts_sig[,c(1,3)]
  Fib_ts_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id_version',"ensembl_gene_id"),filters = "ensembl_transcript_id_version", values = Fib_ts_sig$X, bmHeader = T, )
  Fib_ts_final = merge(Fib_ts_sig, Fib_ts_names, by.x = "X", by.y = "Transcript stable ID version")
  
  ##FvGxI##
  FvGxI_ts = read.csv("DESeq_FvGbyInfect_Full_StringParam_new.csv")
  FvGxI_ts_sig = FvGxI_ts[FvGxI_ts$padj < .1,]  
  FvGxI_ts_sig = FvGxI_ts_sig[rowSums(is.na(FvGxI_ts_sig)) != ncol(FvGxI_ts_sig), ]
  FvGxI_ts_sig = FvGxI_ts_sig[,c(1,3)]
  FvGxI_ts_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id_version',"ensembl_gene_id"),filters = "ensembl_transcript_id_version", values = FvGxI_ts_sig$X, bmHeader = T, )
  FvGxI_ts_final = merge(FvGxI_ts_sig, FvGxI_ts_names, by.x = "X", by.y = "Transcript stable ID version")
  
  ##FvRxI##
  FvRxI_ts = read.csv("DESeq_FvRbyInfect_Full_StringParam_new.csv")
  FvRxI_ts_sig = FvRxI_ts[FvRxI_ts$padj < .1,]  
  FvRxI_ts_sig = FvRxI_ts_sig[rowSums(is.na(FvRxI_ts_sig)) != ncol(FvRxI_ts_sig), ]
  FvRxI_ts_sig = FvRxI_ts_sig[,c(1,3)]
  FvRxI_ts_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id_version',"ensembl_gene_id"),filters = "ensembl_transcript_id_version", values = FvRxI_ts_sig$X, bmHeader = T, )
  FvRxI_ts_final = merge(FvRxI_ts_sig, FvRxI_ts_names, by.x = "X", by.y = "Transcript stable ID version")
  
  ##GBCvF2##
  GvF_ts = read.csv("DESeq_GBCvF2_Full_StringParam_new.csv")
  GvF_ts_sig = GvF_ts[GvF_ts$padj < .1,]  
  GvF_ts_sig = GvF_ts_sig[rowSums(is.na(GvF_ts_sig)) != ncol(GvF_ts_sig), ]
  GvF_ts_sig = GvF_ts_sig[,c(1,3)]
  GvF_ts_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id_version',"ensembl_gene_id"),filters = "ensembl_transcript_id_version", values = GvF_ts_sig$X, bmHeader = T, )
  GvF_ts_final = merge(GvF_ts_sig, GvF_ts_names, by.x = "X", by.y = "Transcript stable ID version")
  
  ##GvRxI##
  GvRxI_ts = read.csv("DESeq_GvRbyInfect_Full_StringParam_new.csv")
  GvRxI_ts_sig = GvRxI_ts[GvRxI_ts$padj < .1,]  
  GvRxI_ts_sig = GvRxI_ts_sig[rowSums(is.na(GvRxI_ts_sig)) != ncol(GvRxI_ts_sig), ]
  GvRxI_ts_sig = GvRxI_ts_sig[,c(1,3)]
  GvRxI_ts_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id_version',"ensembl_gene_id"),filters = "ensembl_transcript_id_version", values = GvRxI_ts_sig$X, bmHeader = T, )
  GvRxI_ts_final = merge(GvRxI_ts_sig, GvRxI_ts_names, by.x = "X", by.y = "Transcript stable ID version")
  
  ##Infection##
  Infect_ts = read.csv("DESeq_Infection_Full_StringParam_new.csv")
  Infect_ts_sig = Infect_ts[Infect_ts$padj < .1,]  
  Infect_ts_sig = Infect_ts_sig[rowSums(is.na(Infect_ts_sig)) != ncol(Infect_ts_sig), ]
  Infect_ts_sig = Infect_ts_sig[,c(1,3)]
  Infect_ts_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id_version',"ensembl_gene_id"),filters = "ensembl_transcript_id_version", values = Infect_ts_sig$X, bmHeader = T, )
  Infect_ts_final = merge(Infect_ts_sig, Infect_ts_names, by.x = "X", by.y = "Transcript stable ID version")

  ##merge with marker data##
  RvG_ts_mark=merge(markers,RvG_ts_final, by.x = "FeatureID", by.y = "Gene stable ID")
  write.csv(RvG_ts_mark, "RvG_markers_Fuess.csv", row.names = FALSE)
  
  Fib_ts_mark=merge(markers,Fib_ts_final, by.x = "FeatureID", by.y = "Gene stable ID")
  write.csv(Fib_ts_mark, "Fibrosis_markers_Fuess.csv", row.names = FALSE)
  
  FvGxI_ts_mark = merge(markers,FvGxI_ts_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(FvGxI_ts_mark, "FvGxI_markers_Fuess.csv", row.names = FALSE)
  
  FvR_ts_mark=merge(markers,FvR_ts_final, by.x = "FeatureID", by.y = "Gene stable ID")
  write.csv(FvR_ts_mark, "FvR_markers_Fuess.csv", row.names = FALSE)
  
  #FvRxI_ts_mark=merge(markers,FvRxI_ts_final, by.x = "FeatureID", by.y = "Gene stable ID")##0 markers
  #write.csv(FvRxI_ts_mark, "FvRxI_markers_Fuess.csv", row.names = FALSE)
  
  GvF_ts_mark = merge(markers,GvF_ts_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(GvF_ts_mark, "GvF_markers_Fuess.csv", row.names = FALSE)
  
  GvRxI_ts_mark = merge(markers,GvRxI_ts_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(GvRxI_ts_mark, "GvRxI_markers_Fuess.csv", row.names = FALSE)
  
  Infect_ts_mark = merge(markers,Infect_ts_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(Infect_ts_mark, "Infect_markers_Fuess.csv", row.names = FALSE)
  
  PxF_ts_mark = merge(markers,PxF_ts_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(PxF_ts_mark, "PxF_markers_Fuess.csv", row.names = FALSE)

##last would be the microbiome MS##
##we use the lists of sig genes for each factor of interest here and include their tau value for later analysis##
##like the tagseq MS we need to convert transcript IDs to gene IDs.
  ##beta##
  beta = read.csv("mbio_beta.csv")
  beta_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = beta$Transcript, bmHeader = T, )
  beta_final = merge(beta, beta_names, by.x = "Transcript", by.y = "Transcript stable ID")
  beta_final = beta_final[,c(3,2)]
  beta_mark = merge(markers,beta_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(beta_mark, "beta_markers_mbio.csv", row.names = FALSE)
  
  ##caulo##
  caulo = read.csv("mbio_caulo.csv")
  caulo_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = caulo$Transcript, bmHeader = T, )
  caulo_final = merge(caulo, caulo_names, by.x = "Transcript", by.y = "Transcript stable ID")
  caulo_final = caulo_final[,c(3,2)]
  caulo_mark = merge(markers,caulo_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(caulo_mark, "caulo_markers_mbio.csv", row.names = FALSE)
  
  ##chitin##
  chitin = read.csv("mbio_chitin.csv")
  chitin_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = chitin$Transcript, bmHeader = T, )
  chitin_final = merge(chitin, chitin_names, by.x = "Transcript", by.y = "Transcript stable ID")
  chitin_final = chitin_final[,c(3,2)]
  chitin_mark = merge(markers,chitin_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(chitin_mark, "chitin_markers_mbio.csv", row.names = FALSE)
  
  ##chlam##
  chlam = read.csv("mbio_chlam.csv")
  chlam_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = chlam$Transcript, bmHeader = T, )
  chlam_final = merge(chlam, chlam_names, by.x = "Transcript", by.y = "Transcript stable ID")
  chlam_final = chlam_final[,c(3,2)]
  chlam_mark = merge(markers,chlam_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(chlam_mark, "chlam_markers_mbio.csv", row.names = FALSE)
  
  ##clostrid##
  clostrid = read.csv("mbio_clostrid.csv")
  clostrid_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = clostrid$Transcript, bmHeader = T, )
  clostrid_final = merge(clostrid, clostrid_names, by.x = "Transcript", by.y = "Transcript stable ID")
  clostrid_final = clostrid_final[,c(3,2)]
  clostrid_mark = merge(markers,clostrid_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(clostrid_mark, "clostrid_markers_mbio.csv", row.names = FALSE)
  
  ##div##
  div = read.csv("mbio_diversity.csv")
  div_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id_version',"ensembl_gene_id"),filters = "ensembl_transcript_id_version", values = div$Transcript.ID, bmHeader = T, )
  div_final = merge(div, div_names, by.x = "Transcript.ID", by.y = "Transcript stable ID version")
  div_final = div_final[,c(3,2)]
  div_mark = merge(markers,div_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(div_mark, "div_markers_mbio.csv", row.names = FALSE)
  
  
  ##geo##
  geo = read.csv("mbio_geo.csv")
  geo_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = geo$Transcript, bmHeader = T, )
  geo_final = merge(geo, geo_names, by.x = "Transcript", by.y = "Transcript stable ID")
  geo_final = geo_final[,c(3,2)]
  geo_mark = merge(markers,geo_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(geo_mark, "geo_markers_mbio.csv", row.names = FALSE)

  
  ##gp10##
  gp10 = read.csv("mbio_gp10.csv")
  gp10_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = gp10$Transcript, bmHeader = T, )
  gp10_final = merge(gp10, gp10_names, by.x = "Transcript", by.y = "Transcript stable ID")
  gp10_final = gp10_final[,c(3,2)]
  gp10_mark = merge(markers,gp10_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(gp10_mark, "gp10_markers_mbio.csv", row.names = FALSE)  
  
  ##halo##
  halo = read.csv("mbio_halo.csv")
  halo_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = halo$Transcript, bmHeader = T, )
  halo_final = merge(halo, halo_names, by.x = "Transcript", by.y = "Transcript stable ID")
  halo_final = halo_final[,c(3,2)]
  halo_mark = merge(markers,halo_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(halo_mark, "halo_markers_mbio.csv", row.names = FALSE)  
  
  ##inc##
  inc = read.csv("mbio_inc.csv")
  inc_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = inc$Transcript, bmHeader = T, )
  inc_final = merge(inc, inc_names, by.x = "Transcript", by.y = "Transcript stable ID")
  inc_final = inc_final[,c(3,2)]
  inc_mark = merge(markers,inc_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(inc_mark, "inc_markers_mbio.csv", row.names = FALSE)  
  
  ##methyl##
  methyl = read.csv("mbio_methyl.csv")
  methyl_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = methyl$Transcript, bmHeader = T, )
  methyl_final = merge(methyl, methyl_names, by.x = "Transcript", by.y = "Transcript stable ID")
  methyl_final = methyl_final[,c(3,2)]
  methyl_mark = merge(markers,methyl_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(methyl_mark, "methyl_markers_mbio.csv", row.names = FALSE)  
  
  
  ##nocard##
  nocard = read.csv("mbio_nocard.csv")
  nocard_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = nocard$Transcript, bmHeader = T, )
  nocard_final = merge(nocard, nocard_names, by.x = "Transcript", by.y = "Transcript stable ID")
  nocard_final = nocard_final[,c(3,2)]
  nocard_mark = merge(markers,nocard_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(nocard_mark, "nocard_markers_mbio.csv", row.names = FALSE) 
  
  ##orb##
  orb = read.csv("mbio_orb.csv")
  orb_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = orb$Transcript, bmHeader = T, )
  orb_final = merge(orb, orb_names, by.x = "Transcript", by.y = "Transcript stable ID")
  orb_final = orb_final[,c(3,2)]
  orb_mark = merge(markers,orb_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(orb_mark, "orb_markers_mbio.csv", row.names = FALSE) 
  
  ##pepto##
  pepto = read.csv("mbio_pepto.csv")
  pepto_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = pepto$Transcript, bmHeader = T, )
  pepto_final = merge(pepto, pepto_names, by.x = "Transcript", by.y = "Transcript stable ID")
  pepto_final = pepto_final[,c(3,2)]
  pepto_mark = merge(markers,pepto_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(pepto_mark, "pepto_markers_mbio.csv", row.names = FALSE) 
  
  ##sparto##
  sparto = read.csv("mbio_sparto.csv")
  sparto_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = sparto$Transcript, bmHeader = T, )
  sparto_final = merge(sparto, sparto_names, by.x = "Transcript", by.y = "Transcript stable ID")
  sparto_final = sparto_final[,c(3,2)]
  sparto_mark = merge(markers,sparto_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(sparto_mark, "sparto_markers_mbio.csv", row.names = FALSE) 
  
  ##rubro##
  rubro = read.csv("mbio_rubro.csv")
  rubro_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = rubro$Transcript, bmHeader = T, )
  rubro_final = merge(rubro, rubro_names, by.x = "Transcript", by.y = "Transcript stable ID")
  rubro_final = rubro_final[,c(3,2)]
  rubro_mark = merge(markers,rubro_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  write.csv(rubro_mark, "rubro_markers_mbio.csv", row.names = FALSE) 
  