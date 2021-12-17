##this code will make the heatmap shown in figure 4g of the manuscript##
##the figure displays correlations between abundance of microbial taxa of interest and expression of neutrophil marker genes##
##data is taken from Fuess et al. 2021, mbio and looks at only sig taxa from that paper##

##package data
library(biomaRt)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(pheatmap)

# select mart and data set for biomart (this will help matching gene name to transcript names)
bm <- useMart("ensembl")
bm <- useDataset("gaculeatus_gene_ensembl", mart=bm)

##first load in our marker lists of interest##

Neut = read.csv("Neutrophil_MarkerData.csv")
Neut_markers = Neut[Neut$Neutrophils.P.Value < .1,]
Neut_markers$type = "Neutrophil"
Neut_markers=Neut_markers[,c(1:2,6)]

##then load in our data from the mbio paper and convert transcript ids to gene ids##
##for each taxa we read in the list of sig transcripts, convert to gene ids
##and then we match with neutrophil markers, select just the columns that have gene name and tau
##and finally we rename the column with taus to the taxon name for later merging into a single combined data frame##
  ##beta##
  beta = read.csv("mbio_beta.csv")
  beta_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = beta$Transcript, bmHeader = T, )
  beta_final = merge(beta, beta_names, by.x = "Transcript", by.y = "Transcript stable ID")
  beta_mark = merge(Neut_markers,beta_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  beta_mark = beta_mark[,c(2,5)]
  names(beta_mark)[names(beta_mark) == 'Taus'] <- 'Betaproteobacteria'
  
  ##chlam##
  chlam = read.csv("mbio_chlam.csv")
  chlam_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = chlam$Transcript, bmHeader = T, )
  chlam_final = merge(chlam, chlam_names, by.x = "Transcript", by.y = "Transcript stable ID")
  chlam_mark = merge(Neut_markers,chlam_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  chlam_mark = chlam_mark[,c(2,5)]
  names(chlam_mark)[names(chlam_mark) == 'Taus'] <- 'Chlamydiales'
  
  ##clostrid##
  clostrid = read.csv("mbio_clostrid.csv")
  clostrid_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = clostrid$Transcript, bmHeader = T, )
  clostrid_final = merge(clostrid, clostrid_names, by.x = "Transcript", by.y = "Transcript stable ID")
  clostrid_mark = merge(Neut_markers,clostrid_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  clostrid_mark = clostrid_mark[,c(2,5)]
  names(clostrid_mark)[names(clostrid_mark) == 'Taus'] <- 'Clostridiaceae'
  
  ##geo##
  geo = read.csv("mbio_geo.csv")
  geo_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = geo$Transcript, bmHeader = T, )
  geo_final = merge(geo, geo_names, by.x = "Transcript", by.y = "Transcript stable ID")
  geo_mark = merge(Neut_markers,geo_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  geo_mark = geo_mark[,c(2,5)]
  names(geo_mark)[names(geo_mark) == 'Taus'] <- 'Geodermatophilaceae'
  
  ##halo##
  halo = read.csv("mbio_halo.csv")
  halo_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = halo$Transcript, bmHeader = T, )
  halo_final = merge(halo, halo_names, by.x = "Transcript", by.y = "Transcript stable ID")
  halo_mark = merge(Neut_markers,halo_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  halo_mark = halo_mark[,c(2,5)]
  names(halo_mark)[names(halo_mark) == 'Taus'] <- 'Halomonadaceae'
  
  ##inc##
  inc = read.csv("mbio_inc.csv")
  inc_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = inc$Transcript, bmHeader = T, )
  inc_final = merge(inc, inc_names, by.x = "Transcript", by.y = "Transcript stable ID")
  inc_mark = merge(Neut_markers,inc_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  inc_mark = inc_mark[,c(2,5)]
  names(inc_mark)[names(inc_mark) == 'Taus'] <- 'Incertae_Sedis_XI'
  
  ##nocard##
  nocard = read.csv("mbio_nocard.csv")
  nocard_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = nocard$Transcript, bmHeader = T, )
  nocard_final = merge(nocard, nocard_names, by.x = "Transcript", by.y = "Transcript stable ID")
  nocard_mark = merge(Neut_markers,nocard_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  nocard_mark = nocard_mark[,c(2,5)]
  names(nocard_mark)[names(nocard_mark) == 'Taus'] <- 'Nocardiaceae'
  
  ##orb##
  orb = read.csv("mbio_orb.csv")
  orb_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = orb$Transcript, bmHeader = T, )
  orb_final = merge(orb, orb_names, by.x = "Transcript", by.y = "Transcript stable ID")
  orb_mark = merge(Neut_markers,orb_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  orb_mark = orb_mark[,c(2,5)]
  names(orb_mark)[names(orb_mark) == 'Taus'] <- 'Orbaceae'
  
  ##pepto##
  pepto = read.csv("mbio_pepto.csv")
  pepto_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = pepto$Transcript, bmHeader = T, )
  pepto_final = merge(pepto, pepto_names, by.x = "Transcript", by.y = "Transcript stable ID")
  pepto_mark = merge(Neut_markers,pepto_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  pepto_mark = pepto_mark[,c(2,5)]
  names(pepto_mark)[names(pepto_mark) == 'Taus'] <- 'Peptostreptococcaceae'
  
  ##rubro##
  rubro = read.csv("mbio_rubro.csv")
  rubro_names <- getBM(mart=bm, attributes=c('ensembl_transcript_id',"ensembl_gene_id"),filters = "ensembl_transcript_id", values = rubro$Transcript, bmHeader = T, )
  rubro_final = merge(rubro, rubro_names, by.x = "Transcript", by.y = "Transcript stable ID")
  rubro_mark = merge(Neut_markers,rubro_final, by.x = "FeatureID", by.y = "Gene stable ID")  
  rubro_mark = rubro_mark[,c(2,5)]
  names(rubro_mark)[names(rubro_mark) == 'Taus'] <- 'Rubrobacteraceae'

##so then we need to merge all these taxon specific data frames into one data frame which we will build the heatmap from##
merge1=merge(beta_mark, chlam_mark, by = "FeatureName", all.x = TRUE, all.y = TRUE)
merge2=merge(merge1, clostrid_mark, by = "FeatureName", all.x = TRUE, all.y = TRUE)
merge3=merge(merge2, geo_mark, by = "FeatureName", all.x = TRUE, all.y = TRUE)
merge4=merge(merge3, halo_mark, by = "FeatureName", all.x = TRUE, all.y = TRUE)
merge5=merge(merge4, inc_mark, by = "FeatureName", all.x = TRUE, all.y = TRUE)
merge6=merge(merge5, nocard_mark, by = "FeatureName", all.x = TRUE, all.y = TRUE)
merge7=merge(merge6, orb_mark, by = "FeatureName", all.x = TRUE, all.y = TRUE)
merge8=merge(merge7, pepto_mark, by = "FeatureName", all.x = TRUE, all.y = TRUE)
merge9=merge(merge8, rubro_mark, by = "FeatureName", all.x = TRUE, all.y = TRUE)
final=merge9
final2 <- final[,-1]
rownames(final2) <- final[,1]

##make the heatmap##
pheatmap(final2, color = viridis(500), 
         show_colnames = TRUE, scale = "none",
         annotation_legend = TRUE)






