#!/usr/bin/env Rscript
library(enrichR)
library(doParallel)
library(tidyverse)
library(gplots)
registerDoParallel(cores = 4)

in_dir = "~/Documents/R/ClockCorr2/diff_expr_tables/LFC_FDR_tables/"
setwd(in_dir)

exemplar_list = c("RORA", "RORB", "RORC", "NFIL3", "FBXL3", "CSNK1D", "CSNK1E", "PER3" , "CIART" ,"NPAS2", "PER2",  "NR1D2", "CLOCK" ,"ARNTL", "CRY2" , "CRY1",  "PER1" , "NR1D1" ,"DBP","TEF")
gene_list = c("PER3" , "CIART" ,"NPAS2", "PER2",  "NR1D2", "CLOCK" ,"ARNTL", "CRY2" , "CRY1",  "PER1" , "NR1D1" ,"DBP","TEF")
circadian_output = c("TSC22D3",	"USP2"	,"TSPAN4"	,"HLF",	"FMO2",	"GRAMD4",	"DTX4",	"CDKN1A",	"MMP14",	"POR",	"WEE1",	"SLC16A1",	"LEO1",	"STK35",	"LONRF3",	"COQ10B",	"BHLHE41",	"TMEM57",	"TNS2",	"HSP90AA1",	"THRA",	"HSPB1",	"CALR",	"GLUL",	"CLDN5",	"FKBP5",	"STIP1",	"HERPUD1",	"LITAF",	"PFKFB3",	"FMO1",	"MTHFD1L",	"VPS13A",	"SRM",	"AHCTF1",	"CMTM6",	"TMEM86A"	,"SLC46A3",	"EDN1",	"FBN1",	"EPHX1",	"TIMP3",	"HDAC4",	"FABP7",	"DAPK1",	"MAP2K7",	"PDCD4",	"PRODH",	"MCAM",	"HSPH1",	"CNTFR",	"REV1",	"CYS1",	"SLC1A5",	"KLF9",	"DGAT2",	"KLF13",	"ACSL1",	"CIRBP",	"SNRK",	"HSPA4L",	"ST6GALNAC3",	"RHOJ",	"AK3",	"ANGPTL2",	"TMEM33",	"RNF145",	"MARVELD1",	"LONRF1",	"FUS",	"MT2A",	"MAPRE2",	"LY6H",	"PRR13",	"CYRIA",	"NAA60")

#enriching the circadian output genes
#GOBP_circOut_enriched = enrichr(circadian_output,c("GO_Biological_Process_2021", "Reactome_2016", "KEGG_2021_Human") )
#all_results = rbind(GOBP_circOut_enriched[[1]], GOBP_circOut_enriched[[2]], GOBP_circOut_enriched[[3]])
#all_results$P.value = all_results$Adjusted.P.value    # quirk of plotEnrich function is that it won't plot BH.Q vals
#plotEnrich(all_results, showTerms = 35, numChar = 400, y = "Count", title = "Over Representation Analysis of Ubiquitous Genes")

#enriching the ubiquitous list:
ubiq_FDR = read.csv("ubiq_FDR_LFC_cutoffs_10_of_25_tissues_genes.csv", row.names = 1)
#ubiq_genes = rownames(ubiq_FDR)
ubiq_genes = ubiq_FDR[-1, 1]

#FOR KEGG: KEGG_2021_Human
#For GOBP: GO_Biological_Process_2021

GOBP_ubiq_enriched = enrichr(ubiq_genes, "KEGG_2021_Human")
#all_results = rbind(GOBP_ubiq_enriched[[1]], GOBP_ubiq_enriched[[2]], GOBP_ubiq_enriched[[3]])  #uncomment for enriching multiple libs
all_results = rbind(GOBP_ubiq_enriched[[1]])
display_cutoff = length(which(all_results$Adjusted.P.value < 0.01))
all_results$P.value = all_results$Adjusted.P.value    # quirk of plotEnrich function is that it won't plot BH.Q vals
plotEnrich(all_results,showTerms = 20 ,numChar = 70, y = "Ratio", title = "Over Representation Analysis of genes FDR<0.05, in 20/25 tissues")


get_all_sig_pathways = function(FDR_tissue_table, threshold = 0.05, pathway = "GO_Biological_Process_2021"){
  res = foreach(i = 1:ncol(FDR_tissue_table), .combine = c) %do% {
     #top_ten_percent = floor(length(which(FDR_tissue_table[,i]!=F))/10) #grab the top 10% most significantly expressed genes
     #ordered_table = FDR_tissue_table[order(FDR_tissue_table[,i]),]
     #sig_genes = rownames(FDR_tissue_table)[which(FDR_tissue_table[,i] != F & FDR_tissue_table[,i] < 0.005)]
     sig_genes = rownames(FDR_tissue_table)[which(FDR_tissue_table[,i] != F)]
     
     enriched_pathways = enrichr(sig_genes, pathway )
     local_pathways = enriched_pathways[[1]][["Adjusted.P.value"]]
     names(local_pathways) = enriched_pathways[[1]][["Term"]]
     return(list(local_pathways))
  }
  return(res)
  
}

#plotting heatmap of significant pathways among all tissues:
#procedures is as follows:
#enrich the list of significant (FDR<0.05) and LFC > 0.5849625 genes in the GOBP pathway in enrichr
#collect all the pathways (which are sig at FDR <0.1) and rank them by how many tissues they're enriched in
#plot a heatmap of the P-value of pathway i for tissue j for all i, j

#all_sig_genes_FDR = read.csv("all_sig_genes_FDR_cutoff_LFC_values.csv", row.names = 1)
all_sig_genes_FDR = read.csv("all_sig_genes_LFC_cutoff_FDR.csv", row.names = 1)

total_pathways = get_all_sig_pathways(all_sig_genes_FDR, pathway = "GO_Biological_Process_2021")

total_enriched_pathways = purrr::reduce(lapply(total_pathways, function(x) x[]), c)
total_enriched_pathways = total_enriched_pathways[total_enriched_pathways<=0.1]
tissue_pathway_freq = table(names(total_enriched_pathways))
tissue_pathway_freq = tissue_pathway_freq[order(tissue_pathway_freq, decreasing = T)]
pathways_of_interest = names(tissue_pathway_freq)[1:40]


heatmap = foreach(i = 1:25, .combine = cbind) %do% {
  log_vec = pathways_of_interest %in% names(total_pathways[[i]]) # logical vector corresponding to pathways_of_interest
  found_in_tissue = pathways_of_interest[log_vec]                   # names of pathways enriched in tissue i
  log_vec[log_vec] = unname(total_pathways[[i]][found_in_tissue])
  
  
  log_vec[which(log_vec==F)] = NA
  return(log_vec)
}
colnames(heatmap) = colnames(all_sig_genes_FDR)
rownames(heatmap) = pathways_of_interest
breaks <- c(0, 0.05, 0.1, 0.15, .2, 0.25)
my_pallet = c("#9900ff", topo.colors(3), "#ffffff")
heatmap.2(as.matrix(heatmap), trace = 'none', col = my_pallet, breaks = breaks, na.color = 'white', key.title ='BH.q value', tracecol = 'black' ,Colv = F, Rowv = F, margins = c(15, 32))


