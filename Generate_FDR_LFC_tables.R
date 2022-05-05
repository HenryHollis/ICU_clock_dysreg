#This file takes in differential expression tables from edgeR and creates relevant
#tables and lists used in our analysis, including those used make the FDR<0.05 list 
#and the FDR<0.5 & abs(L2FC) >0.58

library(tidyverse)
library(doParallel)

BHthreshold = 0.05
LFCcutoff = F      #change to T to enforce log fold change of 0.58
num_tissues = 25   #the number of tissues required to be present in for ubiq lists

in_dir = "~/Documents/R/ClockCorr2/diff_expr_tables/"
out = "~/Documents/R/ClockCorr2/diff_expr_tables/LFC_FDR_tables/"
setwd(in_dir)

registerDoParallel(cores = 4)
edger.files = list.files(pattern = "*_edgeR.csv")
gene_list = c("PER3" , "CIART" ,"NPAS2", "PER2",  "NR1D2", "CLOCK" ,"ARNTL", "CRY2" , "CRY1",  "PER1" , "NR1D1" ,"DBP","TEF")
exemplar_list = c("RORA", "RORB", "RORC", "NFIL3", "FBXL3", "CSNK1D", "CSNK1E", "PER3" , "CIART" ,"NPAS2", "PER2",  "NR1D2", "CLOCK" ,"ARNTL", "CRY2" , "CRY1",  "PER1" , "NR1D1" ,"DBP","TEF")
circadian_output = c("TSC22D3",	"USP2"	,"TSPAN4"	,"HLF",	"FMO2",	"GRAMD4",	"DTX4",	"CDKN1A",	"MMP14",	"POR",	"WEE1",	"SLC16A1",	"LEO1",	"STK35",	"LONRF3",	"COQ10B",	"BHLHE41",	"TMEM57",	"TNS2",	"HSP90AA1",	"THRA",	"HSPB1",	"CALR",	"GLUL",	"CLDN5",	"FKBP5",	"STIP1",	"HERPUD1",	"LITAF",	"PFKFB3",	"FMO1",	"MTHFD1L",	"VPS13A",	"SRM",	"AHCTF1",	"CMTM6",	"TMEM86A"	,"SLC46A3",	"EDN1",	"FBN1",	"EPHX1",	"TIMP3",	"HDAC4",	"FABP7",	"DAPK1",	"MAP2K7",	"PDCD4",	"PRODH",	"MCAM",	"HSPH1",	"CNTFR",	"REV1",	"CYS1",	"SLC1A5",	"KLF9",	"DGAT2",	"KLF13",	"ACSL1",	"CIRBP",	"SNRK",	"HSPA4L",	"ST6GALNAC3",	"RHOJ",	"AK3",	"ANGPTL2",	"TMEM33",	"RNF145",	"MARVELD1",	"LONRF1",	"FUS",	"MT2A",	"MAPRE2",	"LY6H",	"PRR13",	"CYRIA",	"NAA60")

get_tissue_name = function(str){
  tissue_name = str_extract(str, pattern = ".*(?=_dthhrdy)")
  return(tissue_name)
}

genes_for_each_tiss = foreach(i=1:length(edger.files), .combine = c) %dopar% {
  name = get_tissue_name(edger.files[i])
  top.table = (read.csv(edger.files[i], row.names=1))
  sig.table = top.table[which(as.numeric(top.table$FDR) < BHthreshold) ,]        #limit analysis to significant genes
  if (LFCcutoff){
      sig.table = top.table[which((as.numeric(top.table$FDR) < BHthreshold) & abs(as.numeric(top.table$logFC)) > 0.5849625) ,]        #limit analysis to significant genes
  }
  gene_names = rownames(sig.table)
  gene_FC = sig.table$logFC
  gene_FDR = sig.table$FDR
  gene_PVal = sig.table$PValue
  tissue_df = rbind(gene_names, gene_FC, gene_FDR, gene_PVal)
  tissue_df = tissue_df[,order(tissue_df[1,])]
  return(list(tissue_df))
}


genes_universe = unique(purrr::reduce(lapply(genes_for_each_tiss, function(x) x[1,]), c))

logical_table = foreach(i = 1:length(genes_for_each_tiss), .combine = rbind) %do% {     #returns a boolean table numTissues x length(geneuniverse) with true is found, false not found
  tissue_has_gene = genes_universe %in% genes_for_each_tiss[[i]][1,]  #logical vector of length equal to genes_universe
}

gene_counts_across_tiss = apply(logical_table, 2, sum)

####plotting histogram of gene_counts_across_tissues ####
ggplot() + 
  aes(gene_counts_across_tiss)+ 
  geom_histogram(binwidth=1, colour="black", fill="white")+
  xlab("Number of Tissues")+
  ylab("Number of Genes")


ubiq_list = genes_universe[which(gene_counts_across_tiss >=num_tissues)]


generate_summary_tables = function(param = 2, genes){
  ### Params: 2 is Log fold change
  ###         3 is FDR
   genes = sort(genes)
   tab = foreach(i = 1:length(genes_for_each_tiss), .combine = rbind) %dopar% {
      tissue_genes_in_list = genes_for_each_tiss[[i]][1,] %in% genes  #logical vector  indicating which genes in tissue are in "genes" list
      gene_row =  genes %in% genes_for_each_tiss[[i]][1,] #logical vect indicating which genes 
      gene_row[gene_row] = genes_for_each_tiss[[i]][param, tissue_genes_in_list]
      return(gene_row)
   }
   tab = t(tab)
   tab = as.data.frame(tab)
   colnames(tab) = lapply(edger.files, get_tissue_name)
   rownames(tab) = sort(genes)
   return(tab)
}

#with the above function we can create tables of LFC, PVal or FDR for any arbitrary list:
setwd(out)
ubiq_FDR_LFC_22_of_25_tissues_LFC_values = generate_summary_tables(2, ubiq_list)
#write.table(ubiq_FDR_LFC_22_of_25_tissues_LFC_values, "ubiq_FDR_LFC_22_of_25_tissues_LFC_values.csv", sep=',', quote=F, col.names = NA)
ubiq_FDR_cutoff_25_of_25_tissues_LFC_vals = generate_summary_tables(2, ubiq_list)
#write.table(ubiq_FDR_cutoff_25_of_25_tissues_LFC_vals, "ubiq_FDR_cutoff_25_of_25_tissues_LFC_vals.csv", sep=',', quote=F, col.names = NA)
ubiq_FDR = generate_summary_tables(3, ubiq_list)
#write.table(ubiq_FDR, "ubiq_FDR.csv", sep=',', quote=F, col.names = NA)
matrix_clock_LFC = generate_summary_tables(2, gene_list)
#write.table(matrix_clock_LFC, "matrix_clock_LFC.csv", sep=',', quote=F, col.names = NA)
matrix_clock_FDR = generate_summary_tables(2, gene_list)
#write.table(matrix_clock_FDR, "matrix_clock_FDR.csv", sep=',', quote=F, col.names = NA)

exemplar_LFC_vals_LFC_and_FDR_cutoff = generate_summary_tables(2, exemplar_list)
#write.table(exemplar_LFC_vals_LFC_and_FDR_cutoff, "exemplar_LFC_vals_LFC_and_FDR_cutoff.csv", sep=',', quote=F, col.names = NA)
exemplar_LFC_vals_FDR_cutoff = generate_summary_tables(2, exemplar_list)
#write.table(exemplar_LFC_vals_FDR_cutoff, "exemplar_LFC_vals_FDR_cutoff.csv", sep=',', quote=F, col.names = NA)

circadian_output_FDR = generate_summary_tables(3, circadian_output)
#write.table(circadian_output_FDR, "circadian_output_FDR.csv", sep=',', quote=F, col.names = NA)

circadian_output_FDR_LFC_cutoff_LFC_values = generate_summary_tables(2, circadian_output)
#write.table(circadian_output_FDR_LFC_cutoff_LFC_values, "circadian_output_FDR_LFC_cutoff_LFC_values.csv", sep=',', quote=F, col.names = NA)

circadian_output_FDR_cutoff_LFC_values = generate_summary_tables(2, circadian_output)
#write.table(circadian_output_FDR_cutoff_LFC_values, "circadian_output_FDR_cutoff_LFC_values.csv", sep=',', quote=F, col.names = NA)


all_sig_genes_FDR = generate_summary_tables(3, genes_universe)
#write.table(all_sig_genes_FDR, "all_sig_genes_FDR.csv", sep=',', quote=F, col.names = NA)

all_sig_genes_FDR_cutoff_LFC_values = generate_summary_tables(2, genes_universe)
#write.table(all_sig_genes_FDR_cutoff_LFC_values, "all_sig_genes_FDR_cutoff_LFC_values.csv", sep=',', quote=F, col.names = NA)

#making a histogram of significant circadian genes in each tissue:
# circ_frac = c()
# core_frac = c()
# for(i in 1:length(edger.files)){
#   circ_frac = c(circ_frac, sum(circadian_output %in% genes_for_each_tiss[[i]][1,])/length(circadian_output))
#   core_frac = c(core_frac, sum(exemplar_list %in% genes_for_each_tiss[[i]][1,])/length(exemplar_list))
#   
# }
# Tissues = sapply(edger.files, get_tissue_name)
# library(reshape)
# df1 <- data.frame(circ_frac, core_frac, Tissues)
# colnames(df1) = c('Circadian Output', 'Core Clock', 'Tissues')
# df2 <- melt(df1, id.vars='Tissues')
# 
# ggplot(df2, aes(x=Tissues, y=value, fill=variable)) +
#   geom_bar(stat='identity', position='dodge')+
#   scale_fill_manual(values=c("#B2E2F0", "#FFB346"))+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ylab("Fraction of Genes FDR<0.05")
