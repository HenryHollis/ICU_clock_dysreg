# This script does a manual  overrepresentation approach to enriching, similar to EnrichR or DAVID ncbi.
# We are testing for overexpression of the 3 custom clock gene sets which is the motivation for doing this manually...

setwd("~/Documents/R/ClockCorr2/diff_expr_tables2/")
core_clock = c("RORA","RORB","RORC","NFIL3","FBXL3","CSNK1D","CSNK1E","PER3","CIART","NPAS2","PER2","NR1D2","CLOCK","ARNTL","CRY2","CRY1","PER1","NR1D1","DBP","TEF")
#BTCH_8_inclusive = c("PER2","PER1","ARNTL","PER3","NR1D1","DBP","TSC22D3","USP2","NR1D2","TSPAN4","NPAS2","TEF","HLF","FMO2","GRAMD4","DTX4","CDKN1A","MMP14","CLOCK","POR","WEE1","SLC16A1","NFIL3","LEO1","STK35","LONRF3","COQ10B","BHLHE41","TMEM57","CIART","TNS2","HSP90AA1","THRA","HSPB1","CALR","CRY1","GLUL","RORC","CLDN5","FKBP5","STIP1","HERPUD1","LITAF","PFKFB3","FMO1","MTHFD1L","VPS13A","SRM","AHCTF1","CMTM6","TMEM86A","SLC46A3","EDN1","FBN1","EPHX1","TIMP3","HDAC4","FABP7","CRY2","DAPK1","MAP2K7","PDCD4","PRODH","MCAM","HSPH1","CNTFR","REV1","CYS1","SLC1A5","KLF9","DGAT2","KLF13","ACSL1","CIRBP","SNRK","HSPA4L","ST6GALNAC3","RHOJ","AK3","ANGPTL2","TMEM33","RNF145","MARVELD1","LONRF1","FUS","MT2A","MAPRE2","LY6H","PRR13","FAM49A","NAA60")
circadian_output = c("TSC22D3","USP2","TSPAN4","HLF","FMO2","GRAMD4","DTX4","CDKN1A","MMP14","POR","WEE1","SLC16A1","LEO1","STK35","LONRF3","COQ10B","BHLHE41","TMEM57","TNS2","HSP90AA1","THRA","HSPB1","CALR","GLUL","CLDN5","FKBP5","STIP1","HERPUD1","LITAF","PFKFB3","FMO1","MTHFD1L","VPS13A","SRM","AHCTF1","CMTM6","TMEM86A","SLC46A3","EDN1","FBN1","EPHX1","TIMP3","HDAC4","FABP7","DAPK1","MAP2K7","PDCD4","PRODH","MCAM","HSPH1","CNTFR","REV1","CYS1","SLC1A5","KLF9","DGAT2","KLF13","ACSL1","CIRBP","SNRK","HSPA4L","ST6GALNAC3","RHOJ","AK3","ANGPTL2","TMEM33","RNF145","MARVELD1","LONRF1","FUS","MT2A","MAPRE2","LY6H","PRR13","FAM49A","NAA60")
tissue_paths = list.files(pattern = "*_edgeR.csv")

Fisher_exact_on_tiss = function(tissue_file, FDR_cutoff = 0.05, custom_list, LFC = F){
  tissue = read.csv(tissue_file)
  sig_genes = tissue[which(tissue$FDR < FDR_cutoff) ,1] #find the upreg V/A genes
  if(LFC){
    sig_genes = tissue[which((tissue$FDR < FDR_cutoff) & (tissue$logFC > 0.5849625)) ,1] #find the upreg V/A genes
  }
  res = do_fisher_test(custom_list, sig_genes, tissue$X)
  return(res$p.value)
}

do_fisher_test = function(custom_set_genes, VA_sig_genes, background){
  non_sig_lit = setdiff(background, custom_set_genes) #non-significant lit genes
  non_sig_VA = setdiff(background, VA_sig_genes) #non-significant lit genes
  top_left = length(intersect(VA_sig_genes, custom_set_genes)) #number of (sig V/A and sig lit genes)
  bottom_left = length(intersect(non_sig_lit, VA_sig_genes)) # number of (sig V/A genes, nonsig lit genes)
  #top_right = length(custom_set_genes) - top_left                        #number of nonsig VA, sig list genes    
  top_right = length(intersect(non_sig_VA, custom_set_genes))
  #bottom_right = length(non_sig_lit) - bottom_left              #number of nonsig V/A nonsig lit genes
  bottom_right = length(intersect(non_sig_lit, non_sig_VA))
  top = cbind(top_left, top_right)
  tmp = cbind(bottom_left, bottom_right)
  table = rbind(top, tmp)
  #View(table)
  test = fisher.test(table)
  return(test)
}

get_tissue_name = function(str){
  # tissue_name = str_extract(str, pattern = "(?<=\\()(.*)(?=\\))")
  tissue_name = str_extract(str, pattern = "(.*)(?=_dthhrdy)")
  return(tissue_name)
}
core_clock_FDR_only_BHQ = round(p.adjust(sapply(tissue_paths, Fisher_exact_on_tiss, custom_list = core_clock), method = 'BH'), 3)
core_clock_FDR_LFC_BHQ = round(p.adjust(sapply(tissue_paths, Fisher_exact_on_tiss, custom_list = core_clock, LFC = T), method = 'BH'), 3)

circ_output_FDR_only_BHQ = round(p.adjust(sapply(tissue_paths, Fisher_exact_on_tiss, custom_list = circadian_output), method = 'BH'), 3)
circ_output_FDR_LFC_BHQ = round(p.adjust(sapply(tissue_paths, Fisher_exact_on_tiss, custom_list = circadian_output, LFC = T), method = 'BH'), 3)
tissue_names = lapply(tissue_paths, get_tissue_name)

out = cbind(tissue_names, circ_output_FDR_LFC_BHQ, core_clock_FDR_LFC_BHQ)
write.table(out, "../clockcorr3/clock_circ_output_enriching2.csv", sep=',', quote = F, row.names = F)
