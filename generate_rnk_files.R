library(tidyverse)
in_directory = "~/Documents/R/ClockCorr2/diff_expr_tables2/"
setwd(in_directory)
large_tissue_fscores = list.files(pattern = "*_edgeR.csv")

out_directory = "~/Documents/R/ClockCorr2/clockcorr3/literatureCompare/"
get_tissue_name = function(str){
  # tissue_name = str_extract(str, pattern = "(?<=\\()(.*)(?=\\))")
  tissue_name = str_extract(str, pattern = "(.*)(?=_dthhrdy)")
  return(tissue_name)
}
rank_fscores = function(tissue, write = T){     #generates GSEA file for lm analysis files
  large_tissue = read_csv(tissue)
  large_tissue = drop_na(large_tissue)
  #large_tissue$rnk = -log(large_tissue$FDR, 2)   #FIXME, uncomment for log(FDR)
  large_tissue$rnk = abs(large_tissue$logFC)
  large_tissue = arrange(large_tissue, desc(rnk)  )
  if(write){
    setwd(out_directory)
    write.table(select(large_tissue, `...1`, rnk), file = paste0(get_tissue_name(tissue), "LFC_preranked.rnk"), sep = '\t', quote = F , row.names = F, col.names = T)
    setwd(in_directory)
  }
  return(large_tissue)
}

tissue_lm_dfs = sapply(large_tissue_fscores, rank_fscores, write = T)
