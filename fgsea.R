library(fgsea)
library(data.table)
library(ggplot2)
library(stringr)


setwd("~/Documents/R/ClockCorr2/rnk_files/")

pathways <- gmtPathways("../circ_out_exemplar.gmt")
setwd("~/Documents/R/ClockCorr2/clockcorr3/literatureCompare/")

rnk_files = list.files( pattern = "*LFC_preranked.rnk")

run_fgsea = function(file, gsea_param = 1){
  ranks <- read.table(file, header=T, colClasses = c("character", "numeric"))
  ranks <- setNames(ranks$rnk, ranks$...1)
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    minSize  = 1,
                    eps = 0.0,
                    gseaParam = gsea_param,
                    maxSize  = 500)
}

fgseaRes = lapply(rnk_files, run_fgsea, gsea_param=1)
names(fgseaRes) = sapply(rnk_files, str_remove, '.rnk')
pvals = sapply(fgseaRes, `[[`, 2)
pvals = t(pvals)
colnames(pvals) = c("circ_output", "core_clock")
#pvals = apply(pvals, 2, p.adjust, method = 'BH')
pvals = round(pvals, 3)
write.table(pvals, "../clock_circ_output_fGSEA_abs_LFC_ranked.csv", sep=',', quote = F, col.names = NA)

plot_pathway = function(file, pathway, gsea_param = 1){
  name = names(pathway)
  ranks <- read.table(file, header=T, colClasses = c("character", "numeric"))
  ranks <- setNames(ranks$rnk, ranks$...1)
  plotEnrichment(unlist(pathway), gseaParam = gsea_param, ranks) + labs(title=paste(name, str_remove(file, '.rnk')))

}
plot_pathway(rnk_files[1], pathways[1], gsea_param = 1)
