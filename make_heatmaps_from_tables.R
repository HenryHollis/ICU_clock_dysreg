library(gplots)
in_dir = "~/Documents/R/ClockCorr2/diff_expr_tables/LFC_FDR_tables/"
setwd(in_dir)

#making ubiquitous genes (FDR sig in 25 of 25 tissues) LFC plot:
ubiq_LFC = read.csv("ubiq_FDR_cutoff_25_of_25_tissues_LFC_vals.csv", row.names = 1)
ubiq_LFC_mat = sapply(ubiq_LFC, as.numeric)
rownames(ubiq_LFC_mat) = rownames(ubiq_LFC)
breaks = seq(-2, 2, by = .1)
heatmap.2(ubiq_LFC_mat, col = topo.colors(40), breaks = breaks,trace = 'none', tracecol = 'black', key.title = "Log Fold Change", margins = c(15, 10) )

#making ubiquitous genes (FDR sig and LFC high in 22 of 25 tissues) LFC plot:
ubiq_LFC = read.csv("ubiq_FDR_LFC_22_of_25_tissues_LFC_values.csv", row.names = 1)
ubiq_LFC_mat = sapply(ubiq_LFC, as.numeric)
rownames(ubiq_LFC_mat) = rownames(ubiq_LFC)
breaks = seq(-3, 3, by = .1)
heatmap.2(ubiq_LFC_mat, col = topo.colors(length(breaks)-1), breaks = breaks,trace = 'none', tracecol = 'black', key.title = "Log Fold Change", margins = c(15, 10) )

#making circ_output FDR and LFC plot:
circadian_output_FDR_LFC = as.data.frame(read.csv("circadian_output_FDR_LFC_cutoff_LFC_values.csv", row.names = 1))
circadian_output_FDR_LFC[circadian_output_FDR_LFC==FALSE] = NA
circadian_output_FDR_LFC_mat = sapply(circadian_output_FDR_LFC, as.numeric)
rownames(circadian_output_FDR_LFC_mat) = rownames(circadian_output_FDR_LFC)
heatmap.2(circadian_output_FDR_LFC_mat, Rowv = F, Colv = F, col = topo.colors, trace = 'none', tracecol = 'black', key.title = "Log Fold Change", margins = c(15, 5))

#making circ_output FDR and LFC plot:
circadian_output_FDR = as.data.frame(read.csv("circadian_output_FDR_cutoff_LFC_values.csv", row.names = 1))
circadian_output_FDR[circadian_output_FDR==FALSE] = NA
circadian_output_FDR_mat = sapply(circadian_output_FDR, as.numeric)
rownames(circadian_output_FDR_mat) = rownames(circadian_output_FDR)
heatmap.2(circadian_output_FDR_mat, Rowv = F, Colv = F, col = topo.colors, trace = 'none', tracecol = 'black', key.title = "Log Fold Change", margins = c(15, 5))



#Making exemplar LFC_and_FDR plot:
exemplar_LFC = as.data.frame(read.csv("exemplar_LFC_vals_LFC_and_FDR_cutoff.csv", row.names = 1))
exemplar_LFC[exemplar_LFC==FALSE]=NA
exemplar_LFC_mat = sapply(exemplar_LFC, as.numeric)
rownames(exemplar_LFC_mat) = rownames(exemplar_LFC)
heatmap.2(exemplar_LFC_mat, col = topo.colors(40),trace = 'none', tracecol = 'black', key.title = "Log Fold Change", margins = c(15, 10), Colv = F, Rowv = F )

#making exemplar FDR plot:
exemplar_FDR = as.data.frame(read.csv("exemplar_LFC_vals_FDR_cutoff.csv", row.names = 1))
exemplar_FDR[exemplar_FDR==FALSE] = NA
exemplar_FDR_mat = sapply(exemplar_FDR, as.numeric)

rownames(exemplar_FDR_mat) = rownames(exemplar_FDR)
heatmap.2(exemplar_FDR_mat, col = topo.colors(40),trace = 'none', tracecol = 'black', key.title = "Log Fold Change", margins = c(15, 10))
