#This script takes genes expresseion matrices and runs a modified version of the CalcDeltaCCD function from the
#deltaccd library. In the modified function, data gathered from different facilities is permuted seperately
#in order to account for batch effects.

library(tidyverse)
library(doParallel)
library(doRNG)
library(deltaccd)
library(edgeR)
source("~/Documents/R/ClockCorrelation/Rscripts/deltaccd_utils.R")
source("~/Documents/R/ClockCorr2/Corr2scripts/custom_CalcDeltaCCD.R")
registerDoParallel(cores = 4)
in_directory = "~/Documents/R/ClockCorr2/GTEx_counts/"
setwd(in_directory)
gene_list = c("PER3" , "CIART" ,"NPAS2", "PER2",  "NR1D2", "CLOCK" ,"ARNTL", "CRY2" , "CRY1",  "PER1" , "NR1D1" ,"DBP","TEF")

gtex_files=list.files(pattern = "*.csv")
gtex_files = sort(gtex_files)
setwd(in_directory)
count_files = list.files(pattern = "*.csv")
count_files = sort(count_files)
try(if(!assertthat::not_empty(gtex_files)) stop(paste("No .csv files found in", getwd())))

refCor = getRefCor('human', 'pan', FALSE)

get_tissue_name = function(str){
  tissue_name = str_extract(str, pattern = "(?<=\\()(.*)(?=\\))")
  return(tissue_name)
}
collate_data = function(file_name, small = T){
  tissue_data = read.csv(file_name)                           # read csv
  rm_rows = which(!(tissue_data$Description %in% gene_list))   # get row indices of throwaway genes
  small_tissue_data = tissue_data
  if(small){
    small_tissue_data = tissue_data[-rm_rows,] 
    small_tissue_data = bind_rows(tissue_data[1:9, ], small_tissue_data)
  }
  small_tissue_data =
    select(small_tissue_data, -Name) %>%
    unite("Info", Description,Info, sep = "", remove = T) 
  rm(tissue_data)
  
  small_tissue_data = t(small_tissue_data)
  
  if(dim(small_tissue_data)[1] <= 1){return(NA)}
  
  colnames(small_tissue_data) = small_tissue_data[1, ]
  small_tissue_data = as_tibble(small_tissue_data[-1, ])
  
  # small_tissue_data[9,which(small_tissue_data[9,] == "0.0")] = "V"
  # small_tissue_data[9,c(which(small_tissue_data[9,] == "1.0"), which(small_tissue_data[9,] == "2.0"))] = "A"
  # rm_cols = c(which(small_tissue_data[9,] != "3.0"), which(small_tissue_data[9,] != "4.0"))
  # small_tissue_data = small_tissue_data[, -rm_cols]
  
  small_tissue_data = mutate(small_tissue_data, DTHHRDY = replace(DTHHRDY, (DTHHRDY == "0.0" ), "V"))
  small_tissue_data = mutate(small_tissue_data, DTHHRDY = replace(DTHHRDY, (DTHHRDY == "1.0"|DTHHRDY == "2.0" ), "A"))
  small_tissue_data = small_tissue_data %>% filter(DTHHRDY == "V" | DTHHRDY == "A") %>%
    filter(SMCENTER == "B1" | SMCENTER == "C1") %>% 
    unite( 'Mixed_DTHHRDY', SMCENTER,DTHHRDY, sep = "", remove = F)
  
  if(dim(small_tissue_data)[1]==0){return(NA)}   #This is triggered if the sample has multiple testing centeres listed for SMCENTER
  #small_tissue_data = t(small_tissue_data)
  # colnames(small_tissue_data) = small_tissue_data[1, ]
  # small_tissue_data = as_tibble(small_tissue_data[-1, ])
  return(small_tissue_data)
  
} 
age_to_num = function(age_str){
  str_list = str_split(age_str, pattern = "-", n =2)
  num = mean(sapply(str_list, as.numeric))
  return(num)
}
char_to_binary = function(vect, control = "A", treatment = "V"){
  vect = matrix(vect)
  vect[grepl(control,vect )] = 0
  vect[grepl(treatment, vect)] = 1
  vect = sapply(vect, as.numeric)
  return(vect)
}

get_counts_mat= function(collated_data){
  if(is.na(collated_data)){
    return(NA)
  }
  subjects = collated_data[,1]
  counts = t(collated_data[, -c(1:10)])
  counts = data.frame(apply(counts, 2, as.numeric), check.names=F, row.names = rownames(counts))
  colnames(counts) = unlist(subjects)
  return(counts)
}

get_coldata = function(collated_data){
  if(is.na(collated_data)){
    return(NA)
  }
  subjects = collated_data[,1]
  subjects$AGE = apply(collated_data[, 9], 1, age_to_num)
  subjects$SEX = collated_data$SEX
  subjects = mutate(subjects, SEX = replace(SEX, (SEX == "1" ), "M"))
  subjects = mutate(subjects, SEX = replace(SEX, (SEX == "2" ), "F"))
  subjects$dthhrdy = (collated_data$DTHHRDY)
  subjects$center = factor(collated_data$SMCENTER)
  subjects = as.data.frame(subjects)
  subjects = unite( subjects, 'Mixed_DTHHRDY', center,dthhrdy,SEX, sep = "", remove = F)
  rownames(subjects)=subjects$SUBJID; subjects = subjects[,-1]
  return(subjects)
}

run_perms = function(file, min_in_group = 50, custom_func = T, nperm = 1000, correct = 'both', combat = F, plot = F){
  tissue = get_tissue_name(file)
  print(tissue)
  data =collate_data(file, small = F)
  coldata = get_coldata(data)
  if(is.na(coldata)){return(NA)}
  if(any(unname(table(coldata$dthhrdy)) < min_in_group)){
    return(NA)
  }
  cts = get_counts_mat(data)
  if(is.na(cts)){return(NA)}
  
  d0 <- DGEList(cts)
  d <- calcNormFactors(d0, method = "TMM")
  dthhrdy_vec = coldata$Mixed_DTHHRDY
  dthhrdy_bin = char_to_binary(coldata$dthhrdy)
  center = char_to_binary(coldata$center, control = "B1", treatment = "C1")  
  sex = char_to_binary(coldata$SEX, control = "M", treatment = "F")
  age = coldata$AGE
  
  keep <- filterByExpr(d0, group=dthhrdy_bin)          
  d <- d0[keep, , keep.lib.sizes=FALSE]
  tmm <- cpm(d)
  
  
  if(combat){
    print("Using Combat-Seq to adjust for center, gender, and age")
    if (length(unique(sex)) > 1){
      covar_mat = cbind(age, sex)
      tmm = ComBat_seq(tmm, batch =center, covar_mod = covar_mat )
    }else{
      tmm = ComBat_seq(tmm, batch =center, group = age )
      
    }
  }
  
  tmm = tmm[which(rownames(tmm) %in% gene_list),]
  emat = tmm
  
  
  if(plot){
    setwd(corr_plots)
    plotHeatmap(rownames(refCor), emat, groupVec =get_sub_str(dthhrdy_vec, start = 3, end = 3))
    ggsave(paste0(tissue,"CorrelationPlotFiltered.pdf"))
    setwd(in_directory)
  }
  if(custom_func){
    p_test = customCalcDeltaCCD(refCor, emat, nPerm = nperm , dthhrdy_vec, "A", correct = correct)
  }else{
    p_test = calcDeltaCCD(refCor, emat, nPerm = nperm, binarize_labels_cust(dthhrdy_vec, 3, 3), "A" )
  }
  p_test$group1_count = sum(!as.logical(char_to_binary(dthhrdy_vec)))
  p_test$group2_count = sum(as.logical(char_to_binary(dthhrdy_vec)))
  return(p_test)
}


get_simple_ccd = function(file, min_in_group = 50, nPerm = 1000 ){
  tissue = get_tissue_name(file)
  print(tissue)
  data =collate_data(file, small = F)
  coldata = get_coldata(data)
  if(is.na(coldata)){return(NA)}
  if(any(unname(table(coldata$dthhrdy)) < min_in_group)){
    return(NA)
  }
  cts = get_counts_mat(data)
  if(is.na(cts)){return(NA)}

  d0 <- DGEList(cts)
  d <- calcNormFactors(d0, method = "TMM")
  dthhrdy_vec = coldata$Mixed_DTHHRDY
  dthhrdy_bin = char_to_binary(coldata$dthhrdy)
  center = char_to_binary(coldata$center, control = "B1", treatment = "C1")
  sex = char_to_binary(coldata$SEX, control = "M", treatment = "F")
  age = coldata$AGE

  keep <- filterByExpr(d0, group=dthhrdy_bin)
  d <- d0[keep, , keep.lib.sizes=FALSE]
  tmm <- cpm(d)
  #tmm = tmm[which(rownames(tmm) %in% gene_list),]
  emat = tmm
  #rownames(emat) = gene_list
  #emat = emat[checkGenes(emat, refCor), ]
  #calcDeltaCCD_result = calcDeltaCCD(refCor, emat, nPerm = 1000, groupVec = coldata$dthhrdy, groupNormal = 'A' )#, dthhrdy_vec, "A", correct = correct)
  # acute_ccd= calcCCDSimple(refCor, emat[,-as.logical(dthhrdy_bin)]) #acute ccd
  # vent_ccd = calcCCDSimple(refCor, emat[,as.logical(dthhrdy_bin)]) #vent ccd
  acute_ccd= calcCCD(refCor, emat[,-as.logical(dthhrdy_bin)], nPerm = 1000)  #acute ccd
  vent_ccd= calcCCD(refCor, emat[,as.logical(dthhrdy_bin)], nPerm = 1000)  #acute ccd
  
  names= c("A", "V", "A_pval", "V_pval")
  out = c(acute_ccd$CCD,vent_ccd$CCD,  acute_ccd$Pvalue, vent_ccd$Pvalue)
  names(out) = names
  out = as.data.frame(out)
  return(out)
  
}


extract_perm_pvals = function(result, gtex_files){
  pvals = c()
  exlude_tissue_idx =c()
  for(i in 1:length(gtex_files)){
    if(!(is.na(result[[i]]))){
      pvals = rbind(pvals, result[[i]])
    }
    else{exlude_tissue_idx = c(exlude_tissue_idx, i)}
  }
  pvals = as.data.frame(pvals)
  included_gtex = get_tissue_name(gtex_files[-exlude_tissue_idx])
  rownames(pvals) = included_gtex
  return(pvals)
}

#uncomment to start saving files
# gender_shuffle = sapply(gtex_files, run_perms, nperm = 5000, correct = 'gender')
# gender_shuffle_result = extract_perm_pvals(gender_shuffle, gtex_files)
# 
# center_shuffle = sapply(gtex_files, run_perms, nperm = 5000, correct = 'center')
# center_shuffle_result = extract_perm_pvals(center_shuffle, gtex_files)
# 
# GC_shuffle = sapply(gtex_files, run_perms, nperm = 1000, correct = 'both')
# GC_shuffle_result = extract_perm_pvals(GC_shuffle, gtex_files)
# 
# no_shuffle_result = sapply(gtex_files, run_perms, custom_func = 'F')
# no_shuffle_result = extract_perm_pvals(no_shuffle_result, gtex_files)
# 
# colnames(no_shuffle_result)[4] = "no_shuffle_pval"
# no_shuffle_result$bonferroni_no_shuffled = p.adjust(no_shuffle_result$no_shuffle_pval, method = 'bonferroni', n = 25)
# 
# no_shuffle_result$gender_pval = gender_shuffle_result$Pvalue
# no_shuffle_result$bonferoni_gender_pval = p.adjust(no_shuffle_result$gender_pval, method = 'bonferroni', n = 25)
# 
# no_shuffle_result$center_pval = center_shuffle_result$Pvalue
# no_shuffle_result$bonferoni_center_pval = p.adjust(no_shuffle_result$center_pval, method = 'bonferroni', n = 25)
# 
# 
# no_shuffle_result$GC_shuffled_pval = GC_shuffle_result$Pvalue
# no_shuffle_result$bonferoni_GC_shuffled = p.adjust(no_shuffle_result$GC_shuffled_pval, method = 'bonferroni', n = 25)
# 
# 
# 
# tissue_perm_test = as_tibble(no_shuffle_result) %>%
#   select( group1, group1_count, group2, group2_count, everything())
# rownames(tissue_perm_test) = rownames(no_shuffle_result)
# 
# #### Extracting the ccd values for A and V groups####
# simple_ccds = sapply(gtex_files, get_simple_ccd, nPerm = 1000)
