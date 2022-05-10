library(tidyverse)
library(RColorBrewer)
library('doParallel')
library('doRNG')
library(edgeR)
library(limma)
library(openxlsx)
library(sva)
library(doParallel)

cl <- makePSOCKcluster(2)
registerDoParallel(cl)
in_dir = "~/Documents/R/ClockCorr2/GTEx_counts/"
setwd(in_dir)
counts_dir = "~/Documents/R/ClockCorr2/TMMNormedCountsEdgeR_clockcorr3/"
outfile_dir = "~/Documents/R/ClockCorr2/diff_expr_tables2/"
gtex_files =gtex_files=list.files(pattern = "*.csv")
gtex_files=gtex_files[-c(7,8,9,10,11, 12, 13, 14, 15, 16, 17,
                         18, 19, 20,23, 24, 25, 26, 32, 35,
                         36, 39, 42, 43,44,48,49,50,53,54 )] #tissues don't meet inclusion criteria

get_tissue_name = function(str){
  tissue_name = str_extract(str, pattern = "(?<=\\()(.*)(?=\\))")
  return(tissue_name)
}


collate_data = function(file_name){
  tissue_data = read.csv(file_name)      # read csv
  tissue_data =
    select(tissue_data, -Name) %>%
    unite("Info", Description,Info, sep = "", remove = T)
  tissue_data = t(tissue_data)
  
  if(dim(tissue_data)[1] <= 1){return(NA)}
  
  colnames(tissue_data) = tissue_data[1, ]
  tissue_data = as_tibble(tissue_data[-1, ])
  if(dim(tissue_data)[1]==0){return(NA)}   #This is triggered if the sample has multiple testing centers listed for SMCENTER
  
  return(tissue_data)
  
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

get_counts_mat = function(collated_data){
  #Returns only counts from data, which had meta data in it
  if(is.na(collated_data)){
    return(NA)
  }
  subjects = collated_data[,1]
  counts = t(collated_data[, -c(1:9)])
  counts = data.frame(apply(counts, 2, as.numeric), check.names=F, row.names = rownames(counts))
  colnames(counts) = unlist(subjects)
  return(counts)
}

get_coldata = function(collated_data){
  if(is.na(collated_data)){
    return(NA)
  }
  subjects = collated_data[,1]
  subjects$AGE = apply(collated_data[, 8], 1, age_to_num)
  subjects$SEX = factor(collated_data$SEX)
  subjects$dthhrdy = factor(char_to_binary(t(collated_data$DTHHRDY)))
  subjects$center = (collated_data$SMCENTER)
  subjects = as.data.frame(subjects)
  rownames(subjects)=subjects$SUBJID; subjects = subjects[,-1]
  return(subjects)
}

tissue_diff_expr = function(file, min_in_group=50, write = F){
  tiss_name = get_tissue_name(file)
  print(tiss_name)
  tiss_name = str_replace_all(tiss_name, fixed(" "), "")
  start = Sys.time()
  data = collate_data(file)
  coldata = get_coldata(data)
  rm_subjects = which(!(coldata$dthhrdy %in% c(0, 1, 2) & coldata$center %in% c("B1", "C1")))
  
  coldata = coldata[-rm_subjects, ]
  coldata$dthhrdy = as.character(coldata$dthhrdy)
  coldata$dthhrdy[coldata$dthhrdy %in% c('1', '2') ] = "A"
  coldata$dthhrdy[coldata$dthhrdy == '0' ] = "V"
  
  if(all(is.na(coldata))){return(NA)}
  if(any(table(coldata$dthhrdy) < min_in_group)){
    return(NA)
  }
  cts = get_counts_mat(data)
  
  if(is.na(cts)){return(NA)}
  d0 <- DGEList(cts)
  
  #######norming counts to TMMS#####
  keep <- rowSums(cpm(d0)>1) >= 2 #we're only keeping a gene if it has a cpm of 1 or greater for at least two samples
  d <- d0[keep,]
  
  tmm <- calcNormFactors(d, method = "TMM")
  out_tmm = cpm(tmm) #cpm aware of tmm normalization: https://www.biostars.org/p/317701/
  
  
  center = coldata$center
  dthhrdy = coldata$dthhrdy
  age = coldata$AGE
  sex = coldata$SEX
  single_sex = F
  
  if(length(unique(sex))==1){  #single sex tissue
    single_sex = T
    combat_data = ComBat_seq(as.matrix(cts[,-rm_subjects]), batch =center, group = age ) #only one sex
  }else{ #two sex tissue
    covar_mat = cbind(as.numeric(age), as.numeric(coldata$SEX))
    combat_data = ComBat_seq(as.matrix(cts[, -rm_subjects]), batch =center, covar_mod = covar_mat )
  }
  keep <- rowSums(cpm(combat_data)>1) >= 2 #we're only keeping a gene if it has a cpm of 1 or greater for at least two samples
  combat_data_filt <- combat_data[keep,]
  combat_data_filt = DGEList(combat_data_filt)
  combat_tmm <- calcNormFactors(combat_data_filt, method = "TMM")
  combat_tmm = cpm(combat_tmm) #cpm aware of tmm normalization: https://www.biostars.org/p/317701/
  
  
  
  if(write){
   setwd(counts_dir)
   list_of_datasets <- list("TMM" = as.data.frame(out_tmm), "subjectInfo" = as.data.frame(coldata), "combat_corrected" = as.data.frame(combat_tmm))
   tiss_name = str_remove(tiss_name, '-')
   write.xlsx(list_of_datasets, file = paste0(tiss_name, "_TMM_and_subject_info.xlsx"), asTable = T, rowNames=T)
   setwd(in_dir)
  }
  
  ######start differential expr######

  if(single_sex){  #single sex tissue
    single_sex = T
    design <- model.matrix(~center + age + dthhrdy)  #design matrix
  }else{
    design <- model.matrix(~center + age + sex + dthhrdy)  #design matrix
  }

  tmm = tmm[,-rm_subjects]
  y <- estimateDisp(tmm, design)

  fit = glmQLFit(y, design)
  if(single_sex){
    qlf = glmQLFTest(fit, coef=4)
  }else{
    qlf = glmQLFTest(fit, coef=5)
  }
  top.table <- topTags(qlf, n = Inf)

  stop = Sys.time()
  print( difftime(stop,start, units = 'auto'))
  print(paste(length(which(unlist(top.table[,4]) < 0.05)), "genes identified"))
  if(write){
    setwd(outfile_dir)
    write.table(top.table, file = paste0(tiss_name,"_dthhrdy_diff_expr_table_edgeR.csv"), row.names = T, col.names = NA, sep = ",", quote = F)
    setwd(in_dir)
  }
  return(top.table)
}

sapply(gtex_files, tissue_diff_expr, write= T)

