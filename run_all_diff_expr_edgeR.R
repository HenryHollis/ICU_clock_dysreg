library(tidyverse)
library(RColorBrewer)
library('doParallel')
library('doRNG')
library(edgeR)
library(limma)

in_dir = "~/Documents/R/ClockCorr2/GTEx_counts/"
setwd(in_dir)
counts_dir = "~/Documents/R/ClockCorr2/2TMMNormedCountsEdgeR/"
outfile_dir = "~/Documents/R/ClockCorr2/diff_expr_tables2/"
gtex_files =gtex_files=list.files(pattern = "*.csv")
gtex_files=gtex_files[-c(7,8,9,10,11, 12, 13, 14, 15, 16, 17,
                         18, 19, 20,23, 24, 25, 26, 32, 35,
                         36, 39, 42, 43,44,48,49,50,53,54 )] #tissues don't meet inclusion criteria


gene_list = c("PER3" , "CIART" ,"NPAS2", "PER2",  "NR1D2", "CLOCK" ,"ARNTL", "CRY2" , "CRY1",  "PER1" , "NR1D1" ,"DBP","TEF")

get_tissue_name = function(str){
  tissue_name = str_extract(str, pattern = "(?<=\\()(.*)(?=\\))")
  return(tissue_name)
}


collate_data = function(file_name){
  tissue_data = read.csv(file_name)                           # read csv
  tissue_data =
    select(tissue_data, -Name) %>%
    unite("Info", Description,Info, sep = "", remove = T)
  tissue_data = t(tissue_data)
  
  if(dim(tissue_data)[1] <= 1){return(NA)}
  
  colnames(tissue_data) = tissue_data[1, ]
  tissue_data = as_tibble(tissue_data[-1, ])
  #tissue_data = filter(tissue_data, (SMCENTER =='C1' | SMCENTER == 'B1') )
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
  
  
  all_coldata = coldata
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
  keep <- rowSums(cpm(d0)>1) >= 2 #we're only keeping a gene if it has a cpm of 100 or greater for at least two samples
  d <- d0[keep,]
  d = d0
  tmm <- calcNormFactors(d, method = "TMM")
  if(write){
   setwd(counts_dir)
   write.table(tmm, file = paste0(tiss_name,"_TMMcounts.csv"), row.names = T, col.names = NA, sep = ",", quote = F)
   setwd(in_dir)
  }
  
  ######start differential expr######
  center = coldata$center
  dthhrdy = coldata$dthhrdy
  age = coldata$AGE
  sex = coldata$SEX
  single_sex = F
  if(length(unique(sex))==1){  #single sex tissue
    single_sex = T
    design <- model.matrix(~center + age + dthhrdy)  #design matrix
  }else{
    design <- model.matrix(~center + age + sex + dthhrdy)  #design matrix
  }
  
  cts = cts[,-rm_subjects]
 # d0 <- DGEList(cts, )
  #keep <- filterByExpr(d0, group=dthhrdy)
  #d <- d0[keep, , keep.lib.sizes=FALSE]
  y <- estimateDisp(d[, -rm_subjects], design)
  
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

sapply(gtex_files[2:25], tissue_diff_expr, write= T)

