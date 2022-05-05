library(tidyverse)
library(doParallel)

registerDoParallel(cores = 4)
in_directory = "~/Documents/R/ClockCorr2/GTEx_counts/"

setwd(in_directory)
gtex_files=list.files(pattern = "*.csv")
gtex_files=gtex_files[-c(7,8,9,10,11, 12, 13, 14, 15, 16, 17, 18, 19, 20,23, 24, 25, 26, 32, 35, 36, 39, 42, 43,44,48,49,50,53,54 )]
get_tissue_name = function(str){
  tissue_name = str_extract(str, pattern = "(?<=\\()(.*)(?=\\))")
  return(tissue_name)
}
collate_data = function(file_name, small = T){
  tissue_data = read.csv(file_name)                           # read csv
  keep = c(1:9)
  small_tissue_data = tissue_data[keep,]

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
  mat = sapply(str_list, as.numeric)
  num = (colMeans(mat))
  return(num)
}

age_to_num_med = function(age_str){
  str_list = str_split(age_str, pattern = "-", n =2)
  mat = sapply(str_list, as.numeric)
  num = median(colMeans(mat))
  return(num)
}

age_to_num_mad = function(age_str){
  str_list = str_split(age_str, pattern = "-", n =2)
  mat = sapply(str_list, as.numeric)
  num = mad(colMeans(mat))
  return(num)
}

out = foreach(i = 1:25, .combine = rbind) %dopar%{
  data = collate_data(gtex_files[i], small = T)
  print(get_tissue_name(gtex_files[i]))
  a = data %>%
    count(DTHHRDY)%>%print()#, SMCENTER) %>%
  b = data %>%
    count(SEX)%>%print()
  c = data %>%
    count(SMCENTER)%>%print()
    #group_by(SMCENTER) %>%
    #mutate(prop = prop.table(n)) %>% print()
  age = sapply(data$AGE, age_to_num_med)
  d = round(mean(age), 0)
  e = round(dim(filter(data, SEX == 2 & DTHHRDY == 'V'))[1]/dim(filter(data, DTHHRDY == 'V'))[1], 3)*100   #percent female in vent
  f = round(dim(filter(data, SEX == 2 & DTHHRDY == 'A'))[1]/dim(filter(data, DTHHRDY == 'A'))[1], 3)*100   #percent female in acute
  g = round(dim(filter(data, SEX == 2 & SMCENTER == 'B1'))[1]/dim(filter(data, SMCENTER == 'B1'))[1], 3)*100   #percent female in B1
  h = round(dim(filter(data, SEX == 2 & SMCENTER == 'C1'))[1]/dim(filter(data, SMCENTER == 'C1'))[1], 3)*100   #percent female in C1
  i = round(sapply(filter(data, DTHHRDY == 'A' ) %>% select(AGE), age_to_num_med), 0) #median Acute age
  z = round(sapply(filter(data, DTHHRDY == 'A' ) %>% select(AGE), age_to_num_mad), 2) #mad age for acute 

  j = round(sapply(filter(data, DTHHRDY == 'V' ) %>% select(AGE), age_to_num_med), 0) #median vent age
  x = round(sapply(filter(data, DTHHRDY == 'V' ) %>% select(AGE), age_to_num_mad), 2) #mad age for Vent 
  
  k = round(sapply(filter(data, SMCENTER == 'B1' ) %>% select(AGE), age_to_num_med), 0) #median b1 age
  l = round(sapply(filter(data, SMCENTER == 'C1' ) %>% select(AGE), age_to_num_med), 0) #median c1 age
  m = round(dim(filter(data, SMCENTER == 'B1' & DTHHRDY == 'A'))[1]/dim(filter(data, DTHHRDY == 'A'))[1], 3)*100   #percent B1 in vent
  n = round(dim(filter(data, SMCENTER == 'B1' & DTHHRDY == 'V'))[1]/dim(filter(data, DTHHRDY == 'V'))[1], 3)*100   #percent B1 in vent
  
  c(a$n, b$n, c$n, d, e, f, g, h, i,z, j,x,  k, l, m, n)
}

colnames(out) = c("A", "V", "1" , "2", "B1","C1", "Median Age", "%F_in_V", "%F_in_A", "%F_in_B1", "%F_in_C1", "median_A_age","MAD A age", "median_vent_age","MAD V age", "mean_B1_age", "mean_C1_age", "%B1_in_A", "%B1_in_ICU")
rownames(out) = sapply(gtex_files, get_tissue_name )
write.table(out, "../patientStats2.csv", sep=',', quote = F, col.names = NA)
