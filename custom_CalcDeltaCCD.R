#source("~/Documents/R/ClockCorrelation/deltaccd_utils.R")
library('deltaccd')
library('doParallel')
library('doRNG')

customCalcDeltaCCD = function(
  refCor, emat, groupVec, groupNormal, correct = "both", refEmat = NULL, nPerm = 1000,
  geneNames = NULL, dopar = FALSE, scale = FALSE) {
  group2Now = i = NULL
  
  method = 'spearman'
  doOp = if (isTRUE(dopar)) `%dorng%` else `%do%`
  
  refCor = checkRefCor(refCor, refEmat, geneNames, method)
  geneNames = checkGenes(emat, refCor)
  print("Runnning customCalcDeltaCCD")
  if (length(groupVec) != ncol(emat)) {
    stop('Length of groupVec does not match the number of columns in emat.')
  } #else if (!(groupNormal %in% groupVec)) {
  #stop('The supplied value for groupNormal is not present in groupVec.')
  #} 
  # else {
  #   tt = table(groupVec)
  #   if (length(tt) < 2) {
  #     stop('groupVec contains only one unique group.')
  #   } else if (min(tt) < 3) {
  #     stop('Each unique group in groupVec must have at least three samples.')}}
  # 
  checkVar(emat[geneNames, ], groupVec)
  
  binaryGroupVec = binarize_labels_cust(groupVec, start = 3, end = 3)
  
  result = data.table(
    group1 = groupNormal,                                   # Acute group
    group2 = setdiff(sort(unique(binaryGroupVec)), groupNormal))  # Vent group
  
  resultTmp = foreach(i=1:1, .combine = rbind) %do% {
    idx1 = binaryGroupVec %in% c("A", "V")    # true vector of length(dthhrdy)
    idx2 = binaryGroupVec[idx1] == "V"        # boolean, true where ventilator group  
    
    if(correct == "both"){
      BM_indx = which(grepl("B1(A|V)M", groupVec))  # data indices gathered from center B1 and male
      CM_indx = which(grepl("C1(A|V)M", groupVec))  # indices gathered from center C1 and male
      BF_indx = which(grepl("B1(A|V)F", groupVec))  # B1 and female
      CF_indx = which(grepl("C1(A|V)F", groupVec))  # C1 and female
    }else if (correct == "center"){
      B_indx = which(grepl("B1", groupVec))  # data indices gathered from center B1 
      C_indx = which(grepl("C1", groupVec))  # indices gathered from center C1
    }else if (correct == "gender"){
      M_indx = which(grepl("M", groupVec))  # data indices gathered from males
      F_indx = which(grepl("F", groupVec))  # indices gathered from females
    }
    
    ematNow = emat[geneNames, idx1]
    
    deltaCcdObs = calcDeltaCCDSimple(     # returns scalar deltaCCD
      refCor, ematNow, idx2, method = method, scale = scale)
    if (correct == "both"){
      idxPermBM = makePerms(idx2[BM_indx], nPerm = nPerm, dopar = dopar)           # shuffle labels among samples measured from B1 center males
      idxPermCM = makePerms(idx2[CM_indx], nPerm = nPerm, dopar = dopar)           # shuffle labels among samples measured from C1 center males
      idxPermBF = makePerms(idx2[BF_indx], nPerm = nPerm, dopar = dopar)           # B1 females
      idxPermCF = makePerms(idx2[CF_indx], nPerm = nPerm, dopar = dopar)           # C1 females
      if(!any(is.na(c(idxPermBM, idxPermCM, idxPermBF, idxPermCF)) )){
        
        perms = as.data.frame(matrix(0, nrow = nPerm, ncol = length(groupVec)))    # create a matrix of zeros to fill in next step
        perms[,BM_indx] = idxPermBM                                                 # populate "b1"-cols of perms 
        perms[, CM_indx] = idxPermCM                                                 # populate "c1"-cols of perms
        perms[, BF_indx] = idxPermBF                                                 # populate "c1"-cols of perms
        perms[, CF_indx] = idxPermCF                                                 # populate "c1"-cols of perms
      }else{
        return(data.table(DeltaCCD = deltaCcdObs, Pvalue = NA))  #add p-value stats to results table
      }
      
    }else if (correct == "center"){
      idxPermB = makePerms(idx2[B_indx], nPerm = nPerm, dopar = dopar)           # 
      idxPermC = makePerms(idx2[C_indx], nPerm = nPerm, dopar = dopar)           # 
      perms = as.data.frame(matrix(0, nrow = nPerm, ncol = length(groupVec)))    # 
      perms[,B_indx] = idxPermB                                                #  
      perms[, C_indx] = idxPermC                                                 #
      
      
    }else if (correct == "gender"){
      idxPermM = makePerms(idx2[M_indx], nPerm = nPerm, dopar = dopar)           # 
      idxPermF = makePerms(idx2[F_indx], nPerm = nPerm, dopar = dopar)           # 
      perms = as.data.frame(matrix(0, nrow = nPerm, ncol = length(groupVec)))    # 
      perms[,M_indx] = idxPermM                                                  #  
      perms[, F_indx] = idxPermF                                                 # 
    }
    
    perms = apply(perms, 2, as.logical)
    #perms = makePerms(idx2, nPerm = nPerm, dopar = dopar)
    
    deltaCcdRand = doOp(foreach(i = 1:nPerm, .combine = c), { #1000 vector, deltaCCD for all 1000 perms
      calcDeltaCCDSimple(
        refCor, ematNow, perms[i,], method = method, scale = scale)})
    
    #nComb = choose(length(idx2), sum(idx2))
    if(correct == "both"){
      nComb = choose(length(BM_indx), sum(idx2[BM_indx]) ) * choose(length(CM_indx), sum(idx2[CM_indx])) * choose(length(BF_indx), sum(idx2[BF_indx]))* choose(length(CF_indx), sum(idx2[CF_indx]))
    }else if (correct == "center"){
      nComb = choose(length(B_indx), sum(idx2[B_indx]) ) * choose(length(C_indx), sum(idx2[C_indx]))
    }else if (correct == "gender"){
      nComb = choose(length(M_indx), sum(idx2[M_indx]) ) * choose(length(F_indx), sum(idx2[F_indx]))
    }
    pvalue = statmod::permp(      #how many out of your nperms had DeltaCCD >= true labeled DeltaCCD
      sum(deltaCcdRand >= deltaCcdObs), nperm = nPerm, total.nperm = nComb,
      twosided = FALSE, method = 'approximate')
    
    data.table(DeltaCCD = deltaCcdObs, Pvalue = pvalue)  #add p-value stats to results table
    
  }
  
  result = cbind(result, resultTmp)
  return(result)
}

