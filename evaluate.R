evaluate <- function(g1, skel, g2){
  # g1: true DAG
  # skel: estimated skeleton
  # g2: estimated DAG
  
  ### skeleton
  skel_T <- g1 + t(g1)
  skel_T[skel_T == 2] <- 1
  n_T <- sum(skel_T) / 2
  n <- sum(skel) / 2
  
  # TDR: true discovery rate
  # the proportion of true edges discovered in the real network
  TDR_ind <- (skel_T == 1) & (skel == 1)
  TDR <- sum(TDR_ind) / 2 / n_T
  
  # missing/extra edges and rate
  miss_ind <- (skel_T == 1) & (skel == 0)
  miss_num <- sum(miss_ind) / 2
  extra_ind <- (skel_T == 0) & (skel == 1)
  extra_num <- sum(extra_ind) / 2
  miss_extra_num <- miss_num + extra_num
  miss_rate <- miss_num / n_T
  extra_rate <- extra_num / n
  
  
  ### DAG
  n_T <- sum(g1)
  n <- sum(g2)
  
  # DAG: SHD, SHD rate
  shd <- 0
  s1 <- g1 + t(g1)
  s2 <- g2 + t(g2)
  s1[s1 == 2] <- 1
  s2[s2 == 2] <- 1
  ds <- s1 - s2
  ind <- which(ds > 0)
  g1[ind] <- 0
  shd <- shd + length(ind)/2  # missing
  ind <- which(ds < 0)
  g1[ind] <- g2[ind]
  shd <- shd + length(ind)/2  # extra
  d <- abs(g1 - g2)
  shd <- shd + sum((d + t(d)) > 0)/2  # reverse
  shd_rate <- shd / n_T
  
  
  ## F1_TPR
  # TN: true negative
  g1_0 <- g1 == 0
  g2_0 <- g2 == 0
  TN_ind <- g1_0 & g2_0
  TN <- (sum(TN_ind) - length(diag(g1)))
  
  # TP: true positive
  g1_1 <- g1 == 1
  g2_1 <- g2 == 1
  TP_ind <- g1_1 & g2_1
  TP <- sum(TP_ind)
  
  # FP: false positive (re + miss), FN: false negative (extra)
  FN_ind <- g1_0 & g2_1  # extra + re
  FP_ind <- g1_1 & g2_0  # miss + re
  
  re_ind <- t(FN_ind) & FP_ind
  FP <- sum(re_ind)
  re_ind <- re_ind | t(re_ind)
  
  FN_ind <- FN_ind & (!re_ind)
  FP_ind <- FP_ind & (!re_ind)
  
  FN <- sum(FN_ind)
  FP <- sum(FP_ind) + FP
  
  # recall, TPR, precision
  recall <- TP / (TP + FN)
  precision <- TP / (TP + FP)
  F1 <- 2 * (recall * precision) / (recall + precision)
  
  # TRP, FPR
  TPR <- recall
  FPR <- FP / (FP + TN)
  
  result <- list('TDR' = TDR, 
                 'miss_edge_number' = miss_num,
                 'extra_edge_number' = extra_num, 
                 'miss_extra_number' = miss_extra_num,
                 'miss_rate' = miss_rate,
                 'extra_rate' = extra_rate,
                 'SHD' = shd,
                 'SHD_rate' = shd_rate,
                 'recall' = recall,
                 'precision' = precision,
                 'F1' = F1,
                 'TPR' = TPR,
                 'FPR' = FPR)
  return(result)
}
