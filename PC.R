library(pcalg)


# skeleton determination of PC
# if type == discrete, consider G^2 test
# if type == Gaussian, consider Fisher-z
skeleton_pc <- function(data, alpha, method = ("original"), type){
  # type: 'discrete' or 'Gaussian'
  
  n_node <- ncol(data)
  n_score_count <- 0
  skeleton <- matrix(TRUE, nrow = n_node, ncol = n_node)
  diag(skeleton) <- FALSE
  
  if(type == 'Gaussian'){
    cor_mat <- cor(data)
  }
  
  node_seq <- seq_len(n_node)  # now the element means node, and position as well
  d_sep_set <- lapply(node_seq, function(.) vector('list', n_node))
  
  l <- -1
  done <- FALSE
  
  while(!done & any(skeleton)){
    done <- TRUE
    # get the incide of adjacent sets
    ind <- which(skeleton, arr.ind = T)
    ind <- ind[order(ind[, 1]), ]
    n_adj <- nrow(ind)
    # adj for stable
    if(method == 'stable'){
      adj_stable <- split(skeleton, gl(n_node, n_node))
    }
    
    l <- l + 1
    
    for(i in 1:n_adj){
      node_i <- ind[i, 1]
      node_j <- ind[i, 2]
      if(method == 'stable'){
        adj_i <- adj_stable[[node_i]]
      }else{
        adj_i <- skeleton[node_i, ]
      }
      adj_i[node_j] <- FALSE  # delete j
      adj_i <- node_seq[adj_i]
      
      if(length(adj_i) >= l){
        if(length(adj_i) > l) done <- FALSE
        
        K_set <- combn(adj_i, m = l) # obtain all the possible candidate sets
        K_set <- t(K_set)
        if(length(adj_i) == 1) K_set <- as.matrix(adj_i) # when length(adj_i) == 1, combn() consider the input as seq_len(.)
        
        for(j in 1:nrow(K_set)){
          K <- K_set[j, ]
          if(type == 'discrete'){
            p_value <- gSquareDis(node_i, node_j, S = K, dm = as.matrix(data))
          }else if(type == 'Gaussian'){
            suffStat <- list(C = cor_mat, n = nrow(data))
            p_value <- gaussCItest(node_i, node_j, S = K, suffStat)
          }
          
          n_score_count <- n_score_count + 1
          if(p_value >= alpha){
            skeleton[node_i, node_j] <- skeleton[node_j, node_i] <- F
            d_sep_set[[node_i]][[node_j]] <- K
            break
          }else{
            if(length(K) == 0) break  # if K is empty
            if(j == nrow(K_set)) break   # if K is the last set
          }
        }
      }
    }
  }
  result <- list(skeleton = skeleton, d_sep_set = d_sep_set,
                 n_score_count = n_score_count)
  return(result)
}



# orientation determination of PC

orientate_pc <- function(result_pc){
  # result_pc is the output of skeleton_pc
  # contains the skeleton and the conditional independent sets
  
  skeleton <- result_pc$skeleton
  d_sep_set <- result_pc$d_sep_set
  
  PDAG <- skeleton + 0
  
  # determine the v-structure
  ind <- which(PDAG == 1, arr.ind = T)
  for(m in 1:nrow(ind)){
    node_i <- ind[m, 1]
    node_k <- ind[m, 2]
    ind_j <- which(PDAG[node_k, ] == 1)
    ind_j <- setdiff(ind_j, node_i)
    for(node_j in ind_j){
      # if i, j are not adjacent
      # if k is not the Sepset of i and j
      # then v-structure, i -> k <- j
      if(!(node_k %in% d_sep_set[[node_i]][[node_j]] | 
           node_k %in% d_sep_set[[node_j]][[node_i]]) &
         (PDAG[node_i, node_j] == 0 & PDAG[node_j, node_i] == 0)){
        PDAG[node_i, node_k] <- PDAG[node_j, node_k] <- 1
        PDAG[node_k, node_i] <- PDAG[node_k, node_j] <- 0
      }
    }
  }
  
  
  ind <- which(PDAG == 1 & t(PDAG) == 1, arr.ind = TRUE)
  if(length(ind) > 0){
    for(m in 1:nrow(ind)){
      # rule 1
      node_j <- ind[m, 1]
      node_k <- ind[m, 2]
      ind_i <- which(PDAG[, node_j] == 1 & PDAG[node_j, ] == 0)
      for(node_i in ind_i){
        if(PDAG[node_i, node_k] == 0 & PDAG[node_k, node_i] == 0){
          PDAG[node_j, node_k] <- 1
          PDAG[node_k, node_j] <- 0
        }
      }
      
      # rule 2
      node_i <- ind[m, 1]
      node_j <- ind[m, 2]
      ind_k <- which(PDAG[node_i, ] == 1 & PDAG[, node_i] == 0 &
                       PDAG[node_j, ] == 0 & PDAG[, node_j] == 1)
      if(length(ind_k) > 0){
        PDAG[node_i, node_j] <- 1
        PDAG[node_j, node_i] <- 0
      }
      
      # rule 3
      node_i <- ind[m, 1]
      node_j <- ind[m, 2]
      ind_lk <- which(PDAG[node_i, ] == 1 & PDAG[, node_i] == 1 &
                        PDAG[node_j, ] == 1 & PDAG[, node_j] == 0)
      if(length(ind_lk) >= 2){
        ind_lk <- t(combn(ind_lk, 2))    
        for(z in 1:nrow(ind_lk)){
          node_k <- ind_lk[z, 1]
          node_l <- ind_lk[z, 2]
          if(PDAG[node_k, node_l] == 0 & PDAG[node_l, node_k] == 0){
            PDAG[node_i, node_j] <- 1
            PDAG[node_j, node_i] <- 0
          }
        }
      }
    }
  }
  
  return(PDAG)
}
