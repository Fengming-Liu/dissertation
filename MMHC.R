MMPC <- function(data, alpha = .05, type = 'Gaussian'){
  p <- ncol(data)
  skeleton <- matrix(0, ncol = p, nrow = p)
  n_score_cal <- 0
  n_forward <- 0
  n_backward <- 0
  
  # data transformation for discrete
  data <- as.data.frame(data)
  if(type == 'discrete'){
    data[] <- lapply(data, factor)
  }
  
  
  # forward
  CPC <- rep(list(NULL), p)  # candidate parent and child set
  dep_set <- matrix(rep(seq(1, p), p), nrow = p, byrow = T)
  diag(dep_set) <- NA
  
  for(i in 1:p){
    dep_set_i <- na.omit(dep_set[i, ])
    CPC_i <- CPC[[i]]
    p_store <- matrix(alpha, ncol = length(dep_set_i))
    colnames(p_store) <- dep_set_i
    p_store <- p_store[-1,]
    set_test <- c()
    
    while(length(dep_set_i) != 0){
      n_forward <- n_forward + 1
      
      p_store_iter <- c()
      if(is.null(CPC_i)){
        for(j in dep_set_i){
          if(type == 'discrete'){
            p_value <- ci.test(data[, i], data[, j], test = 'x2')$p.value
          }else if(type == 'Gaussian'){
            p_value <- ci.test(data[, i], data[, j], test = 'zf')$p.value
          }
          
          if(p_value > alpha){  # independent, pop out
            p_store <- p_store[, -which(colnames(p_store) == j), drop = F]
            dep_set_i <- dep_set_i[-which(dep_set_i == j)]
          }else{  # dependent, store the p
            p_store_iter <- c(p_store_iter, p_value)
          }
        }
        if(length(p_store_iter) != 0){
          p_store_iter <- as.matrix(t(p_store_iter))
          name_p_store <- c(rownames(p_store), 'independent test')
          p_store <- rbind(p_store, p_store_iter)
          rownames(p_store) <- name_p_store
          set_test <- name_p_store
        }
        
      }else{
        for(m in 1:length(CPC_i)){
          if(length(CPC_i) == 1){
            comb_set <- CPC_i
            comb_set <- as.matrix(comb_set)
          }else{
            comb_set <- combn(as.vector(CPC_i), m)
            comb_set <- as.matrix(comb_set)
          }
          for(k in 1:ncol(comb_set)){
            p_store_iter <- c()
            CPC_test <- comb_set[, k]
            name_CPC_test <- paste(CPC_test, collapse = '')
            if(!(name_CPC_test %in% set_test)){
              for(j in dep_set_i){
                if(type == 'discrete'){
                  p_value <- ci.test(data[, i], data[, j], data[, CPC_test], test = 'x2')$p.value
                }else if(type == 'Gaussian'){
                  p_value <- ci.test(data[, i], data[, j], data[, CPC_test], test = 'zf')$p.value
                }
                
                if(p_value > alpha){  # independent, pop out
                  p_store <- p_store[, -which(colnames(p_store) == j), drop = F]
                  dep_set_i <- dep_set_i[-which(dep_set_i == j)]
                }else{  # dependent, store the p
                  p_store_iter <- c(p_store_iter, p_value)
                }
              }
              if(length(p_store_iter) != 0){
                p_store_iter <- as.matrix(t(p_store_iter))
                name_p_store <- c(rownames(p_store), name_CPC_test)
                p_store <- as.matrix(rbind(p_store, p_store_iter))
                rownames(p_store) <- name_p_store
                set_test <- name_p_store
              }
            }
          }
        }
        
      }
      
      
      # compare and update
      if(nrow(p_store) != 0){
        # the node with the min association
        node_col_max <- apply(p_store, 2, max)
        node_min_ass <- as.numeric(colnames(p_store)[which.min(node_col_max)])
        CPC_i <- c(CPC_i, node_min_ass)
        p_store <- p_store[, -which(dep_set_i == node_min_ass), drop = F]
        dep_set_i <- dep_set_i[-which(dep_set_i == node_min_ass)]
        dep_set[i, node_min_ass] <- NA
      }else{
        break  # if all dep_set remain are independent
      }
      
      
      if(!is.null(CPC_i)){
        CPC[[i]] <- CPC_i
      }
    }
  }
  # backward
  for(i in 1:p){
    CPC_i <- CPC[[i]]
    if(!is.null(CPC_i)){
      for(j in CPC_i){
        
        dsep_node <- CPC_i[-which(CPC_i == j)]
        if(length(dsep_node) != 0){
          for(m in 1:length(dsep_node)){
            n_backward <- n_backward + 1
            if(length(dsep_node) == 1){
              dsep_set <- dsep_node
            }else{
              dsep_set <- combn(dsep_node, m = m)
            }
            
            n_score_cal <- n_score_cal + ncol(as.matrix(dep_set))
            if(type == 'discrete'){
              p_value <- apply(as.matrix(dsep_set), 2, function(x) ci.test(data[, i], data[, j], data[, x], test = 'x2')$p.value)
            }else if(type == 'Gaussian'){
              p_value <- apply(as.matrix(dsep_set), 2, function(x) ci.test(data[, i], data[, j], data[, x], test = 'zf')$p.value)
            }
            
            p_test <- p_value > alpha
            if(sum(p_test) != 0){
              CPC_i <- CPC_i[-which(CPC_i == j)]
            }
          }
        }
      }
      CPC[[i]] <- CPC_i
    }
  }
  
  # turn the CPC into skeleton
  for(i in 1:length(CPC)){
    CPC_i <- CPC[[i]]
    skeleton[i, CPC_i] <- 1
  }
  skeleton <- skeleton + t(skeleton)
  skeleton[skeleton != 2] <- 0
  skeleton[skeleton == 2] <- 1
  
  result <- list(skeleton = skeleton, n_score_count = n_score_cal,
                 n_forward = n_forward, n_backward = n_backward)
  return(result)
}



MMHC <- function(data, alpha = .05, type){
  result_mmpc <- MMPC(data, alpha = alpha, type = 'Gaussian')
  skeleton <- result_mmpc$skeleton
  
  result_hc <- hill_climbing(data, skeleton = skeleton, 
                             type = 'Gaussian')
  graph <- result_hc$graph
  
  result <- list(graph = graph, skeleton = skeleton,
                 mmpc_n_score_count = result_mmpc$n_score_count,
                 mmpc_n_forward = result_mmpc$n_forward,
                 mmpc_n_backward = result_mmpc$n_backward,
                 hc_change = result_hc$change,
                 hc_n_mt_update = result_hc$n_mt_update,
                 hc_n_BIC_update = result_hc$n_BIC_update)
  return(result)
}