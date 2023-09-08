# discrete BIC calculation
BIC_selfwritten <- function(node, graph, data){
  # obtain the levels of each v.r.
  level <- apply(data, 2, unique)
  level <- unlist(as.numeric(apply(as.matrix(level), 2, max), drop = F))
  
  pr <- which(graph[, node] == 1)  # parents
  pr_lv <- level[pr]
  if(length(pr) == 0){
    q_i <- 1
  }else{
    q_i <- prod(level[pr]) # number of parents combination
  }
  r_i <- level[node]
  qr <- q_i * r_i
  
  # generate the combination of node with each parent
  prod_each <- 1
  comb_pr_node <- data.frame()
  pr_node_lv <- c(pr_lv, r_i)
  for(i in pr_node_lv){
    comb_each <- rep(1:i, each = prod_each, length.out = qr)
    comb_pr_node <- rbind(comb_pr_node, comb_each)
    prod_each <- prod_each * i
  }
  comb_pr_node <- t(comb_pr_node)
  
  # calculate the number of each combination
  N_ijk <- c()
  data_pr_node <- data[, c(pr, node)]
  for(i in 1:nrow(comb_pr_node)){
    num_comb_each <- apply(as.matrix(data_pr_node), 1, function(x) all(x == comb_pr_node[i, ]))
    N_ijk <- c(N_ijk, sum(num_comb_each))
  }
  comb_pr_node <- cbind(comb_pr_node, N_ijk)
  
  
  # generate the combination of parens
  if(length(pr) == 0){
    N_ij <- nrow(data)
  }else{
    prod_each <- 1
    comb_pr <- data.frame()
    for(i in pr_lv){
      comb_each <- rep(1:i, each = prod_each, length.out = q_i)
      comb_pr <- rbind(comb_pr, comb_each)
      prod_each <- prod_each * i
    }
    comb_pr <- t(comb_pr)
    
    # calculate the number of each combination
    n_pr <- length(pr)
    N_ij <- c()
    for(i in 1:nrow(comb_pr)){
      pos <- which(apply(as.matrix(comb_pr_node[, 1:n_pr]), 1, function(x) all(x == comb_pr[i, ])))
      num_comb_each <- sum(N_ijk[pos])
      N_ij <- c(N_ij, num_comb_each)
    }
  }
  
  # avoid the 0
  N_ijk <- N_ijk + 1e-20
  N_ij <- N_ij + 1e-20
  
  cal1 <- sum(N_ijk * log(N_ijk / N_ij))
  cal2 <- log(nrow(data)) / 2 * (r_i - 1) * q_i
  BIC <- cal1 - cal2
  
  return(BIC)
}



# continuous BIC calculation
BIC_g <- function(node, graph, data){
  data <- as.data.frame(data)
  name_data <- names(data)
  name_node <- name_data[node]
  n <- nrow(data)
  p <- ncol(data)
  
  # obtain parents of node
  pr <- which(graph[, node] == 1)
  name_pr <- name_data[pr]
  
  if(length(pr) == 0){
    logLike <- with(data, sum(dnorm(data[,node], mean(data[,node]),
                                    sd(data[,node]), log=TRUE)))
    penalty <- 0.5 * log(n) * 2
    
  }else{
    forml <- paste(name_pr, collapse = '+')
    forml <- paste(name_node,'~',forml, collapse = '')
    forml <- as.formula(forml)
    model <- lm(forml, as.data.frame(data))
    logLike <- with(data, sum(dnorm(data[,node],
                                    fitted(model), 
                                    stats::sigma(model), 
                                    log=TRUE)))
    penalty <- 0.5 * log(n) * (2 + length(pr))
  }
  
  BIC_g <- logLike - penalty
  
  return(BIC_g)
}


# cycle detection
topological_sorting <- function(graph_test){
  while(1){
    del <- which(apply(as.matrix(graph_test), 2, sum) == 0)
    if(length(del) == 0){
      ifCycle <- TRUE
      break
    }
    graph_test <- graph_test[-del, -del]
    graph_test <- as.matrix(graph_test)
    if(nrow(graph_test) == 0){
      ifCycle <- FALSE
      break
    }
  }
  return(ifCycle)
}


hill_climbing <- function(data, type = 'Gaussian', skeleton = NULL){
  p <- ncol(data)
  n <- nrow(data)
  graph <- matrix(0, nrow = p, ncol = p)
  n_BIC_update <- 0 # count the number of BIC update
  
  if(is.null(skeleton)){ # when no skeleton is provided in advance
    add_diff_mat <- matrix(0, nrow = p, ncol = p)
    add_diff_mat[upper.tri(add_diff_mat)] <- 1
  }else{
    add_diff_mat <- skeleton
    add_diff_mat[lower.tri(add_diff_mat)] <- 0
    blacklist_diff <- skeleton == 0
  }
  
  for(i in 1:(p-1)){
    add_set <- which(add_diff_mat[i, ] == 1)
    for(j in add_set){
      add_graph <- graph
      add_graph[i, j] <- 1
      if(type == 'discrete'){
        add_score <- BIC_selfwritten(node = j, graph = add_graph, data)
        add_diff <- add_score - BIC_selfwritten(node = j, graph = graph, data)
      }else if(type == 'Gaussian'){
        add_score <- BIC_g(node = j, graph = add_graph, data)
        add_diff <- add_score - BIC_g(node = j, graph = graph, data)
      }
      
      add_diff_mat[i, j] <- add_diff
      n_BIC_update <- n_BIC_update + 1
    }
  }
  add_diff_mat <- add_diff_mat + t(add_diff_mat)
  
  re_diff_mat <- matrix(0, nrow = p, ncol = p)
  
  improve <- T
  change <- c()
  add_record <- c()
  n_mt_update <- 0 # count the number of matrix update
  
  while(improve){
    improve <- F
    
    # select the highest score
    max_select <- function(implement, mat){
      diff_max <- max(mat)
      if(diff_max != 0){
        diff_max_pos <- which(mat == diff_max, arr.ind = T)
        if(nrow(diff_max_pos) == 2){
          pos <- diff_max_pos[1,1] < diff_max_pos[2,1]
          if(pos){
            diff_max_pos <- diff_max_pos[1, ]
          }else{
            diff_max_pos <- diff_max_pos[2, ]
          }
        }else{
          diff_max_pos <- diff_max_pos[1, ]
        }
      }else{ 
        diff_max_pos <- c(0, 0)
      }
      result <- c(implement, diff_max, diff_max_pos)
      
      return(result)
    }
    
    # select the maximum and detect the cycle
    max_select_cycle_detect <- function(add_diff_mat, re_diff_mat){
      ifCycle <- T
      add_diff_cycle <- add_diff_mat
      re_diff_cycle <- re_diff_mat
      while(ifCycle){
        add_max <- max_select(1, add_diff_cycle) # implement 1, means add
        re_max <- max_select(3, re_diff_cycle) # implement 3, means reverse
        diff_max <- as.data.frame(rbind(add_max, re_max))
        names(diff_max) <- c('implement', 'diff', 'from', 'to')
        diff_max <- diff_max[which.max(diff_max$diff), ]
        
        # cycle detect
        cycle_graph <- graph
        node_from <- diff_max$from
        node_to <- diff_max$to
        if(diff_max$implement == 1){
          cycle_graph[node_from, node_to] <- 1
          ifCycle <- topological_sorting(cycle_graph)
          if(ifCycle){
            add_diff_cycle[node_from, node_to] <- 0
          }else{
            break
          }
        }else if(diff_max$implement == 3){
          cycle_graph[node_from, node_to] <- 0
          cycle_graph[node_to, node_from] <- 1
          ifCycle <- topological_sorting(cycle_graph)
          if(ifCycle){
            re_diff_cycle[node_from, node_to] <- 0
          }else{
            break
          }
        }
      }
      return(diff_max)
    }
    
    diff_max <- max_select_cycle_detect(add_diff_mat, re_diff_mat)
    
    if(diff_max$diff > 0){ # if improve
      improve = T
      
      if(diff_max$implement == 1){
        n_mt_update <- n_mt_update + 1
        
        # update graph
        add_pos <- c(diff_max$from, diff_max$to)
        add_record <- rbind(add_record, add_pos)
        node_from <- add_pos[1]
        node_to <- add_pos[2]
        graph[node_from, node_to] <- 1
        change <- rbind(change, diff_max)
        
        # update diff
        update_diff <- c()
        update_set <- seq_len(p)[-node_to]
        update_change_node <- change$from[change$to == node_to]
        for(i in update_set){
          update_graph <- graph
          if(i %in% update_change_node){
            update_graph[i, node_to] <- 0
          }else{
            update_graph[i, node_to] <- 1
          }
          if(type == 'discrete'){
            update_diff_cal <- BIC_selfwritten(node_to, update_graph, data)
            update_diff_cal <- update_diff_cal - BIC_selfwritten(node_to, graph, data)
          }else if(type == 'Gaussian'){
            update_diff_cal <- BIC_g(node_to, update_graph, data)
            update_diff_cal <- update_diff_cal - BIC_g(node_to, graph, data)
          }
          update_diff <- c(update_diff, update_diff_cal)
          n_BIC_update <- n_BIC_update + 1
        }
        update_detail <- as.data.frame(cbind(update_set, update_diff))
        
        
        # update add_diff_mat
        add_diff_mat[update_set, node_to] <- update_diff
        for(i in 1:nrow(add_record)){
          zero_pos1 <- add_record[i, 1]
          zero_pos2 <- add_record[i, 2]
          add_diff_mat[zero_pos1, zero_pos2] <- 0
        }
        if(!is.null(skeleton)){
          add_diff_mat[blacklist_diff] <- 0  # delete the ones not in skeleton
        }
        
        # update re_diff_mat
        update_re_pos <- add_record
        for(i in 1:nrow(update_re_pos)){
          re_from <- update_re_pos[i, 1]
          re_to <- update_re_pos[i, 2]
          update_graph_re1 <- graph
          update_graph_re1[re_from, re_to] <- 0
          update_graph_re2 <- update_graph_re1
          update_graph_re2[re_to, re_from] <- 1
          if(type == 'discrete'){
            re_diff_pos1 <- BIC_selfwritten(re_to, update_graph_re1, data) -
              BIC_selfwritten(re_to, graph, data)
            re_diff_pos2 <- BIC_selfwritten(re_from, update_graph_re2, data) -
              BIC_selfwritten(re_from, graph, data)
          }else if(type == 'Gaussian'){
            re_diff_pos1 <- BIC_g(re_to, update_graph_re1, data) -
              BIC_g(re_to, graph, data)
            re_diff_pos2 <- BIC_g(re_from, update_graph_re2, data) -
              BIC_g(re_from, graph, data)
          }
          
          re_diff_pos <- re_diff_pos1 + re_diff_pos2
          re_diff_mat[re_from, re_to] <- re_diff_pos
          n_BIC_update <- n_BIC_update + 2
        }
        
        
      }else if(diff_max$implement == 3){ # reverse
        # update graph
        re_pos <- c(diff_max$from, diff_max$to)
        add_record <- rbind(add_record, add_pos)
        node_to <- re_pos[1]
        node_from <- re_pos[2]
        graph[node_from, node_to] <- 1
        graph[node_to, node_from] <- 0
        change <- rbind(change, diff_max)
        
        # update diff of node to
        update_diff <- c()
        update_set <- seq_len(p)[-node_to]
        update_change_node <- change$from[change$to == node_to]
        for(i in update_set){
          update_graph <- graph
          if(i %in% update_change_node){
            update_graph[i, node_to] <- 0
          }else{
            update_graph[i, node_to] <- 1
          }
          if(type == 'discrete'){
            update_diff_cal <- BIC_selfwritten(node_to, update_graph, data)
            update_diff_cal <- update_diff_cal - BIC_selfwritten(node_to, graph, data)
          }else if(type == 'Gaussian'){
            update_diff_cal <- BIC_g(node_to, update_graph, data)
            update_diff_cal <- update_diff_cal - BIC_g(node_to, graph, data)
          }
          update_diff <- c(update_diff, update_diff_cal)
          n_BIC_update <- n_BIC_update + 1
        }
        update_detail <- as.data.frame(cbind(update_set, update_diff))
        
        # update diff of node from
        update_diff_from <- c()
        update_set_from <- seq_len(p)[-node_from]
        update_change_node <- change$from[c(change$to[-length(change$to)] == node_from, F)]
        for(i in update_set_from){
          update_graph_from <- graph
          if(i %in% update_change_node){
            update_graph_from[i, node_from] <- 0
          }else{
            update_graph_from[i, node_from] <- 1
          }
          if(type == 'discrete'){
            update_diff_from_cal <- BIC_selfwritten(node_from, 
                                                    update_graph_from, data)
            update_diff_from_cal <- update_diff_from_cal - 
              BIC_selfwritten(node_from, graph, data)
          }else if(type == 'Gaussian'){
            update_diff_from_cal <- BIC_g(node_from,
                                          update_graph_from, data)
            update_diff_from_cal <- update_diff_from_cal - 
              BIC_g(node_from, graph, data)
          }
          update_diff_from <- c(update_diff_from, update_diff_from_cal)
        }
        update_detail_from <- as.data.frame(cbind(update_set_from, update_diff_from))
        
        # update add_diff_mat about node_to
        add_diff_mat[update_set, node_to] <- update_diff
        add_diff_mat[update_set_from, node_from] <- update_diff_from
        for(i in 1:nrow(add_record)){
          zero_pos1 <- add_record[i, 1]
          zero_pos2 <- add_record[i, 2]
          add_diff_mat[zero_pos1, zero_pos2] <- 0
        }
        if(!is.null(skeleton)){
          add_diff_mat[blacklist_diff] <- 0
        }
        
        
        # update re_diff_mat
        update_re_pos <- as.matrix(change[, c(3,4)])
        re_pos_ch <- which(change$implement == 3)
        for(i in re_pos_ch){
          update_re_dir <- apply(update_re_pos, 1, function(x) paste(x, collapse = ''))
          re_dir <- change[i, c(3,4)]
          re_diff_mat[as.numeric(re_dir[1]), as.numeric(re_dir[2])] <- 0
          re_dir_rec <- paste(re_dir, collapse = '')
          update_re_pos <- update_re_pos[update_re_dir != re_dir_rec, ]
          update_re_pos <- rbind(update_re_pos, rev(as.numeric(re_dir)))
        }
        
        for(i in 1:nrow(update_re_pos)){
          re_from <- update_re_pos[i, 1]
          re_to <- update_re_pos[i, 2]
          update_graph_re1 <- graph
          update_graph_re1[re_from, re_to] <- 0
          update_graph_re2 <- update_graph_re1
          update_graph_re2[re_to, re_from] <- 1
          if(type == 'discrete'){
            re_diff_pos1 <- BIC_selfwritten(re_to, update_graph_re1, data) -
              BIC_selfwritten(re_to, graph, data)
            re_diff_pos2 <- BIC_selfwritten(re_from, update_graph_re2, data) -
              BIC_selfwritten(re_from, graph, data)
          }else if(type == 'Gaussian'){
            re_diff_pos1 <- BIC_g(re_to, update_graph_re1, data) -
              BIC_g(re_to, graph, data)
            re_diff_pos2 <- BIC_g(re_from, update_graph_re2, data) -
              BIC_g(re_from, graph, data)
          }
          re_diff_pos <- re_diff_pos1 + re_diff_pos2
          re_diff_mat[re_from, re_to] <- re_diff_pos
          n_BIC_update <- n_BIC_update + 2
        }
      }
      
    }
  }
  result <- list(graph = graph, change = change, n_mt_update = n_mt_update,
                 n_BIC_update = n_BIC_update)
  return(result)
}

