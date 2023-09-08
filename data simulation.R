## self-written
data_simulation <- function(n, p, level, n_neighbor, type){
  # n: sample size
  # p: number of variable
  # level: vector, numbers of levels that each variable allows to take
  # n_neighbor: neighbor size
  # type: 'discrete' or 'Gaussian'
  
  
  # generate skeleton (PDAG)
  regenerate <- T
  while(regenerate){
    skeleton <- matrix(rbinom(p*p, 1, (n_neighbor / (p-1))), nrow = p)
    # consider the upper triangle
    skeleton[lower.tri(skeleton)] <- 0
    diag(skeleton) <- 0
    
    # if there are no roots, regenerate the skeleton
    ifRegenerate <- apply(as.matrix(skeleton), 2, sum)
    ifRegenerate <- ifRegenerate == 0
    if(sum(ifRegenerate) != 0){
      regenerate = F
    }
  }
  
  
  if(type == 'discrete'){
    # generate probabilities
    prob <- list()
    for(i in 1:ncol(skeleton)){  # i is the node
      i_lv_num <- level[i]
      i_pr_num <- sum(skeleton[, i])  # the number of parents of i
      i_pr <- which(skeleton[, i] == 1)
      pr_lv_num <- level[i_pr]
      
      lv_num <- c(i_lv_num, pr_lv_num)
      prod_num <- tail(cumprod(lv_num), 1)
      prob_i <- array(runif(prod_num), dim = lv_num)
      prob <- c(prob, list(prob_i))
    }
    
    # generate data
    generate <- T
    
    while(generate){
      data <- matrix(0, nrow = n, ncol = ncol(skeleton))
      
      for(j in 1:n){
        ancestor <- which(apply(skeleton, 2, sum) == 0)
        pr_lv <- c()
        for(t in ancestor){
          lv <- sample(level[t], 1, prob = prob[[t]])
          pr_lv <- c(pr_lv, lv)
        }
        data[j, ancestor] <- pr_lv
        pr_num <- ancestor
        
        parent <- data.frame(pr_num, pr_lv)
        node_for_gen <- parent  # used for data generation
        node_for_gen$lv <- level[node_for_gen$pr_num]
        
        while(nrow(parent) != 0){
          
          
          ch_num <- which(skeleton[parent$pr_num, ] == 1, arr.ind = T)
          if(is.matrix(ch_num)){   # if there are more than 2 parents, then children need to be unique()
            ch_num <- unique(ch_num[, 2])
          }
          # delete the ch that did generate before
          ch_num <- ch_num[! ch_num %in% node_for_gen$pr_num]
          ch_lv <- c()
          for(i in ch_num){
            ch_prob <- prob[[i]]
            parent_direct <- which(skeleton[, i] == 1)
            
            # turn the probability array into list
            ch_prob <- apply(ch_prob, 1, c)
            # calculate the position
            pos <- node_for_gen[node_for_gen$pr_num %in% parent_direct, ]
            # pos <- parent[parent$pr_num == parent_direct, 2]
            if(nrow(pos) == 1){
              pos <- pos$pr_lv
            }else{
              pos1 <- pos$pr_lv[1]
              for(t in 2:nrow(pos)){
                pos1 <- pos1 + (pos$pr_lv[t] - 1) * prod(pos$lv[1:t-1])
                
              }
              pos <- pos1
            }
            ch_prob <- ch_prob[pos, ]
            ch_sample <- sample(level[i], 1, prob = ch_prob)
            
            data[j, i] <- ch_sample
            ch_lv <- c(ch_lv, ch_sample)
            node_for_gen <- rbind(node_for_gen, c(i, ch_sample, level[i]))
            node_for_gen <- node_for_gen[order(node_for_gen$pr_num),]
          }
          
          parent <- data.frame(pr_num = ch_num, pr_lv = ch_lv)
          
          # if one path reaches the end, then delete the node from parent directly
          posi_child <- apply(matrix(skeleton[parent$pr_num, ], ncol = ncol(skeleton)), 1, sum)
          end_node <- which(posi_child == 0)
          if(length(end_node) != 0){
            parent <- parent[-end_node, ]
          }
        }
      }
      
      # all the data column should have all the levels
      generate <- F
      ifEveryLevel <- apply(as.matrix(data), 2, unique)
      ifEveryLevel <- length(unlist(ifEveryLevel))
      if(ifEveryLevel < sum(level)){
        generate <- T
      }
    }
  }
  
  
  if(type == 'Gaussian'){
    # generate the weight matrix
    wgt_pos <- which(skeleton == 1, arr.ind = T)
    wgt <- skeleton
    for(i in 1:nrow(wgt_pos)){
      wgt_row <- wgt_pos[i, 1]
      wgt_col <- wgt_pos[i, 2]
      wgt[wgt_row, wgt_col] <- runif(1, min = .1, max = 1)
    }
    wgt <- t(wgt)
    
    # generate multinormal data
    Sigma_para <- diag(1, nrow = p)
    Sigma <- solve(diag(1, nrow = p) - wgt) %*% solve(t(diag(1, nrow = p) - wgt))
    mu <- rep(0, p)
    data <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
  }
  
  PDAG <- skeleton
  output <- list(PDAG = PDAG, data = data)
  return(output)
}





## simulate discrete and continuous data respectively
## discrete one uses the package 'bnlearn'
## quicker
library(gtools)
library(bnlearn)


# discrete data simulation
data_discrete_simulation <- function(N, n_set, p_set, n_neighbor,
                                     type = 'discrete'){
  # N: repeat times
  # n_set: vector, list of sample sizes
  # p_set: vector, list of variable sizes
  
  dg.dir <- 'E:data generation\\'
  
  for(n in n_set){
    for(p in p_set){
      level <- rep(2, p)  # binomial level=2
      name.dir <- paste('n=',n,',','p=',p,sep = '')
      name.dir <- paste(dg.dir, name.dir, sep = '')
      setwd(name.dir)
      
      cat('n=',n,'p=',p)
      
      for(n_sample in 1:N){
        # generate PDAG
        regenerate <- T
        while(regenerate){
          PDAG <- matrix(rbinom(p*p, 1, (n_neighbor / (p-1))), nrow = p)
          # consider the upper triangle
          PDAG[lower.tri(PDAG)] <- 0
          diag(PDAG) <- 0
          
          # if there are no roots, regenerate the skeleton
          ifRegenerate <- apply(as.matrix(PDAG), 2, sum)
          ifRegenerate <- ifRegenerate == 0
          if(sum(ifRegenerate) != 0){
            regenerate = F
          }
        }
        
        # generate model
        PDAG_string <- c()
        for(i in 1:ncol(PDAG)){
          parent <- which(PDAG[, i] != 0)
          node_string <- as.character(i)
          if(length(parent) > 0){
            dist_string <- paste(parent, collapse = ':')
            node_string <- paste(i, dist_string, sep = '|')
          }
          node_string <- paste('[', node_string, ']', sep = '')
          PDAG_string <- c(PDAG_string, node_string)
        }
        PDAG_string <- paste0(PDAG_string, collapse = '')
        bn <- model2network(PDAG_string)
        
        # generate prob & data
        prob <- list()
        for(i in 1:ncol(PDAG)){  # i is the node
          i_lv_num <- level[i]
          i_pr_num <- sum(PDAG[, i])  # the number of parents of i
          i_pr <- which(PDAG[, i] == 1)
          pr_lv_num <- level[i_pr]
          
          lv_num <- c(i_lv_num, pr_lv_num)
          prod_num <- tail(cumprod(lv_num), 1)
          
          # set all prob > .1, avoid near-unfaithfulness
          prob_i_pre <- rdirichlet(prod_num/2, alpha = rep(1, i_lv_num))
          prob_i <- c()
          for(j in 1:nrow(prob_i_pre)){
            while(any(prob_i_pre[j, ] < .1)){
              prob_i_pre[j, ] <- rdirichlet(1, alpha = rep(1, i_lv_num))
            }
            prob_i <- c(prob_i, prob_i_pre[j, ])
          }
          
          prob_i <- array(prob_i, dim = lv_num)
          prob <- c(prob, list(prob_i))
        }
        names(prob) <- as.character(c(1:p))
        fit <- custom.fit(bn, dist = prob)
        data <- rbn(fit, n)
        data <- as.vector(unlist(lapply(data, as.numeric)))
        data <- matrix(data, nrow = n, ncol = p, byrow = F)
        
        data_sim <- list(PDAG = PDAG, data = data)
        name <- paste('data_sim',n_sample,'.RData', sep = '')
        
        save(data_sim, file = name)
      }
    }
  }
}






# continous data simulation
library(MASS)
data_continuous_simulation <- function(N, n_set, p_set, n_neighbor,
                                       type = 'Gaussian'){
  save_dir <- 'E:\\data generation continuous\\'
  for(n in n_set){
    for(p in p_set){
      name.dir <- paste('n=',n,',','p=',p,sep = '')
      name.dir <- paste(save_dir, name.dir, sep = '')
      setwd(name.dir)
      
      cat('n=',n,'p=',p)
      
      for(j in 1:N){
        # generate PDAG
        regenerate <- T
        while(regenerate){
          PDAG <- matrix(rbinom(p*p, 1, (n_neighbor / (p-1))), nrow = p)
          # consider the upper triangle
          PDAG[lower.tri(PDAG)] <- 0
          diag(PDAG) <- 0
          
          # if there are no roots, regenerate the skeleton
          ifRegenerate <- apply(as.matrix(PDAG), 2, sum)
          ifRegenerate <- ifRegenerate == 0
          if(sum(ifRegenerate) != 0){
            regenerate = F
          }
        }
        
        if(type == 'Gaussian'){
          # generate the weight matrix
          wgt_pos <- which(PDAG == 1, arr.ind = T)
          wgt <- PDAG
          for(i in 1:nrow(wgt_pos)){
            wgt_row <- wgt_pos[i, 1]
            wgt_col <- wgt_pos[i, 2]
            wgt[wgt_row, wgt_col] <- runif(1, min = .1, max = 1)
          }
          wgt <- t(wgt)
          
          # generate multinormal data
          Sigma_para <- diag(1, nrow = p)
          Sigma <- solve(diag(1, nrow = p) - wgt) %*% solve(t(diag(1, nrow = p) - wgt))
          mu <- rep(0, p)
          data <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
        }
        
        data_sim <- list(data = data, PDAG = PDAG)
        file_name <- paste('data_sim',j, '.RData',sep='')
        save(data_sim, file = file_name)
      }
    }
  }
}



# an example of the parameters
n_set <- c(50, 250, 500, 1000)
p_set <- c(10, 50, 250, 500)
n_neighbor <- 4
N = 500
