#########################################################################################
# Analysis on hierarchy structure of complex networks (Corominas-Murtra et al. 2013 PNAS)
# Written by Masato Abe, 2020/7/9
#########################################################################################

###############
### Library ###
###############
library(igraph)    

################
### Function ###
################

# transform to Gc
GET_GC_FUNC <- function(g){
  scc <- components(g, mode = c("strong"))
  gmat2 <- matrix(0, ncol = scc$n, nrow = scc$n)
  if(scc$no < vcount(g)){
    elist <- get.edgelist(g)
    for(i in 1:nrow(elist)){
      gmat2[scc$membership[elist[i, 1]], scc$membership[elist[i, 2]]] <- 1
    }
    diag(gmat2) <- 0
    v <- graph.adjacency(gmat2)
  }else{
    v <- g
  }
  return(v)
}


entropym <- function(g){
  ent <- 0
  indeg <- degree(g, mode = c("in"))
  outdeg <- degree(g, mode = c("out"))
  indeg0.num <- which(indeg == 0)
  outdeg0.num <- which(outdeg == 0)

  n <- vcount(g)
  m <- length(indeg0.num)
  mu <- length(outdeg0.num)
  adj <- get.adjacency(g)
  B <- matrix(ncol = n - m, nrow = n - m)
  hloc <- indeg[indeg > 0]
  L <- 0
  matsum <- 1
  seki <- adj
  while(matsum > 0){
    seki <- seki %*% adj
    matsum <- sum(seki)
    L <- L + 1
  }
  
  k <- 1
  l <- 1
  for(i in 1:n){
    l <- 1
    for(j in 1:n){
      hantei <- 0
      for(h in 1:length(indeg0.num)){
        if(i == indeg0.num[h]){
          hantei <- 1
        }
        if(j == indeg0.num[h]){
          hantei <- 1
        }
      }
      if(hantei == 0){
        B[k,l] <- adj[i,j]/sum(adj[,j])
        l <- l + 1
      }
      
    }
    if(l>1){
      k <- k + 1
    }
  }
  
  phi <- t(B) + diag(n - m)
  
  for(i in 2:L){
    tB <- t(B)
    for(j in 1:(i - 1)){
      tB <- tB %*% (t(B))
    }
    phi <- phi + tB
  }	
  
  ent <- 0	
  for(i in 1:(length(which(outdeg[-indeg0.num] == 0)))){
    for(j in 1:(n - m)){
      ent <- ent + phi[which(outdeg[-indeg0.num] == 0)[i], j] * log(hloc[j])	
    }
  }
  ent <- ent/mu
  return(ent)
}


entropyM <- function(g){
  
  ent <- 0
  indeg <- degree(g,mode = c("in"))
  outdeg <- degree(g,mode = c("out"))
  indeg0.num <- which(indeg == 0)
  outdeg0.num <- which(outdeg == 0)
  n <- vcount(g)
  m <- length(indeg0.num)
  mu <- length(outdeg0.num)
  adj <- get.adjacency(g)
  B <- matrix(ncol = n - mu, nrow = n - mu)
  hloc <- outdeg[outdeg>0]
  L <- 0
  matsum <- 1
  seki <- adj
  while(matsum > 0){
    seki <- seki %*% adj
    matsum <- sum(seki)
    L <- L + 1
  }
  
  k <- 1
  l <- 1
  for(i in 1:n){
    l <- 1
    for(j in 1:n){
      hantei <- 0
      for(h in 1:length(outdeg0.num)){
        if(i == outdeg0.num[h]){
          hantei <- 1
        }
        if(j == outdeg0.num[h]){
          hantei <- 1
        }
      }
      if(hantei == 0){
        B[k, l] <- adj[i, j]/sum(adj[i, ])
        l <- l + 1
      }
      
    }
    if(l > 1){
      k <- k + 1
    }
  }
  
  phi <- B + diag(n - mu)
  
  for(i in 2:L){
    tB <- B
    for(j in 1:(i-1)){
      tB <- tB %*% (B)
    }
    phi <- phi + tB
  }	
  
  ent <- 0	
  for(i in 1:(length(which(indeg[-outdeg0.num] == 0)))){
    for(j in 1:(n - mu)){
      ent <- ent + phi[which(indeg[-outdeg0.num]==0)[i],j] * log(hloc[j])	
    }
  }
  ent <- ent/m
  return(ent)
}

fg <- function(g){
  a <- 0
  composize <- clusters(g)$no
  
  if(composize == 1){
    {if(max(entropyM(g), entropym(g)) == 0){
      a <- 0
    }
      else{
        a <- (entropyM(g) - entropym(g))/max(entropyM(g), entropym(g))
      }
    }
    
    {if(abs(a) < 0.00000001){
      return(0)
    }
      else {
        return(a)
      }}
  }
  
  else{
    numsing <- 0
    for(i in 1:composize){
      mem <- which(clusters(g)$membership == i)
      if(length(mem) > 1){
        subg <- get.adjacency(g)[mem, mem]
        subg <- graph.adjacency(subg)
        {if(max(entropyM(subg), entropym(subg)) == 0){
          b <- 0
        }
          else{
            b <- (entropyM(subg) - entropym(subg))/max(entropyM(subg), entropym(subg))
          }
        }
      }
      
      else{
        b <- 0
        numsing <- numsing + 1
      }
      a <- a + b * length(mem)
    }
    a <- a/(vcount(g))
    return(a)
  }
}

### Treeness
T_FUNC <- function(g){
  newg <- GET_GC_FUNC(g)
  
  if(vcount(newg) > 1){
  v <- fg(newg)
  
  gd <- newg
  gu <- newg
  hantei <- 0
  w <- 1
  
  while(hantei == 0){
    gd <- delete.vertices(gd, (which(degree(gd, mode=c("in")) == 0)))
    if(sum(degree(gd, mode=c("out"))) == 0){
      hantei <- 1
    }
    w <- w + 1
  }
  
  gd <- newg
  gu <- newg
  if(w > 2){
    for(i in 1:(w - 2)){
      gd <- delete.vertices(gd,(which(degree(gd, mode=c("in")) == 0)))
      gu <- delete.vertices(gu,(which(degree(gu, mode=c("out")) == 0)))
      v <- v + fg(gd) + fg(gu)
    }
  }

  v <- v/(2*w - 3)
  }else{
    v <- 0
  }
  return(v)
}


### Orderability
O_FUNC <- function(g){
  scc <- components(g, mode = c("strong"))
  v <- sum(scc$csize == 1) / vcount(g)
  return(v)
}


CALC_G <- function(newg, w){
  id.M <- which(degree(newg, mode = "in") == 0)
  id.m <- which(degree(newg, mode = "out") == 0)
  v <- c(0, 0)
  for(i in 1:length(id.M)){
    for(j in 1:length(id.m)){
      p <- all_simple_paths(newg, from = id.M[i], to = id.m[j])
      if(length(p) > 0){
      for(k in 1:length(p)){
      v[1] <- v[1] + length(p[[k]]) / sum(w[p[[k]]])
      v[2] <- v[2] + 1
      }
      }
    }
  }
  return(v)
}

### Feedforwardness
F_FUNC <- function(g){
  scc <- components(g, mode = c("strong"))
  if(scc$no > 1){
    gc <- c(0, 0)
    if(scc$no < vcount(g)){
      newg <- GET_GC_FUNC(g)
      w <- scc$csize
      
      while(ecount(newg) > 0){
      gc <- gc + CALC_G(newg, w)
      id.m <- which(degree(newg, mode = c("out")) == 0)  # bottom up leaf removal
      newg <- delete.vertices(newg, id.m)
      w <- w[-id.m]
      }
      
      v <- gc[1] / gc[2]
      
    }else{
      v <- 1
    }
  }else{
    v <- 0
  }
  return(v)
}




############
### MAIN ###
############