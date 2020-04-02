# functions

ricker <- function(J, tmax, r.mn, r.sd, K) {
  # values: J, tmax
  # vectors: r.mn, r.sd, K (length nSpp)
  
  N <- matrix(nrow=tmax, ncol=J)
  N[1,] <- K/2
  for(k in 2:tmax) {
    r.k <- rnorm(J, r.mn*(1 - N[k-1,]/K), r.sd)
    # bug if r.mn < 0 and (1 - N/K) < 0 -- produces a positive r
    if(any(r.mn < 0 & N[k-1,]/K > 1)) r.k[r.mn < 0 & N[k-1,]/K > 1] <- -1
    N[k,] <- N[k-1,]*exp(r.k)
  }
  return(N)
}



simulate_communities <- function(J, els, pDet, beta, beta.sd, tmax) {
  # J: number of species
  # els: vector of elevations
  # pDet: dataframe with species and dbeta parameters shp1, shp2
  # beta: vector of slopes (r.mn.ij = b0.j + b1.j * el.i + b2.j * el.i^2)
  # beta.sd: vector for sd in slopes among species
  # tmax: number of years to simulate
  # 
  # returns: 
  #  N: abundances [yr, j, el]
  #  Z: occurrences [yr, j, el]
  #  rng: ranges [j, minBin/maxBin, yr]
  
  r.sd <- runif(J, 0, 1)
  K <- runif(J, 1e2, 1e4)
  
  els.std <- scale(els)
  b <- matrix(rnorm(J*length(beta), beta, beta.sd), ncol=J)
  X <- cbind(1, els.std, els.std^2)
  
  N <- purrr::map(seq_along(els), ~ricker(J, tmax, X[.,] %*% b, r.sd, K)) %>%
    unlist %>% array(dim=c(tmax, J, length(els)))
  Z <- N > 1
  
  rng <- map(1:tmax, ~t(apply(Z[.,,], 1, function(x) {
    if(any(x)) {
      range(which(x))
    } else {
      range(NA)
    }
  }))) %>% unlist %>% array(dim=c(J, 2, tmax))
  
  pDet.par <- pDet %>% select(shp1, shp2) %>% as.matrix
  
  pDet.el <- map(1:n.els, ~rbeta(J, pDet$shp1, pDet$shp2)) %>%
    do.call('rbind', .) 
  
  return(list(N=N, Z=Z, rng=rng, pDet.par=pDet.par, pDet.el=pDet.el))
}











sample_community <- function(J, els, yrs.obs, comm.true, Ytot) {
  # assign sample years to each bin randomly
  yrs.samp <- map(Ytot, ~sample(yrs.obs, size=., replace=T) %>% sort)
  
  # generate detections from species pool
  # pr(obs) = N*pDet
  Y <- matrix(0, nrow=length(els), ncol=J)
  for(i in seq_along(els)) {
    if(Ytot[i] > 0) {
      for(k in yrs.samp[[i]]) {
        spp.obs <- sample(1:J, 1, prob=comm.true$N[k,,i]*comm.true$pDet.el[i,])
        Y[i,spp.obs] <- Y[i,spp.obs] + 1
      }
    }
  }
  
  # compute ranges, interpPatchy & binsAway
  rng.obs <- matrix(0, ncol=J, nrow=2)
  binsAway <- interp.rng <- matrix(0, nrow=length(els), ncol=J)
  interpPatchy <- rep(0, J)
  for(j in 1:J) {
    # observed interpolated range
    if(any(Y[,j]>0)) {
      rng.obs[,j] <- range(which(Y[,j]>0))
      interp.rng.j <- rng.obs[1,j]:rng.obs[2,j]
      interp.rng[,j] <- 1:length(els) %in% interp.rng.j
      # within interpolated range, proportion of bins with Y=0
      interpPatchy[j] <- sum(Y[interp.rng.j,j]==0)/length(interp.rng.j)
      for(i in (1:length(els))[-interp.rng.j]) {
        # number of bins away from last detection (outside interpolated range only)
        binsAway[i,j] <- min(abs(i - rng.obs[,j])) 
      }
    }
  }
  return(list(Y=Y, spAbund=colSums(Y), rng=rng.obs, binsAway=binsAway, 
              interp.rng=interp.rng, interpPatchy=interpPatchy))
}


