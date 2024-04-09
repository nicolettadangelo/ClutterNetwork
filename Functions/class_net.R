class_net <- function(X, K = NULL, n_feat = 1, n_entr = 1:35, verbose = T){
  
  if(n_feat > 3) stop("n_feat > 3 not allowed")
  
  X0 <- X
  K0 <- K
  
  Z <- as.lpp(unmark(X), marks = factor(1:npoints(X)))

  class_track <- list()
  K_track <- vector(l = n_feat)
  
  if(is.null(K0)){
    o_track <- list()
    deltas_track <- list()
  }
  
  for(n in 1:n_feat){
    
    if(verbose) cat("Classifying feature at iteration", n, "\n", "\n")
    
    if(is.null(K0)){
      
      if(verbose)  cat("Selecting K ... ")
      
      deltas <- vector()
      
      for(i in n_entr){
        deltas[i] <- class_net_Engine(unmark(Z), K = i)$entropy
      }
      
      x0 <- n_entr
      z0 <- - x0
      out.lm <- lm(deltas ~ 1)
      o <- try(segmented(out.lm, ~ z0), silent = T) 
      K <- - round(o$psi[2])
      
      o_track[[n]] <- o
      deltas_track[[n]] <- deltas
      
      if(verbose) cat(paste("Estimated K =", K, "\n", "\n"))
    }
    
    kthNND_n <- nndist(Z, k = K)
    vols <- numeric(npoints(X))#
    for(j in 1:npoints(X)) {#
      vols[j] <- lineardisclength(L=Z$domain,x=X$data[j,,drop=TRUE],r=kthNND_n[j])#
    }#
    kthNND_n <- vols#
    em_n <- nncleanEngine2(kthNND = kthNND_n, k = K, d = 2, verbose = F) 
    pp_n <- c(em_n$delta1, em_n$delta2)
    entr <-  - sum(pp_n * log2(pp_n))
    zz_n <- em_n$z
    zz_n <- factor(zz_n, levels = c(0, 1))
    levels(zz_n) <- c("noise", "feature")
    
    clas <- lpp(cbind(Z$data, zz_n), L = Z$domain)

    class_track[[n]] <- clas
    K_track[n] <- K
      
      Z <- clas[clas$data[, n + 5, drop = T] == "feature"]
      
      if(n == n_feat){
        id <- as.numeric(Z$data$marks)
        v <- vector(l = npoints(X0))
        v[id] <- "feature"
        v[ - id] <- "clutter"
      } 
      
  }
  
  if(verbose) cat(paste("Procedure ended. \n"))
  
  if(is.null(K0)){
        out <- list(n_feat = n_feat, K_track = K_track, K0 = K0,
                    X = X0, id = v, o_track = o_track,
                    class_track = class_track, n_entr = n_entr,
                    deltas_track = deltas_track)
    
  } else{
    out <- list(n_feat = n_feat, K_track = K_track, K0 = K0, 
                X = X0, id = v,
                 class_track = class_track, n_entr = n_entr)
  }
  class(out) <- "localdetect"
  return(out)
  
}





