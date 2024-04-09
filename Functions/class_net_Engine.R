class_net_Engine <- function(X, K){
  
  kthNND_n <- nndist(X, k = K)
  vols <- numeric(npoints(X))#
  for(j in 1:npoints(X)) {#
    vols[j] <- lineardisclength(L=domain(X),x=X$data[j,,drop=TRUE],r=kthNND_n[j])#
  }#
  kthNND_n <- vols#
  em_n <- nncleanEngine2(kthNND = kthNND_n, k = K, d = 2, verbose = F) 
  pp_n <- c(em_n$delta1, em_n$delta2)
  pp_n[which(pp_n == 0)] <- 0.000000001
  entr <-  - sum(pp_n * log2(pp_n))
  zz_n <- em_n$z
  zz_n <- factor(zz_n, levels = c(0, 1))
  levels(zz_n) <- c("noise", "feature")
  clas <- as.lpp(unmark(X), marks = zz_n)
  
  if(is.multitype(X)){
    out <- list(rates = calculate_statistics(table(clas$data$marks, X$data$pattern)),
                entropy = entr, X = X, clas = clas, K = K)
  } else {
    out <- list(entropy = entr, X = X, clas = clas, K = K)
  }
  
  class(out) <- "localdetectE"
  return(out)
  
}


