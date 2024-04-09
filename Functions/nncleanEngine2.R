nncleanEngine2 <-
  function(kthNND, k, d,
           tol = 0.001, maxit = 50,
           plothist = FALSE, lineargs = list(), 
           verbose = TRUE, Xname = "X"){
    ## Adapted from statlib file NNclean.q
    ## Authors: Simon Byers and Adrian Raftery
    ## Adapted for spatstat by Adrian Baddeley
    
    n <- length(kthNND)
    
    ## Error handler by Adrian
    if(k >= n) {
      if(verbose)
        cat(paste("Cannot compute neighbours of order k =", k,
                  "for a pattern of", n, "data points;",
                  "treating all points as noise"),
            call. = FALSE)
      return(list(z = rep(0, n),
                  probs = rep(0, n),
                  lambda1 = NA, lambda2 = NA, p = 0,
                  kthNND = kthNND, d = d, n = n, k = k,
                  niter = 0, maxit = maxit,
                  converged = TRUE,
                  hist = NULL))
    }
    
    ## Undocumented extension by Adrian Baddeley 2014
    ## Allow different dimensions in feature and noise.
    ## d[1] is cluster dimension.
    
    d <- ensure2vector(d)
    # alpha.d <- (2. * pi ^ (d / 2.)) / (d * gamma(d / 2.))
    alpha.d <- (2. * 2 ^ (d / 2.)) / (d * gamma(d / 2.))
    
    # raise to power d for efficiency
    # kNNDpowd1 <- kthNND ^ (d[1])
    # kNNDpowd2 <- kthNND ^ (d[2])
    kNNDpowd1 <- kthNND
    kNNDpowd2 <- kthNND
    
    #
    # Now use kthNND in E-M algorithm
    # First set up starting guesses.
    #
    #
    probs <- numeric(n)
    thresh <- (min(kthNND) + diff(range(kthNND)) / 3.)
    high <- (kthNND > thresh)
    delta <- as.integer(high)
    p <- 0.5
    lambda1 <- k / (alpha.d[1] * mean(kNNDpowd1[!high]))
    lambda2 <- k / (alpha.d[2] * mean(kNNDpowd2[high]))
    loglik.old <- 0.
    loglik.new <- 1.
    #
    # Iterator starts here, 
    #
    Z <- !kthNND
    niter <- 0
    while(abs(loglik.new - loglik.old)/(1 + abs(loglik.new)) > tol) {
      if(niter >= maxit) {
        warning(paste("E-M algorithm failed to converge in",
                      maxit, ngettext(maxit, "iteration", "iterations")),
                call. = FALSE)
        break
      }
      niter <- niter + 1
      # E - step
      # f1 <- dknn(kthNND[!Z], lambda = lambda1, k = k, d = d[1])
      # f2 <- dknn(kthNND[!Z], lambda = lambda2, k = k, d = d[2])
      f1 <- dknn2(kthNND[!Z], lambda = lambda1, k = k, d = d[1])
      f2 <- dknn2(kthNND[!Z], lambda = lambda2, k = k, d = d[2])
      delta[!Z] <- (p * f1) / (p * f1 + (1 - p) * f2)
      delta[Z] <- 0
      # M - step
      sumdelta <- sum(delta)
      negdelta <- 1. - delta
      p <- sumdelta / n
      lambda1 <- (k * sumdelta) / (alpha.d[1] * sum(kNNDpowd1 * delta))
      lambda2 <- (k * (n - sumdelta)) / (alpha.d[2] * sum(kNNDpowd2 * negdelta))
      # evaluate marginal loglikelihood
      loglik.old <- loglik.new
      loglik.new <- sum( - p * lambda1 * alpha.d[1] * (kNNDpowd1 * delta)
                         - (1. - p) * lambda2 * alpha.d[2] * (kNNDpowd2 * negdelta)
                         + delta * k * log(lambda1 * alpha.d[1]) +
                           negdelta * k * log(lambda2 * alpha.d[2]))
      if(verbose) 
        cat(paste("Iteration", niter, "\tlogLik =", loglik.new,
                  "\tp =", signif(p,4), "\n"))
    }
    #
    # delta1 <- dknn(kthNND[!Z], lambda = lambda1, k = k, d = d[1])
    # delta2 <- dknn(kthNND[!Z], lambda = lambda2, k = k, d = d[2])
    delta1 <- dknn2(kthNND[!Z], lambda = lambda1, k = k, d = d[1])
    delta2 <- dknn2(kthNND[!Z], lambda = lambda2, k = k, d = d[2])
    probs[!Z] <- delta1 / (delta1 + delta2)
    probs[Z] <- 1
    #
    if(verbose) {
      cat("Estimated parameters:\n")
      cat(paste("p [cluster] =", signif(p, 5), "\n"))
      cat(paste("lambda [cluster] =", signif(lambda1, 5), "\n"))
      cat(paste("lambda [noise]   =", signif(lambda2, 5), "\n"))
    }
    #
    # z will be the classifications. 1= in cluster. 0= in noise. 
    #
    return(list(delta1 = delta1,
                delta2 = delta2, 
                z = round(probs),
                probs = probs,
                lambda1 = lambda1, lambda2 = lambda2, p = p,
                kthNND = kthNND, d = d, n = n, k = k,
                niter = niter, maxit = maxit,
                converged = (niter >= maxit),
                hist = if(plothist) H else NULL))
  }



