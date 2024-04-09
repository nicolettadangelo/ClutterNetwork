dknn2 <- function (x, k = 1, d = 2, lambda = 1) {
  validposint(k, "dknn")
  validposint(d, "dknn")
  # alpha.d <- (2 * pi ^ (d / 2))/(d * gamma(d / 2))
  alpha.d <- (2 * 2 ^ (d / 2))/(d * gamma(d / 2))
  # y <- dgamma(x ^ d, shape = k, rate = lambda * alpha.d)
  y <- dgamma(x, shape = k, rate = lambda * alpha.d)
  # y <- y * d * x ^ (d - 1)
  return(y)
}
