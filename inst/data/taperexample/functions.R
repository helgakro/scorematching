## Taper function we're using

wendland2.1 <- function(d, taprange){
  d <- d/taprange
  return((1-d)^4*(4*d+1)*(d<1))
}

## Get the storage information we'll need to create a sparse matrix

make.tapersetup <- function(d, f.tap, taprange)
{
  n <- nrow(d)
  inrange <- d < taprange
  good.dists <- d[inrange]
  taps <- f.tap(good.dists, taprange)
  names(taps) <- names(good.dists) <- NULL
  ja <- row(d)[inrange]
  ia <- as.integer(c(1, 1 + cumsum(apply(inrange, 1, sum))))
  index <- (col(d)[inrange] - 1) * n + ja
  return(list(n = n, good.dists = good.dists, taps = taps,
              ja = ja, ia = ia, index = index))
}

## Likelihood functions -- these all have timing built into them
## For general use, you may want to remove the timings

nll <- function(x){
  print(x)
  times <- rep(NA, 5)
  times[1] <- system.time(corr.matrix <- exp(-d/x))[1]
  times[2] <- system.time(Q <- chol(corr.matrix))[1]
  times[3] <- system.time(logdet <- 2 * sum(log(diag(Q))))[1]
  times[4] <- system.time(bs <- backsolve(Q, z, transpose = TRUE))[1]
  times[5] <- system.time(distval <- drop(crossprod(bs)))[1]
  y <- (n * log(2 * pi) + n * log(distval / n) + logdet + n)/2
  attr(y, "times") <- times
  return(y)
}

## The next two functions require the spam package

nll.1taper <- function(x, setup){
  print(x)
  times <- rep(NA, 5)
  times[1] <- system.time(corr.matrix.Taper <- new("spam",
                                                   entries = exp(-setup$good.dists/x) *
                                                   setup$taps,
                                                   colindices = setup$ja,
                                                   rowpointers = setup$ia,
                                                   dimension = as.integer(rep(n, 2))))[1]
  times[2] <- system.time(Q <- chol(corr.matrix.Taper))[1]
  times[3] <- system.time(logdet <- 2 * as.numeric(determinant(Q, log = TRUE)$modulus))[1]
  times[4] <- system.time(bs <- forwardsolve(Q, z))[1] # Same as usual backsolve(, transpose=TRUE)!
  times[5] <- system.time(distval <- drop(crossprod(bs)))[1]
  y <- (n * log(2 * pi) + n * log(distval / n) + logdet + n)/2
  attr(y, "times") <- times
  return(y)
}

nll.2taper <- function(x, setup){
  print(x)
  times <- rep(NA, 6)
  times[1] <- system.time(corr.matrix.Taper <- new("spam",
                                                   entries = exp(-setup$good.dists/x) *
                                                   setup$taps,
                                                   colindices = setup$ja,
                                                   rowpointers = setup$ia,
                                                   dimension = as.integer(rep(n, 2))))[1]
  times[2] <- system.time(Q <- chol(corr.matrix.Taper))[1]
  times[3] <- system.time(logdet <- 2 * as.numeric(determinant(Q, log = TRUE)$modulus))[1]
  times[4] <- system.time(inv <- backsolve(Q, forwardsolve(Q, diag(n))))[1]
  times[5] <- system.time(slot(corr.matrix.Taper, "entries") <- inv[setup$index] *
                          setup$taps)[1]
  times[6] <- system.time(distval <- drop(t(z) %*% corr.matrix.Taper %*% z))[1]
  y <- (n * log(2 * pi) + n * log(distval / n) + logdet + n)/2
  attr(y, "times") <- times
  return(y)
}
