######################################
## Maximize the profile likelihoods ##
######################################

library(fields)
library(spam)
source("functions.R")
load("anom1962.RData")
load("likcalc.RData")

## Minimize the profile negative log-likelihood
## Take range from minimizing value over the grid

interval <- rho.seq[(1:length(rho.seq))[nll.seq == min(nll.seq)] + c(-1, 1)]
rho.mle <- optimize(f = nll, interval = interval)$minimum

setup <- make.tapersetup(d, wendland2.1, taprange = 50)
rm(d); gc()

## Minimize the 1taper profile negative log-likelihood

interval <- rho.seq[(1:length(rho.seq))[nll.1taper.seq == min(nll.1taper.seq)] + c(-1, 1)]
rho.mle.1taper <- optimize(f = nll.1taper, interval = interval, setup = setup)$minimum

## Minimize the 2taper profile negative log-likelihood

interval <- rho.seq[(1:length(rho.seq))[nll.2taper.seq == min(nll.2taper.seq)] + c(-1, 1)]
rho.mle.2taper <- optimize(f = nll.2taper, interval = interval, setup = setup)$minimum

save(rho.mle, rho.mle.1taper, rho.mle.2taper, file = "maxlik.RData")
