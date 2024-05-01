###################################################
## Calculate the profile likelihoods over a grid ##
###################################################

library(fields)
library(spam)
source("functions.R")
load("anom1962.RData")

## Calculate the distance matrix and save it for later use

d <- rdist.earth(loc)

## If you don't have enough memory to calculate the distance matrix
## all at once, try this:

##d <- matrix(NA, n, n)
##index <- round(seq(1, n, length = 5))
##for(i in 1:4){
##  d[index[i]:index[i+1],] <- rdist.earth(loc[index[i]:index[i+1],], loc)
##}

save(z, loc, n, d, file = "anom1962.RData")

## Calculate the profile nll's over a grid

rho.seq <- seq(20, 80, length = 30)
nll.seq <- sapply(rho.seq, nll)

## You can change the taper range (gamma in the paper) here

setup <- make.tapersetup(d, wendland2.1, taprange = 50)
rm(d); gc() # Free up some memory

nll.1taper.seq <- sapply(rho.seq, nll.1taper, setup = setup)
nll.2taper.seq <- sapply(rho.seq, nll.2taper, setup = setup)

save(rho.seq, nll.seq, nll.1taper.seq, nll.2taper.seq,
     file = "likcalc.RData")
