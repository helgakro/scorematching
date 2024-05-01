#########################################
## Timing for one likelihood evalation ##
#########################################

library(fields)
library(spam)
library(xtable) # If you want to put the results in a latex table
source("functions.R")
load("anom1962.RData")

k <- 10 # Number of repetitions to use for timing

t0 <- t1 <- matrix(NA, nrow = k, ncol = 5)
t2 <- matrix(NA, nrow = k, ncol = 6)

## Timing for full likelihood

for(i in 1:k) t0[i,] <- attr(nll(50), "times") # Pick an arbitrary value to evaluate
t0 <- apply(t0, 2, mean)

## Timing for initialization of taper matrix

setuptime <- system.time(setup <- make.tapersetup(d, wendland2.1, taprange = 50))
rm(d); gc()

## Timing for tapered versions

for(i in 1:k) t1[i,] <- attr(nll.1taper(50, setup = setup), "times")
t1 <- apply(t1, 2, mean)

for(i in 1:k) t2[i,] <- attr(nll.2taper(50, setup = setup), "times")
t2 <- apply(t2, 2, mean)

times <- cbind(c(t0[1:4], rep(NA, 2), t0[5]),
               c(t1[1:4], rep(NA, 2), t1[5]),
               c(t2[1:3], NA, t2[4:6]))

## Tapered Cholesky -- should be averaged, since same calculation

times[2, 2] <- times[2, 3] <- mean(times[2,2:3])
times <- rbind(times, apply(times, 2, sum, na.rm = TRUE))

rownames(times) <- c("Gamma or Gamma circ T",
                     "Cholesky decomposition",
                     "Log determinant",
                     "Backsolve",
                     "Full solve",
                     "Second taper",
                     "Quadratic form",
                     "Total")
colnames(times) <- c("No taper", "One taper", "Two tapers")
xtable(times)

save(times, setuptime, file = "timing.RData")
