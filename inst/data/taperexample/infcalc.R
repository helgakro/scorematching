###################################################################
## Calculate confidence intervals for parameters                 ##
## This uses explicit forms based on the exponential covariance  ##
###################################################################

library(fields)
library(spam)
source("functions.R")
load("anom1962.RData")
load("maxlik.RData")

## Calculate sigma2.mle and c.mle

g <- exp(-d / rho.mle)
Qg <- chol(g)
sigma2.mle <- drop(crossprod(backsolve(Qg, z, transpose = TRUE))) / n
c.mle <- sigma2.mle/rho.mle

## Calculate the information matrix for sigma^2 and rho

ginvgd <- backsolve(Qg, backsolve(Qg, g * d, transpose = TRUE))
rm(Qg); gc()
tr.ginvgd <- sum(diag(ginvgd))
tr.ginvgd2 <- sum(diag(ginvgd %*% ginvgd))
tr.iginvgd2 <- sum(diag((diag(n) + ginvgd/rho.mle) %*% (diag(n) + ginvgd/rho.mle)))
tr.iginvgdr <- sum(diag((diag(n) + ginvgd/rho.mle) %*% ginvgd))

rm(g, ginvgd)
gc()

FI <- matrix(NA, nrow = 2, ncol = 2)
FI[1, 1] <- n / (sigma2.mle^2 * 2)
FI[1, 2] <- FI[2, 1] <- tr.ginvgd / (sigma2.mle * rho.mle^2 * 2)
FI[2, 2] <- tr.ginvgd2 / (rho.mle^4 * 2)
FIinv <- solve(FI)

## Calculate the infomation matrix for sigma^2 and c

FIc <- matrix(NA, nrow = 2, ncol = 2)
FIc[1, 1] <- tr.iginvgd2 / (2 * sigma2.mle^2)
FIc[1, 2] <- FIc[2, 1] <- -  tr.iginvgdr / (2 * sigma2.mle^2)
FIc[2, 2] <- tr.ginvgd2 / (2 * sigma2.mle^2)
FIcinv <- solve(FIc)

rm(list = ls(pattern = "tr")); gc()

## Calculate what we need for the tapering estimators

setup <- make.tapersetup(d, wendland2.1, taprange = range)
g <- exp(-d / rho.mle.2taper)
rm(d); gc()

## Calculate sigma2.mle.1taper

gt <- new("spam", entries = exp(- setup$good.dists / rho.mle.1taper) * setup$taps,
          colindices = setup$ja, rowpointers = setup$ia,
          dimension = as.integer(rep(n, 2)))
Qgt <- chol(gt)
sigma2.mle.1taper <- drop(crossprod(forwardsolve(Qgt, z))) / n
c.mle.1taper <- sigma2.mle.1taper / rho.mle.1taper

rm(gt, Qgt); gc()

## Calculate sigma2.mle.2taper

gt <- new("spam", entries = exp(- setup$good.dists / rho.mle.2taper) * setup$taps,
          colindices = setup$ja, rowpointers = setup$ia,
          dimension = as.integer(rep(n, 2)))
Qgt <- chol(gt)
gtinvt <- gt
slot(gtinvt, "entries") <- backsolve(Qgt, forwardsolve(Qgt, diag(n)))[setup$index] *
  setup$taps
sigma2.mle.2taper <- drop(t(z) %*% gtinvt %*% z) / n
c.mle.2taper <- sigma2.mle.2taper / rho.mle.2taper

## Calculate the robust information matrix for sigma2 and rho

gtd <- gt
slot(gtd, "entries") <- slot(gt, "entries") * setup$good.dists

gtinvgtd <- backsolve(Qgt, forwardsolve(Qgt, gtd)) # Notice syntax is different

trace1 <- sum(diag(gtinvgtd))
trace2 <- sum(diag(gtinvgtd %*% gtinvgtd))

gtinvtg <- as.matrix(gtinvt %*% g)
gc()

trace3 <- sum(diag(gtinvtg %*% gtinvtg))
gc()

big <- new("spam", entries = (gtinvgtd %*% backsolve(Qgt, forwardsolve(Qgt, diag(n))))[setup$index] * setup$taps,
           colindices = setup$ja, rowpointers = setup$ia,
          dimension = as.integer(rep(n, 2)))
big <- as.matrix(big %*% g)
gc()

trace4 <- sum(diag(gtinvtg %*% big))

rm(gtinvtg); gc()

trace5 <- sum(diag(big %*% big))
gc()

EdG <- EGG <- matrix(NA, nrow = 2, ncol = 2)

EdG[1, 1] <- - n / (2 * sigma2.mle.2taper^2)
EdG[1, 2] <- EdG[2, 1] <- - trace1 / (2 * sigma2.mle.2taper * rho.mle.2taper^2)
EdG[2, 2] <- - trace2 / (2 * rho.mle.2taper^4)

EGG[1, 1] <- trace3 / (2 * sigma2.mle.2taper^2)
EGG[1, 2] <- EGG[2, 1] <- trace4 / (2 * sigma2.mle.2taper * rho.mle.2taper^2)
EGG[2, 2] <- trace5 / (2 * rho.mle.2taper^4)

RI <- t(EdG) %*% solve(EGG, EdG)
RIinv <- solve(RI)

## Calculate the robust information matrix for sigma2 and c

EdGc <- EGGc <- matrix(NA, nrow = 2, ncol = 2)

trace6 <- sum(diag( (diag(n) + gtinvgtd/rho.mle.2taper) %*%
                    (diag(n) + gtinvgtd/rho.mle.2taper) ))
gc()

trace7 <- sum(diag( (diag(n) + gtinvgtd/rho.mle.2taper) %*% gtinvgtd ))
gc()

big2 <- new("spam", entries = ((diag(n) + gtinvgtd / rho.mle.2taper) %*%
                               backsolve(Qgt, forwardsolve(Qgt, diag(n))))[setup$index] *
            setup$taps,
            colindices = setup$ja, rowpointers = setup$ia,
            dimension = as.integer(rep(n, 2)))
rm(Qgt, gtinvgtd); gc()
big2 <- as.matrix(big2 %*% g)
rm(g); gc()

trace9 <- sum(diag(big2 %*% big))
rm(big); gc()

trace8 <- sum(diag(big2 %*% big2))
rm(big2); gc()

EdGc[1, 1] <- - trace6 / (2 * sigma2.mle.2taper^2)
EdGc[1, 2] <- EdGc[2, 1] <- trace7 / (2 * sigma2.mle.2taper^2)
EdGc[2, 2] <- - trace2 / (2 * sigma2.mle.2taper^2)

EGGc[1, 1] <- trace8 / (2 * sigma2.mle.2taper^2)
EGGc[1, 2] <- EGGc[2, 1] <- -trace9 / (2 * sigma2.mle.2taper^2)
EGGc[2, 2] <- trace5 / (2 * sigma2.mle.2taper^2)

RIc <- t(EdGc) %*% solve(EGGc, EdGc)
RIcinv <- solve(RIc)

#load("infcalc.RData")

print(rho.mle + c(-2, 0, 2) * sqrt(FIinv[2, 2]))
print(rho.mle.2taper + c(-2, 0, 2) * sqrt(RIinv[2, 2]))
print(rho.mle.1taper)

print(sigma2.mle + c(-2, 0, 2) * sqrt(FIinv[1, 1]))
print(sigma2.mle.2taper + c(-2, 0, 2) * sqrt(RIinv[1, 1]))
print(sigma2.mle.1taper)

print(c.mle + c(-2, 0, 2) * sqrt(FIcinv[1, 1]))
print(c.mle.2taper + c(-2, 0, 2) * sqrt(RIcinv[1, 1]))
print(c.mle.1taper)

save(sigma2.mle, rho.mle, c.mle, sigma2.mle.1taper, rho.mle.1taper, c.mle.1taper,
     sigma2.mle.2taper, rho.mle.2taper, c.mle.2taper, FI, FIinv, FIc, FIcinv,
     RI, RIinv, RIc, RIcinv, file = "infcalc.RData")
