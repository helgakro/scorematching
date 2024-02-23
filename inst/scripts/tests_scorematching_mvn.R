library(MASS)
library(ggplot2)
library(matlib)

# sroot_mvn <- function(y,params){
#   mu<-params[["mu"]]
#   sigma<-params[["sigma"]]
#   n <- length(mu)
#   ep <- sigma*sqrt(pi/2)*orthopolynom::glaguerre.polynomials() #Laguerre polynomial ég veit ekki alveg hvernig þetta á að koma inn.
#
#   epp <- sigma*sqrt(2)*gamma((n+1)/2)/gamma(n/2)
#   score<-sroot_normal(y,mu,sigma)
#
#   return(score)
# }


# sroot_mvn <- function(y,params){
#   mu<-params[["mu"]]
#   sigma<-params[["sigma"]]
#   n <- length(mu)
#   x=c(-10:0.1:10) #hversu langt niður á að fara?
#   ep
#
#   epp <- sigma*sqrt(2)*gamma((n+1)/2)/gamma(n/2)
#   score<-sroot_normal(y,mu,sigma)
#
#   return(score)
# }


conditional_mvn <- function(i,x,params){
  mu <- params[["mu"]]
  precmat <- params[["prec"]]
  mu_new <- mu[i]-sum(precmat[i,-i]*(x[-i]-mu[-i]))/precmat[i,i]
  sigma_new <- 1/sqrt(precmat[i,i])
  return(list("mu"=mu_new,"sigma"=sigma_new))
}

loo_score <- function(obs, mu, precmat){
  n_obs <- nrow(obs)
  n_dim <- nrow(Sigma)
  score_vec <- matrix(0,n_obs,n_dim)
  for(i in c(1:n_dim)){
    for (k in c(1:n_obs)){
      param <- conditional_mvn(i,m[k,],list("mu"=mu,"prec"=precmat))
      score_vec[k,i] <- mean(sroot_normal(m[,i],param$mu,param$sigma))
    }
  }
  return(mean(score_vec))
}

mu <- rep(0,2)
Sigma <- matrix(c(10,3,3,2),2,2)
Sigma
precmat <- inv(Sigma)
precmat
m<-mvrnorm(n = 10000, mu,Sigma)
#m<-as.data.frame(m)
#str(m)
#ggplot(m, aes(x=V1, y=V2))+
#  geom_point(alpha = .2) +
#  geom_density_2d()+
#  theme_bw()

Sigma_indep <- matrix(c(10,0,0,2),2,2)


loo_score(m, mu, precmat)
loo_score(m,mu, inv(Sigma_indep))
loo_score(m,rep(1,2),precmat)
