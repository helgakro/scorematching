library(ggplot2)
library(matlib)
library(mvtnorm)
library(tictoc)


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


conditional_mvn_old <- function(i,x,params){
  mu <- params$mu
  precmat <- params$prec
  mu_new <- mu[i]-sum(precmat[i,-i]*(x[-i]-mu[-i]))/precmat[i,i]
  sigma_new <- 1/sqrt(precmat[i,i])
  return(list("mu"=mu_new,"sigma"=sigma_new))
}

loo_score_old <- function(obs, mu, precmat){
  n_obs <- nrow(obs)
  n_dim <- nrow(precmat)
  score_vec <- matrix(0,n_obs,n_dim)
  for(i in c(1:n_dim)){
    for (k in c(1:n_obs)){
      param <- conditional_mvn_old(i,obs[k,],list("mu"=mu,"prec"=precmat))
      score_vec[k,i] <- mean(sroot_normal(obs[k,i],param$mu,param$sigma))
    }
  }
  return(mean(score_vec))
}


conditional_mvn <- function(i,x,params){
  mu <- params$mu
  precmat <- params$prec
  tmp<-t(precmat[i,-i]*t(x[,-i]-mu[-i]))
  if(is.matrix(tmp)){
    mu_new <- mu[i]-rowSums(tmp)/precmat[i,i]
  }else{
    mu_new <- mu[i]-tmp/precmat[i,i]
  }
  sigma_new <- 1/sqrt(precmat[i,i])
  return(list("mu"=mu_new,"sigma"=sigma_new))
}

loo_score <- function(obs, mu, precmat){
  n_obs <- nrow(obs)
  n_dim <- nrow(precmat)
  score_vec <- nrow(n_dim)
  for(i in c(1:n_dim)){
      param <- conditional_mvn(i,obs,list("mu"=mu,"prec"=precmat))
      score_vec[i] <- mean(sroot_normal(obs[,i],param$mu,param$sigma))
  }
  return(mean(score_vec))
}

loo_score_2 <- function(obs, mu, precmat){
  n_obs <- nrow(obs)
  obs <- t(obs)
  n_dim <- nrow(precmat)
  score_vec <- nrow(n_dim)
  Qmu <- precmat%*%mu #matrix multiplication
  Qx <- precmat%*%obs
  for(i in c(1:n_dim)){
    score_vec[i] <- mean(sroot_normal(obs[i,],(Qmu[i]-Qx[i,])/precmat[i,i]+obs[i,],1/sqrt(precmat[i,i])))
  }
  return(mean(score_vec))
}




#####################################
# Small example comparison
#####################################

mu <- rep(0,2)
Sigma <- matrix(c(10,3,3,2),2,2)
Sigma
precmat <- inv(Sigma)
precmat
m<-rmvnorm(n = 1000, mu,Sigma)
#m<-as.data.frame(m)
#str(m)
#ggplot(m, aes(x=V1, y=V2))+
#  geom_point(alpha = .2) +
#  geom_density_2d()+
#  theme_bw()

Sigma_indep <- matrix(c(10,0,0,2),2,2)


loo_score(m, mu, precmat)
loo_score(m,mu, inv(Sigma_indep))
loo_score(m,rep(2,2),precmat)





########################################
# Parameter estimation
########################################

get_prec_mat <- function(rho,n){
  Qmat <- diag(x=1+rho^2,n)
  Qmat[1,1]=1
  Qmat[n,n]=1
  for(i in c(1:n-1)){
    Qmat[i,i+1]=-rho
    Qmat[i+1,i]=-rho
  }
  return(Qmat)
}

get_cov_mat <- function(rho,n){
  Sigma <- diag(x=1,n)
  for(i in c(1:n-1)){
    for (j in c(1:n-1)){
      if(i+j<=n){
        Sigma[i,i+j]=rho^j
        Sigma[i+j,i]=rho^j
      }
    }
  }
  return(Sigma/(1-rho^2))
}


narr <- 2^c(1:10)
times_score <- rep(0,length(narr))
times_log <- rep(0,length(narr))
o_score_list <- vector("list", length(narr))
o_log_list <- vector("list", length(narr))
times_score_2 <- rep(0,length(narr))
o_score_list_2 <- vector("list", length(narr))


for(i_n in c(1:length(narr))){

print(paste("Iteration",i_n))
rhotrue <- 0.5
n <- narr[i_n]
mu<- rep(0,n)
Qmat <- get_prec_mat(rhotrue,n)
m<-rmvnorm(n = 1000, mu,inv(Qmat))


my_obj_func_old <- function(par){
  rho <- par[1]
  #mu <- par[-1]
  mu <- rep(par[2],n)
  return(loo_score_old(m,mu,get_prec_mat(rho,ncol(m))))
}

my_obj_func <- function(par){
  rho <- par[1]
  #mu <- par[-1]
  mu <- rep(par[2],n)
  return(loo_score(m,mu,get_prec_mat(rho,ncol(m))))
}

my_obj_func_2 <- function(par){
  rho <- par[1]
  #mu <- par[-1]
  mu <- rep(par[2],n)
  return(loo_score_2(m,mu,get_prec_mat(rho,ncol(m))))
}

# my_log_obj_func <- function(par){
#   rho <- par[1]
#   mu <- par[-1]
#   return(-mean(log(dmvnorm(m,mu,get_cov_mat(rho,ncol(m))))))
# }

my_log_obj_func <- function(par){
  rho <- par[1]
  #mu <- par[-1]
  mu <- rep(par[2],n)
  return(-mean(log(dmvnorm(m,mu,get_cov_mat(rho,ncol(m))))))
  #return(-mean(log(dmvnorm(m,mu,inv(get_prec_mat(rho,ncol(m)))))))
}

rho0<-0
#mu0 <- rep(1,n)
mu0<-1
starttime <- Sys.time()
o2<-optim(par=c(rho0,mu0),my_log_obj_func,control=list(maxit=50000))
endtime <- Sys.time()
times_log[i_n]<-difftime(endtime,starttime, units="secs")
o_log_list[[i_n]]<-o2

starttime <- Sys.time()
o1<-optim(par=c(rho0,mu0),my_obj_func,control=list(maxit=50000))
endtime <- Sys.time()
times_score[i_n]<-difftime(endtime,starttime, units="secs")
o_score_list[[i_n]]<-o1
# tic("Score optim old")
# optim(par=c(rho0,mu0),my_obj_func_old,control=list(maxit=2000))
# toc()

starttime <- Sys.time()
o1_2<-optim(par=c(rho0,mu0),my_obj_func_2,control=list(maxit=50000))
endtime <- Sys.time()
times_score_2[i_n]<-difftime(endtime,starttime, units="secs")
o_score_list_2[[i_n]]<-o1_2




print(paste("Times: score",times_score[i_n],"score_2",times_score_2[i_n],"log",times_log[i_n]))
print(paste("Convergence: score",o1$convergence,"score_2",o1_2$convergence,"log",o2$convergence))
print(paste("Counts: score",o1$counts[[1]],"score_2",o1_2$counts[[1]],"log",o2$counts[[1]]))
print(o1$par)
print(o1_2$par)
print(o2$par)
}

times_score_2/times_log

plot(times_score_2[1:8]/times_log[1:8])
plot(log(times_log[1:8]),log(times_score_2[1:8]))

(log(times_log[8])-log(times_log[5]))/(log(times_score_2[8])-log(times_score_2[5]))
times_score
times_log

o_log_list
o_score_list

##############
tic()
loo_score(m,mu,Qmat)
toc()
tic()
loo_score_2(m,mu,Qmat)
toc()
tic()
loo_score_old(m,mu,Qmat)
toc()
