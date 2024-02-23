sroot(c(1,2,3),list("kappa"=1, "tau"=4))


sroot(c(1,2,3),list("mu"=1, "sigma"=4),type="norm")

############################
# Normal inference example
############################

n=1000
mutrue=1
sigmatrue=2

y=rnorm(n,mean=mutrue,sd=sigmatrue)
hist(y)

mu0=0
sigma0=1

scorefcn <- function(theta) sum(sroot(y,list("mu"=theta[1], "sigma"=theta[2]),type="norm"))

optim(c(mu0,sigma0),scorefcn)


# Bera saman við likelihood

logliknorm<- function(theta) { sum(-log(dnorm(y,mean=theta[1],sd=theta[2]))) }

optim(c(mu0,sigma0),logliknorm)



############ plot scorefcn

library(plotly)
library(plot3D)
xlist<-seq(0, 5, length.out = 100)
ylist<-seq(1, 6, length.out = 100)
M <- mesh(xlist,ylist)
u<-M$x
v<-M$y
z<-0*u
z2<-0*u
for(i in c(1:length(z))){
  if(i%%10==0){
    print(i)
  }
  z[i]<-scorefcn(c(u[i],v[i]))
  z2[i]<-logliknorm(c(u[i],v[i]))
}


fig <- plot_ly(x=xlist,y=ylist,z = t(z), type = "contour")
fig

fig2 <- plot_ly(x=xlist,y=ylist,z = t(z2), type = "contour")
fig2


#################################
# Repeated experiments
#################################

compare_inference_norm <- function(){
n=1000
mutrue=1
sigmatrue=2

y=rnorm(n,mean=mutrue,sd=sigmatrue)
hist(y)

mu0=0
sigma0=1

scorefcn <- function(theta) sum(sroot(y,list("mu"=theta[1], "sigma"=theta[2]),type="norm"))

o1<-optim(c(mu0,sigma0),scorefcn)


# Bera saman við likelihood

logliknorm<- function(theta) { sum(-log(dnorm(y,mean=theta[1],sd=theta[2]))) }

o2<-optim(c(mu0,sigma0),logliknorm)
return(list("sroot"=o1, "logs"=o2))
}

n_rep <- 100
res_list_sroot <- list()
res_list_logs <- list()
for(i in c(1:n_rep)){
  res_inference_norm <- compare_inference_norm()
  res_list_sroot <- rbind(res_list_sroot,res_inference_norm$sroot$par)
  res_list_logs <- rbind(res_list_logs,res_inference_norm$logs$par)
}

plot(res_list_sroot[,1],res_list_logs[,1],asp=1)
abline(a=0,b=1)
plot(res_list_sroot[,2],res_list_logs[,2],asp=1)
abline(a=0,b=1)






##############################
# Truncated normal example
##############################

library(truncnorm)
n=1000#10000000
mutrue=1
sigmatrue=2
a=0
b=5

y=rtruncnorm(n,mean=mutrue,sd=sigmatrue,a=a,b=b)
#plot(y)
hist(y)

mu0=1
sigma0=1

scorefcn <- function(theta) {
  #sum(scoringRules::crps_tnorm(y,location=theta[1], scale=theta[2], lower=a,upper=b))
  sum(sroot(y,list("mu"=theta[1], "sigma"=theta[2]),type="tnorm", a=a,b=b))
}

optim(c(mu0,sigma0),scorefcn)$par

start_time <- Sys.time()
optim(c(mu0,sigma0),scorefcn)
end_time <- Sys.time()
end_time - start_time


# MC score

# scorefcn_mc <- function(theta) {
#   #sum(scoringRules::crps_tnorm(y,location=theta[1], scale=theta[2], lower=a,upper=b))
#   score <- sum(sroot(y,list("mu"=theta[1], "sigma"=theta[2]),type="mctnorm", a=a,b=b))
#   print(score)
#   return(score)
# }
#
# optim(c(mu0,sigma0),scorefcn_mc)$par
#
# start_time <- Sys.time()
# optim(c(mu0,sigma0),scorefcn_mc)
# end_time <- Sys.time()
# end_time - start_time

# Bera saman við likelihood

logliktrunc <- function(theta) { sum(-log(dtruncnorm(y,a=a,b=b,mean=theta[1],sd=theta[2]))) }

start_time <- Sys.time()
optim(c(mu0,sigma0),logliktrunc)
end_time <- Sys.time()
end_time - start_time


#bera saman við Hyvarinen

library(scoringRules)





############ plot scorefcn

library(plotly)
library(plot3D)
xlist<-seq(-1, 5, length.out = 15)
ylist<-seq(1, 9, length.out = 10)
M <- mesh(xlist,ylist)
u<-M$x
v<-M$y
z<-0*u
z2<-0*u
for(i in c(1:length(z))){
  z[i]<-scorefcn(c(u[i],v[i]))
  z2[i]<-logliktrunc(c(u[i],v[i]))
}


fig <- plot_ly(x=xlist,y=ylist,z = t(z), type = "contour")
fig

fig2 <- plot_ly(x=xlist,y=ylist,z = t(z2), type = "contour")
fig2



#################################
# Repeated experiments
#################################

compare_inference_tnorm <- function(){
  n=1000#10000000
  mutrue=1
  sigmatrue=2
  a=0
  b=5

  y=rtruncnorm(n,mean=mutrue,sd=sigmatrue,a=a,b=b)

  mu0=1
  sigma0=1

  scorefcn <- function(theta) {
    #sum(scoringRules::crps_tnorm(y,location=theta[1], scale=theta[2], lower=a,upper=b))
    sum(sroot(y,list("mu"=theta[1], "sigma"=theta[2]),type="tnorm", a=a,b=b))
  }

  o1 <- optim(c(mu0,sigma0),scorefcn)


  # Bera saman við likelihood

  logliktrunc <- function(theta) { sum(-log(dtruncnorm(y,a=a,b=b,mean=theta[1],sd=theta[2]))) }

  o2<-optim(c(mu0,sigma0),logliktrunc)


  return(list("sroot"=o1, "logs"=o2))
}

n_rep <- 1000
res_tnorm_list_sroot <- list()
res_tnorm_list_logs <- list()
for(i in c(1:n_rep)){
  res_inference_tnorm <- compare_inference_tnorm()
  res_tnorm_list_sroot <- rbind(res_tnorm_list_sroot,res_inference_tnorm$sroot$par)
  res_tnorm_list_logs <- rbind(res_tnorm_list_logs,res_inference_tnorm$logs$par)
}

plot(res_tnorm_list_sroot[,1],res_tnorm_list_logs[,1],asp=1)
abline(a=0,b=1)
plot(res_tnorm_list_sroot[,2],res_tnorm_list_logs[,2],asp=1)
abline(a=0,b=1)



######################################
# Compare exact and mc sroot score
######################################

sroot_tnorm(1,0,1,0,1)
sapply(c(1,2),sroot_tnorm_mc,mu=0,s=1,l=0,u=1)


#######################################
# Multivariate?
#######################################


sum(scoringRules::crps_tnorm(y,location=0, scale=1, lower=a,upper=b))
sum(mycrps_tnorm(y,0,1,a,b))
