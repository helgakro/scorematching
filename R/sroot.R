#' Compute s_root score
#'
#' @param y observations
#' @param params model parameters
#' @param type model type, can be one of "r-Pareto", "norm"
#'
#' @return s_root score
#' @export
#'
#' @examples sroot(c(1,2,3),list("kappa"=1, "tau"=4))
sroot <- function(y,params,type="r-Pareto",a=0,b=1){
  if(type=="r-Pareto"){
    kappa <- params[["kappa"]]
    tau <- params[["tau"]]
    score <- sroot_rpareto(y,kappa,tau)
  }else if(type=="norm"){
    mu<-params[["mu"]]
    sigma<-params[["sigma"]]
    score<-sroot_normal(y,mu,sigma)
  }else if(type=="tnorm"){
    mu<-params[["mu"]]
    sigma<-params[["sigma"]]
    score<-sroot_tnorm(y,mu,sigma,a,b)
  }else if(type=="mctnorm"){
    mu<-params[["mu"]]
    sigma<-params[["sigma"]]
    if(length(y)>1){
      score<-y*0
      for(i in c(1:length(y))){
        score[i]<-sroot_tnorm_mc(y[i],mu=mu,s=sigma,l=a,u=b)
      }
    } else {
      score<-sroot_tnorm_mc(y,mu=mu,s=sigma,l=a,u=b)
    }
  }
  return(score)
}


sroot_rpareto <- function(y,kappa,tau){
  ep <- kappa #todo: formula
  epp <- tau # todo: formula
  score <- ep/sqrt(epp)
  return(score)
}

#sroot for normal distribution
sroot_normal <- function(y,mu,sigma){
  ep <- 2*sigma*dnorm((mu-y)/sigma)+(mu-y)*(2*pnorm((mu-y)/sigma)-1)
  epp <- 2*sigma/sqrt(pi)
  score <- ep/sqrt(epp)
  return(score)
}

#crps for normal distribution
crps_normal <- function(y,mu,sigma){
  ep <- 2*sigma*dnorm((mu-y)/sigma)+(mu-y)*(2*pnorm((mu-y)/sigma)-1)
  epp <- 2*sigma/sqrt(pi)
  score <- ep-epp/2
  return(score)
}


#crps for normal distribution
scrps_normal <- function(y,mu,sigma){
  ep <- 2*sigma*dnorm((mu-y)/sigma)+(mu-y)*(2*pnorm((mu-y)/sigma)-1)
  epp <- 2*sigma/sqrt(pi)
  score <- ep/epp+log(epp)/2
  return(score)
}

#sroot for truncated normal distribution
sroot_tnorm <- function(y,mu,sigma,l,u){
  y<-(y-mu)/sigma
  u<-(u-mu)/sigma
  l<-(l-mu)/sigma
  # ep <- sigma*(2*dnorm(y)+y*(2*pnorm(y)-pnorm(u)-pnorm(l)))
  # epp <- sigma*2*(pnorm(u*sqrt(2))-pnorm(l*sqrt(2)))/sqrt(pi)
  ep <- sigma*(2*dnorm(y)+y*(2*pnorm(y)-pnorm(u)-pnorm(l))-(dnorm(u)+dnorm(l)))
  epp <- sigma*2*((pnorm(u*sqrt(2))-pnorm(l*sqrt(2)))/sqrt(pi)-(pnorm(u)-pnorm(l))*(dnorm(u)+dnorm(l)))
  score <- ep/sqrt(epp)
  return(score)
}


#crps for truncated normal distribution
mycrps_tnorm <- function(y,mu,sigma,l,u){
  y<-(y-mu)/sigma
  u<-(u-mu)/sigma
  l<-(l-mu)/sigma
  score <- sigma/(pnorm(u)-pnorm(l))*(2*dnorm(y)+y*(2*pnorm(y)-pnorm(u)-pnorm(l))) - sigma/(pnorm(u)-pnorm(l))^2*(pnorm(u*sqrt(2))-pnorm(l*sqrt(2)))/sqrt(pi)
  return(score)
}


sroot_norm_mc <- function(y,mu,s){
  n<-100000
  x <- rnorm(n,mu,s)
  x1 <- rnorm(n,mu,s)
  x2 <- rnorm(n,mu,s)
  ep <- mean(abs(x-y))
  epp <- mean(abs(x1-x2))
  score <- ep/sqrt(epp)
  return(score)
}

sroot_tnorm_mc <- function(y,mu,s,l,u){
  n<-1000000
  x <- rnorm(n,mu,s)
  x <- x[x>l & x<u]
  x1 <- rnorm(n,mu,s)
  x1 <- x1[x1>l & x1<u]
  x2 <- rnorm(n,mu,s)
  x2 <- x2[x2>l & x2<u]
  if(length(x1)<length(x2)){
    x2 <- x2[c(1:length(x1))]
  } else {
    x1 <- x1[c(1:length(x2))]
  }
  ep <- mean(abs(x-y))
  epp <- mean(abs(x1-x2))
  score <- ep/sqrt(epp)
  return(score)
}

