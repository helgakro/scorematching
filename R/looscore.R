
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


loo_score_sapply <- function(obs, mu, precmat){
  n_obs <- nrow(obs)
  obs <- t(obs)
  n_dim <- nrow(precmat)
  score_vec <- nrow(n_dim)
  Qmu <- precmat%*%mu #matrix multiplication
  Qx <- precmat%*%obs
  return(mean(sapply(c(1:n_dim), function(i)sroot_normal(obs[i,],(Qmu[i]-Qx[i,])/precmat[i,i]+obs[i,],1/sqrt(precmat[i,i])))))
}


loo_score_vectorised <- function(obs, mu, precmat){
  n_obs <- nrow(obs)
  obs <- t(obs)
  n_dim <- nrow(precmat)
  score_vec <- nrow(n_dim)
  return(mean(sroot_normal(c(obs),as.vector(precmat%*%((mu-obs))/diag(precmat))+c(obs),rep(1/sqrt(diag(precmat)),n_obs))))
}

loo_log_score <- function(obs, mu, precmat){
  n_obs <- nrow(obs)
  obs <- t(obs)
  n_dim <- nrow(precmat)
  score_vec <- nrow(n_dim)
  return(-mean(log(dnorm(c(obs),as.vector(precmat%*%((mu-obs))/diag(precmat))+c(obs),rep(1/sqrt(diag(precmat)),n_obs)))))
}


log_dmvn <- function(obs,mu,precmat){
  n_dim <- nrow(precmat)
  n_obs <- nrow(obs)
  if(is.null(n_obs)){
    return(-n_dim/2*log(2*pi)+1/2*log(det(precmat))-1/2*t(obs-mu)%*%precmat%*%(obs-mu))
  }else{
    return(-n_dim/2*log(2*pi)+1/2*log(det(precmat))-1/2*mean(sapply(c(1:n_obs), function(i) as.numeric(t(obs[i,]-mu)%*%precmat%*%(obs[i,]-mu)))))
    #return(-n_dim/2*log(2*pi)+1/2*log(det(precmat))-1/2*mean(diag((obs-mu)%*%precmat%*%t(obs-mu))))
  }
}

loo_score_vectorised_eps <- function(obs, mu, precmat,sigma){
  n_obs <- nrow(obs)
  obs <- t(obs)
  n_dim <- nrow(precmat)
  score_vec <- nrow(n_dim)
  return(mean(sroot_normal(c(obs),as.vector(precmat%*%((mu-obs))/diag(precmat))+c(obs),sqrt(rep(1/diag(precmat),n_obs)+sigma^2))))
}

loo_log_score_eps <- function(obs, mu, precmat,sigma){
  n_obs <- nrow(obs)
  obs <- t(obs)
  n_dim <- nrow(precmat)
  score_vec <- nrow(n_dim)
  return(-mean(log(dnorm(c(obs),as.vector(precmat%*%((mu-obs))/diag(precmat))+c(obs),rep(1/sqrt(diag(precmat)),n_obs)+sigma^2))))
}


log_dmvn_eps <- function(obs,mu,precmat,sigma){
  n_dim <- nrow(precmat)
  n_obs <- nrow(obs)
  if(sigma>0){
    precmat <- precmat + Diagonal(n_dim)/sigma^2
  }
  if(is.null(n_obs)){
    return(-n_dim/2*log(2*pi)+1/2*log(det(precmat))-1/2*t(obs-mu)%*%precmat%*%(obs-mu))
  }else{
    return(-n_dim/2*log(2*pi)+1/2*log(det(precmat))-1/2*mean(sapply(c(1:n_obs), function(i) as.numeric(t(obs[i,]-mu)%*%precmat%*%(obs[i,]-mu)))))
    #return(-n_dim/2*log(2*pi)+1/2*log(det(precmat))-1/2*mean(diag((obs-mu)%*%precmat%*%t(obs-mu))))
  }
}

