
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
  tmp<-Matrix::t(precmat[i,-i]*Matrix::t(x[,-i]-mu[-i]))
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
  obs <- Matrix::t(obs)
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
  obs <- Matrix::t(obs)
  n_dim <- nrow(precmat)
  score_vec <- nrow(n_dim)
  Qmu <- precmat%*%mu #matrix multiplication
  Qx <- precmat%*%obs
  return(mean(sapply(c(1:n_dim), function(i)sroot_normal(obs[i,],(Qmu[i]-Qx[i,])/precmat[i,i]+obs[i,],1/sqrt(precmat[i,i])))))
}


loo_score_vectorised <- function(obs, mu, precmat){
  n_obs <- nrow(obs)
  obs <- Matrix::t(obs)
  n_dim <- nrow(precmat)
  score_vec <- nrow(n_dim)
  return(mean(sroot_normal(as.vector(obs),as.vector(precmat%*%((mu-obs))/Matrix::diag(precmat))+as.vector(obs),rep(1/sqrt(Matrix::diag(precmat)),n_obs))))
}

loo_log_score <- function(obs, mu, precmat){
  n_obs <- nrow(obs)
  obs <- Matrix::t(obs)
  n_dim <- nrow(precmat)
  score_vec <- nrow(n_dim)
  return(mean(log(sqrt(2*pi)*rep(1/sqrt(Matrix::diag(precmat)),n_obs))+0.5*((-as.vector(precmat%*%((mu-obs))/Matrix::diag(precmat)))/rep(1/sqrt(Matrix::diag(precmat)),n_obs))^2))
  #return(-mean(log(stats::dnorm(c(obs),as.vector(precmat%*%((mu-obs))/Matrix::diag(precmat))+c(obs),rep(1/sqrt(Matrix::diag(precmat)),n_obs)))))
}


log_dmvn <- function(obs,mu,precmat){
  n_dim <- nrow(precmat)
  n_obs <- nrow(obs)
  if(is.null(n_obs)){
    return(-n_dim/2*log(2*pi)+1/2*log(Matrix::det(precmat))-1/2*Matrix::t(obs-mu)%*%precmat%*%(obs-mu))
  }else{
    return(-n_dim/2*log(2*pi)+1/2*log(Matrix::det(precmat))-1/2*mean(sapply(c(1:n_obs), function(i) as.numeric(Matrix::t(obs[i,]-mu)%*%precmat%*%(obs[i,]-mu)))))
    #return(-n_dim/2*log(2*pi)+1/2*log(Matrix::det(precmat))-1/2*mean(Matrix::diag((obs-mu)%*%precmat%*%t(obs-mu))))
  }
}

# loo_score_vectorised_eps <- function(obs, mu, precmat,sigma,A){
#   n_obs <- nrow(obs)
#   obs <- t(obs)
#   n_dim <- nrow(precmat)
#   score_vec <- nrow(n_dim)
#   Qeps <- Matrix::Diagonal(n_dim)/sigma^2
#   Qxy <- precmat + t(A)%*%Qeps%*%A
#   invQxy <- solve(Qxy)
#   invQxyi<- inv_sherman_morrison(Q,A[i,],A[i,])/sigma^2
#   return(mean(sroot_normal(c(obs),as.vector(Qeps%*%((mu-obs))/diag(precmat))+c(obs)/sigma^2+(1-1/sigma^2)*c(mu),sqrt(rep(1/diag(precmat),n_obs)+sigma^2))))
# }

loo_score_vectorised_eps_old <- function(obs, mu, Qx,sigma,A){
  n_obs <- nrow(obs) #how many observations of the field
  obs <- Matrix::t(obs) #each column represents an observation of the field
  n_dim <- nrow(Qx) #number of "observed" points in field
  score_vec <- nrow(n_dim)
  Qeps <- Matrix::Diagonal(n_dim)/sigma^2 #noise precision
  Qxy <- Qx + Matrix::t(A)%*%Qeps%*%A #full conditional precision matrix
  invQxy <- solve(Qxy) #inverse of invQxy, to be used by the S-M formula
  Qeps <- Matrix::Diagonal(n_dim-1)/sigma^2 #noise precision for n-1 dimension
  invQe <- solve(Qeps) #inverse noise precision for n-1 dimension
  #invQxyi<- lapply(c(1:n_dim),function(i) inv_sherman_morrison(invQxy,A[i,,drop=FALSE]/sigma,-A[i,,drop=FALSE]/sigma))
  invQxyi<- lapply(c(1:n_dim),function(i) solve(Qx + Matrix::t(A[-i,,drop=FALSE])%*%Qeps%*%A[-i,,drop=FALSE])) #inverse precision matrix, independent on observations
  muxyi <- lapply(c(1:n_dim),function(i) mu + invQxyi[[i]]%*%Matrix::t(A[-i,,drop=FALSE])%*%invQe%*%(obs[-i]-A[-i,,drop=FALSE]%*%mu))
  params<-sapply(c(1:n_dim),function(i) c(as.numeric(A[i,,drop=FALSE]%*%muxyi[[i]]),sqrt(as.numeric(A[i,,drop=FALSE]%*%invQxyi[[i]]%*%Matrix::t(A[i,,drop=FALSE])+sigma^2))))
  return(mean(sroot_normal(c(obs),params[1,],params[2,])))
}

loo_score_vectorised_eps <- function(obs, mu, Qx,sigma,A){
  n_obs <- nrow(obs) #how many observations of the field
  obs <- Matrix::t(obs) #each column represents an observation of the field
  n_dim <- nrow(Qx) #number of "observed" points in field
  score_vec <- nrow(n_dim)
  Qeps <- Matrix::Diagonal(n_dim)/sigma^2 #noise precision
  Qxy <- Qx + Matrix::t(A)%*%Qeps%*%A #full conditional precision matrix
  invQxy <- solve(Qxy) #inverse of invQxy, to be used by the S-M formula
  Qeps <- Matrix::Diagonal(n_dim-1)/sigma^2 #noise precision for n-1 dimension
  invQe <- solve(Qeps) #inverse noise precision for n-1 dimension
  #invQxyi<- lapply(c(1:n_dim),function(i) inv_sherman_morrison(invQxy,A[i,,drop=FALSE]/sigma,-A[i,,drop=FALSE]/sigma))
  invQxyi<- lapply(c(1:n_dim),function(i) solve(Qx + Matrix::t(A[-i,,drop=FALSE])%*%Qeps%*%A[-i,,drop=FALSE])) #inverse precision matrix, independent on observations
  # invQxyi <- array(as.numeric(unlist(invQxyi)),dim=c(n_dim,n_dim,n_dim))
  muxyi <- lapply(c(1:n_dim),function(i) mu + invQxyi[[i]]%*%Matrix::t(A[-i,,drop=FALSE])%*%invQe%*%(obs[-i,]-A[-i,,drop=FALSE]%*%mu))
  params<-do.call(cbind,lapply(c(1:n_dim),function(i) rbind(as.numeric(A[i,,drop=FALSE]%*%muxyi[[i]]),rep(sqrt(as.numeric(A[i,,drop=FALSE]%*%invQxyi[[i]]%*%Matrix::t(A[i,,drop=FALSE])+sigma^2)),n_obs))))
  return(mean(sroot_normal(c(obs),params[1,],params[2,])))
}

loo_score_vectorised_eps <- function(obs, mu, Qx,sigma,A){
  n_obs <- nrow(obs) #how many observations of the field
  obs <- Matrix::t(obs) #each column represents an observation of the field
  n_dim <- nrow(Qx) #number of "observed" points in field
  score_vec <- nrow(n_dim)
  Qeps <- Matrix::Diagonal(n_dim)/sigma^2 #noise precision
  Qxy <- Qx + Matrix::t(A)%*%Qeps%*%A #full conditional precision matrix
  invQxy <- solve(Qxy) #inverse of invQxy, to be used by the S-M formula
  Qeps <- Matrix::Diagonal(n_dim-1)/sigma^2 #noise precision for n-1 dimension
  invQe <- solve(Qeps) #inverse noise precision for n-1 dimension
  #invQxyi<- lapply(c(1:n_dim),function(i) inv_sherman_morrison(invQxy,A[i,,drop=FALSE]/sigma,-A[i,,drop=FALSE]/sigma))
  invQxyi<- lapply(c(1:n_dim),function(i) solve(Qx + Matrix::t(A[-i,,drop=FALSE])%*%Qeps%*%A[-i,,drop=FALSE])) #inverse precision matrix, independent on observations
  # invQxyi <- array(as.numeric(unlist(invQxyi)),dim=c(n_dim,n_dim,n_dim))
  muxyi <- lapply(c(1:n_dim),function(i) mu + invQxyi[[i]]%*%Matrix::t(A[-i,,drop=FALSE])%*%invQe%*%(obs[-i,]-A[-i,,drop=FALSE]%*%mu))
  params<-do.call(cbind,lapply(c(1:n_dim),function(i) rbind(as.numeric(A[i,,drop=FALSE]%*%muxyi[[i]]),rep(sqrt(as.numeric(A[i,,drop=FALSE]%*%invQxyi[[i]]%*%Matrix::t(A[i,,drop=FALSE])+sigma^2)),n_obs))))
  return(mean(sroot_normal(c(obs),params[1,],params[2,])))
}



inv_sherman_morrison<- function(invQ,u,v){
  #return(invQ-(invQ%*%u%*%Matrix::t(v)%*%invQ)/as.numeric(1+Matrix::t(v)%*%invQ%*%u))
  return(invQ-(invQ%*%(Matrix::t(u)%*%v)%*%invQ)/as.numeric(1+v%*%invQ%*%Matrix::t(u)))
}

inv_sherman_morrison_I<- function(invQ,i,sigma){
  return(invQ+(invQ[,i]%*%Matrix::t(invQ[i,])/sigma^2)/(1+invQ[i,i]))
}

loo_log_score_eps <- function(obs, mu, Qx,sigma,A){
  n_obs <- nrow(obs) #how many observations of the field
  obs <- Matrix::t(obs) #each column represents an observation of the field
  n_dim <- nrow(Qx) #number of "observed" points in field
  score_vec <- nrow(n_dim)
  Qeps <- Matrix::Diagonal(n_dim)/sigma^2 #noise precision
  Qxy <- Qx + Matrix::t(A)%*%Qeps%*%A #full conditional precision matrix
  invQxy <- solve(Qxy) #inverse of invQxy, to be used by the S-M formula
  Qeps <- Matrix::Diagonal(n_dim-1)/sigma^2 #noise precision for n-1 dimension
  invQe <- solve(Qeps) #inverse noise precision for n-1 dimension
  #invQxyi<- lapply(c(1:n_dim),function(i) inv_sherman_morrison(invQxy,A[i,,drop=FALSE]/sigma,-A[i,,drop=FALSE]/sigma))
  invQxyi<- lapply(c(1:n_dim),function(i) solve(Qx + Matrix::t(A[-i,,drop=FALSE])%*%Qeps%*%A[-i,,drop=FALSE])) #inverse precision matrix, independent on observations
  # invQxyi <- array(as.numeric(unlist(invQxyi)),dim=c(n_dim,n_dim,n_dim))
  muxyi <- lapply(c(1:n_dim),function(i) mu + invQxyi[[i]]%*%Matrix::t(A[-i,,drop=FALSE])%*%invQe%*%(obs[-i,]-A[-i,,drop=FALSE]%*%mu))
  params<-do.call(cbind,lapply(c(1:n_dim),function(i) rbind(as.numeric(A[i,,drop=FALSE]%*%muxyi[[i]]),rep(sqrt(as.numeric(A[i,,drop=FALSE]%*%invQxyi[[i]]%*%Matrix::t(A[i,,drop=FALSE])+sigma^2)),n_obs))))
  #return(-mean(log(stats::dnorm(c(obs),params[1,],params[2,]))))
  return(mean(log(sqrt(2*pi)*params[2,])+0.5*((c(obs)-params[1,])/params[2,])^2))
}


log_dmvn_eps <- function(obs,mu,precmat,sigma){
  n_dim <- nrow(precmat)
  n_obs <- nrow(obs)
  if(sigma>0){
    precmat <- precmat + Matrix::Diagonal(n_dim)/sigma^2
  }
  if(is.null(n_obs)){
    return(-n_dim/2*log(2*pi)+1/2*log(Matrix::det(precmat))-1/2*Matrix::t(obs-mu)%*%precmat%*%(obs-mu))
  }else{
    return(-n_dim/2*log(2*pi)+1/2*log(Matrix::det(precmat))-1/2*mean(sapply(c(1:n_obs), function(i) as.numeric(Matrix::t(obs[i,]-mu)%*%precmat%*%(obs[i,]-mu)))))
    #return(-n_dim/2*log(2*pi)+1/2*log(Matrix::det(precmat))-1/2*mean(Matrix::diag((obs-mu)%*%precmat%*%t(obs-mu))))
  }
}

