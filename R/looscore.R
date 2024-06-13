
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



loo_score_vectorised <- function(obs, mu, precmat,score="sroot"){
  n_obs <- nrow(obs)
  obs <- Matrix::t(obs)
  n_dim <- nrow(precmat)
  score_vec <- nrow(n_dim)
  #return(mean(sroot_normal(as.vector(obs),as.vector(precmat%*%((mu-obs))/Matrix::diag(precmat))+as.vector(obs),rep(1/sqrt(Matrix::diag(precmat)),n_obs))))
  if(score=="sroot"){
    return(mean(sroot_normal(as.vector(obs),as.vector(precmat%*%((mu-obs))/Matrix::diag(precmat))+as.vector(obs),rep(1/sqrt(Matrix::diag(precmat)),n_obs))))
  }else if(score=="crps"){
    return(mean(crps_normal(as.vector(obs),as.vector(precmat%*%((mu-obs))/Matrix::diag(precmat))+as.vector(obs),rep(1/sqrt(Matrix::diag(precmat)),n_obs))))
  }else if(score=="scrps"){
    return(mean(scrps_normal(as.vector(obs),as.vector(precmat%*%((mu-obs))/Matrix::diag(precmat))+as.vector(obs),rep(1/sqrt(Matrix::diag(precmat)),n_obs))))
  }else if(score=="rcrps"){
    return(mean(rcrps_normal(as.vector(obs),as.vector(precmat%*%((mu-obs))/Matrix::diag(precmat))+as.vector(obs),rep(1/sqrt(Matrix::diag(precmat)),n_obs),rep(2,n_obs)))) #todo: insert c in a better way
  }else{
    print(paste(score," is not a supported score type"))
    return(NULL)
  }
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
    #return(-n_dim/2*log(2*pi)+1/2*log(Matrix::det(precmat))-1/2*Matrix::t(obs-mu)%*%precmat%*%(obs-mu))
    return(-n_dim/2*log(2*pi)+sum(log(Matrix::diag(Matrix::chol(precmat))))-1/2*Matrix::t(obs-mu)%*%precmat%*%(obs-mu))
  }else{
    return(-n_dim/2*log(2*pi)+sum(log(Matrix::diag(Matrix::chol(precmat))))-1/2*mean(sapply(c(1:n_obs), function(i) as.numeric(Matrix::t(obs[i,]-mu)%*%precmat%*%(obs[i,]-mu)))))
    #return(-n_dim/2*log(2*pi)+1/2*log(Matrix::det(precmat))-1/2*mean(sapply(c(1:n_obs), function(i) as.numeric(Matrix::t(obs[i,]-mu)%*%precmat%*%(obs[i,]-mu)))))
    #return(-n_dim/2*log(2*pi)+1/2*log(Matrix::det(precmat))-1/2*mean(Matrix::diag((obs-mu)%*%precmat%*%t(obs-mu))))
  }
}



loo_score_vectorised_eps <- function(obs, mu, Qx,sigma,A,score="sroot"){
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
  if(score=="sroot"){
    return(mean(sroot_normal(c(obs),params[1,],params[2,])))
  }else if(score=="crps"){
    return(mean(crps_normal(c(obs),params[1,],params[2,])))
  }else if(score=="scrps"){
    return(mean(crps_normal(c(obs),params[1,],params[2,])))
  }else{
    print(paste(score," is not a supported score type"))
    return(NULL)
  }

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



get_rmse <- function(obs, mu, precmat){
  n_obs <- nrow(obs)
  obs <- Matrix::t(obs)
  n_dim <- nrow(precmat)
  score_vec <- nrow(n_dim)
  return(sqrt(mean((as.vector(obs)-(as.vector(precmat%*%((mu-obs))/Matrix::diag(precmat))+as.vector(obs)))^2)))
}
