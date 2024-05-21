inference_fix_resp <- function(m,spde,n_mesh,Q,n_outlier=0,outlier_val=NULL,A=NULL){

  n_x <-  n_mesh #narr_rep[i_n]
  if(is.null(A)){
    A<-Diagonal(n_x)
  }
  n_y <- nrow(A)
  I<-Diagonal(n_y)
  mu<- rep(0,n_x)


  my_obj_func_3 <- function(par){
    print(par)

    theta <- par
    mu <- rep(0,n_x)
    #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
    Qx <- inla.spde.precision(spde, theta=theta)
    Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))

    Qtheta <- Qx
    muy <- A%*%mu
    #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
    score <- loo_score_vectorised(m,muy,Qtheta)
    print(score)
    return(score)
  }

  # my_log_score_obj_func <- function(par){
  #   theta <- par
  #   mu <- rep(0,n)
  #   Qx <- inla.spde.precision(spde, theta=theta)
  #   mux <- mu
  #   score <- loo_log_score_eps(m,mux,Qx,sigma_val,A)
  #   print(score)
  #   return(score)
  # }

  my_log_score_obj_func <- function(par){
    print(par)

    theta <- par
    mu <- rep(0,n_x)
    #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
    Qx <- inla.spde.precision(spde, theta=theta)
    Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))

    Qtheta <- Qx
    muy <- A%*%mu
    #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
    score <- loo_log_score(m,muy,Qtheta)
    print(score)
    return(score)
  }

  #
  # my_log_obj_func <- function(par){
  #   theta <- par
  #   mu <- rep(0,n)
  #   Qxy <- inla.spde.precision(spde, theta=theta)+I*sigma_val^2
  #   muxy <- mu
  #   if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
  #   return(-log_dmvn_eps(m,muxy,Qxy,sigma_val))
  # }


  my_log_obj_func <- function(par){
    print(par)

    theta <- par
    mu <- rep(0,n_x)
    #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
    Qx <- inla.spde.precision(spde, theta=theta)
    Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))
    Qtheta <- Qx
    muy <- A%*%mu
    #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
    score <- -log_dmvn(m,muy,Qtheta)
    print(score)
    return(score)
  }


  kappa0<--2
  tau0<-0.1


  starttime <- Sys.time()
  o1<-optim(par=c(kappa0,tau0),my_obj_func_3,control=list(maxit=50000))

  endtime <- Sys.time()
  times_score_rep<-difftime(endtime,starttime, units="secs")
  o_score_list_rep<-o1

  o2<-NULL
  times_log_rep<-NULL
  # starttime <- Sys.time()
  # o2<-optim(par=c(kappa0,tau0),my_log_obj_func,control=list(maxit=50000))
  # endtime <- Sys.time()
  # times_log_rep<-difftime(endtime,starttime, units="secs")
  # o_log_list_rep<-o2
  #
  # score_ll_est <- my_obj_func_3(o2$par)

  o3<-NULL
  times_log_score_rep<-NULL
  score_ll_est<-NULL
  # starttime <- Sys.time()
  # o3<-optim(par=c(kappa0,tau0),my_log_score_obj_func,control=list(maxit=50000))
  # endtime <- Sys.time()
  # times_log_score_rep<-difftime(endtime,starttime, units="secs")
  # o_log_score_list_rep<-o3


  return(list(o1=o1,o2=o2,o3=o3,t1=times_score_rep,t2=times_log_rep,t3=times_log_score_rep,score_ll=score_ll_est))
}


inference_norm_resp <- function(m,spde,n_mesh,Q,n_outlier=0,outlier_val=NULL,sigma_val=1,A=NULL,sroot=TRUE,ll=TRUE,slog=TRUE,crps=FALSE,scrps=FALSE){
  n_x <-  n_mesh #narr_rep[i_n]
  if(is.null(A)){
    A<-Diagonal(n_x)
  }
  n_y <- nrow(A)
  I<-Diagonal(n_y)
  mu<- rep(0,n_x)


  my_obj_func_3 <- function(par,scoretype="sroot"){
    print(par)

    sig <- exp(par[1])
    theta <- par[-1]
    mu <- rep(0,n_x)
    #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
    Qx <- inla.spde.precision(spde, theta=theta)
    Qeps <- I/sig^2

    Qtheta <- Qeps-Qeps%*%A%*%solve(Qx+Matrix::t(A)%*%Qeps%*%A)%*%Matrix::t(A)%*%Qeps
    muy <- A%*%mu
    #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
    score <- loo_score_vectorised(m,muy,Qtheta,scoretype)
    print(score)
    return(score)
  }

  # my_log_score_obj_func <- function(par){
  #   theta <- par
  #   mu <- rep(0,n)
  #   Qx <- inla.spde.precision(spde, theta=theta)
  #   mux <- mu
  #   score <- loo_log_score_eps(m,mux,Qx,sigma_val,A)
  #   print(score)
  #   return(score)
  # }

  my_log_score_obj_func <- function(par){
    print(par)

    sig <- exp(par[1])
    theta <- par[-1]
    mu <- rep(0,n_x)
    #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
    Qx <- inla.spde.precision(spde, theta=theta)
    Qeps <- I/sig^2

    Qtheta <- Qeps-Qeps%*%A%*%solve(Qx+Matrix::t(A)%*%Qeps%*%A)%*%Matrix::t(A)%*%Qeps
    muy <- A%*%mu
    #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
    score <- loo_log_score(m,muy,Qtheta)
    print(score)
    return(score)
  }

  #
  # my_log_obj_func <- function(par){
  #   theta <- par
  #   mu <- rep(0,n)
  #   Qxy <- inla.spde.precision(spde, theta=theta)+I*sigma_val^2
  #   muxy <- mu
  #   if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
  #   return(-log_dmvn_eps(m,muxy,Qxy,sigma_val))
  # }


  my_log_obj_func <- function(par){
    print(par)
    sig <- exp(par[1])
    theta <- par[-1]
    mu <- rep(0,n_x)
    #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
    Qx <- inla.spde.precision(spde, theta=theta)
    Qeps <- I/sig^2
    Qtheta <- Qeps-Qeps%*%A%*%solve(Qx+Matrix::t(A)%*%Qeps%*%A)%*%Matrix::t(A)%*%Qeps
    muy <- A%*%mu
    #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
    score <- -log_dmvn(m,muy,Qtheta)
    print(score)
    return(score)
  }

  sig0 <- log(0.1)
  kappa0<--2
  tau0<-0.1

  if(sroot){
  starttime <- Sys.time()
  o1<-optim(par=c(sig0,kappa0,tau0),my_obj_func_3,control=list(maxit=50000))

  endtime <- Sys.time()
  times_score_rep<-difftime(endtime,starttime, units="secs")
  o_score_list_rep<-o1
  }else{
    o1<-NULL
    times_score_rep<-NULL
  }
  if(ll){
  starttime <- Sys.time()
  o2<-optim(par=c(sig0,kappa0,tau0),my_log_obj_func,control=list(maxit=50000))
  endtime <- Sys.time()
  times_log_rep<-difftime(endtime,starttime, units="secs")
  o_log_list_rep<-o2

  score_ll_est <- my_obj_func_3(o2$par)
  }else{
  o2<-NULL
  times_log_rep<-NULL
  score_ll_est<-NULL
  }
  if(slog){
  starttime <- Sys.time()
  o3<-optim(par=c(sig0,kappa0,tau0),my_log_score_obj_func,control=list(maxit=50000))
  endtime <- Sys.time()
  times_log_score_rep<-difftime(endtime,starttime, units="secs")
  o_log_score_list_rep<-o3
  }else{
  o3<-NULL
  times_log_score_rep<-NULL
  }

  if(crps){
    starttime <- Sys.time()
    o4<-optim(par=c(sig0,kappa0,tau0),my_obj_func_3,scoretype="crps",control=list(maxit=50000))

    endtime <- Sys.time()
    times_score_rep_crps<-difftime(endtime,starttime, units="secs")
    o_score_rep_crps<-o4
  }else{
    o4<-NULL
    times_score_rep_crps<-NULL
  }

  if(scrps){
    starttime <- Sys.time()
    o5<-optim(par=c(sig0,kappa0,tau0),my_obj_func_3,scoretype="scrps",control=list(maxit=50000))

    endtime <- Sys.time()
    times_score_rep_scrps<-difftime(endtime,starttime, units="secs")
    o_score_rep_scrps<-o5
  }else{
    o5<-NULL
    times_score_rep_scrps<-NULL
  }

  return(list(o1=o1,o2=o2,o3=o3,o4=o4,o5=o5,t1=times_score_rep,t2=times_log_rep,t3=times_log_score_rep,score_ll=score_ll_est,t4=times_score_rep_crps,t5=times_score_rep_scrps))
}

inference_norm_resp_old <- function(m,spde,n_mesh,Q,n_outlier=0,outlier_val=NULL,sigma_val=1,A=NULL){

    n_x <-  n_mesh #narr_rep[i_n]
    if(is.null(A)){
      A<-Diagonal(n_x)
    }
    n_y <- nrow(A)
    I<-Diagonal(n_y)
    mu<- rep(0,n_x)


    my_obj_func_3 <- function(par){
      print(par)

      sig <- exp(par[1])
      theta <- par[-1]
      mu <- rep(0,n_x)
      #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
      Qx <- inla.spde.precision(spde, theta=theta)
      Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))
      Qeps <- I/sig^2

      Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
      muy <- A%*%mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- loo_score_vectorised(m,muy,Qtheta)
      print(score)
      return(score)
    }

    # my_log_score_obj_func <- function(par){
    #   theta <- par
    #   mu <- rep(0,n)
    #   Qx <- inla.spde.precision(spde, theta=theta)
    #   mux <- mu
    #   score <- loo_log_score_eps(m,mux,Qx,sigma_val,A)
    #   print(score)
    #   return(score)
    # }

    my_log_score_obj_func <- function(par){
      print(par)

      sig <- exp(par[1])
      theta <- par[-1]
      mu <- rep(0,n_x)
      #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
      Qx <- inla.spde.precision(spde, theta=theta)
      Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))
      Qeps <- I/sig^2

      Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
      muy <- A%*%mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- loo_log_score(m,muy,Qtheta)
      print(score)
      return(score)
    }

    #
    # my_log_obj_func <- function(par){
    #   theta <- par
    #   mu <- rep(0,n)
    #   Qxy <- inla.spde.precision(spde, theta=theta)+I*sigma_val^2
    #   muxy <- mu
    #   if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
    #   return(-log_dmvn_eps(m,muxy,Qxy,sigma_val))
    # }


    my_log_obj_func <- function(par){
      print(par)
      sig <- exp(par[1])
      theta <- par[-1]
      mu <- rep(0,n_x)
      #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
      Qx <- inla.spde.precision(spde, theta=theta)
      Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))
      Qeps <- I/sig^2
      Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
      muy <- A%*%mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- -log_dmvn(m,muy,Qtheta)
      print(score)
      return(score)
    }

    sig0 <- log(1)
    kappa0<--2
    tau0<-0.1


    starttime <- Sys.time()
    o1<-optim(par=c(sig0,kappa0,tau0),my_obj_func_3,control=list(maxit=50000))

    endtime <- Sys.time()
    times_score_rep<-difftime(endtime,starttime, units="secs")
    o_score_list_rep<-o1

    starttime <- Sys.time()
    o2<-optim(par=c(sig0,kappa0,tau0),my_log_obj_func,control=list(maxit=50000))
    endtime <- Sys.time()
    times_log_rep<-difftime(endtime,starttime, units="secs")
    o_log_list_rep<-o2

    score_ll_est <- my_obj_func_3(o2$par)

    starttime <- Sys.time()
    o3<-optim(par=c(sig0,kappa0,tau0),my_log_score_obj_func,control=list(maxit=50000))
    endtime <- Sys.time()
    times_log_score_rep<-difftime(endtime,starttime, units="secs")
    o_log_score_list_rep<-o3




  return(list(o1=o1,o2=o2,o3=o3))
}


inference_norm_resp_old2 <- function(m,spde,n_mesh,Q,n_outlier=0,outlier_val=NULL,sigma_val=1,A=NULL){

    n_x <-  n_mesh #narr_rep[i_n]
    if(is.null(A)){
      A<-Diagonal(n_x)
    }
    n_y <- nrow(A)
    I<-Diagonal(n_y)
    mu<- rep(0,n_x)


    my_obj_func_3 <- function(par){
      print(par)

      sig <- exp(par[1])
      theta <- par[-1]
      mu <- rep(0,n_x)
      #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
      Qx <- inla.spde.precision(spde, theta=theta)
      Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))
      Qeps <- I/sig^2

      Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
      muy <- A%*%mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- loo_score_vectorised(m,muy,Qtheta)
      print(score)
      return(score)
    }

    # my_log_score_obj_func <- function(par){
    #   theta <- par
    #   mu <- rep(0,n)
    #   Qx <- inla.spde.precision(spde, theta=theta)
    #   mux <- mu
    #   score <- loo_log_score_eps(m,mux,Qx,sigma_val,A)
    #   print(score)
    #   return(score)
    # }

    my_log_score_obj_func <- function(par){
      print(par)

      sig <- exp(par[1])
      theta <- par[-1]
      mu <- rep(0,n_x)
      #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
      Qx <- inla.spde.precision(spde, theta=theta)
      Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))
      Qeps <- I/sig^2

      Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
      muy <- A%*%mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- loo_log_score(m,muy,Qtheta)
      print(score)
      return(score)
    }

    #
    # my_log_obj_func <- function(par){
    #   theta <- par
    #   mu <- rep(0,n)
    #   Qxy <- inla.spde.precision(spde, theta=theta)+I*sigma_val^2
    #   muxy <- mu
    #   if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
    #   return(-log_dmvn_eps(m,muxy,Qxy,sigma_val))
    # }


    my_log_obj_func <- function(par){
      print(par)
      sig <- exp(par[1])
      theta <- par[-1]
      mu <- rep(0,n_x)
      #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
      Qx <- inla.spde.precision(spde, theta=theta)
      Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))
      Qeps <- I/sig^2
      Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
      muy <- A%*%mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- -log_dmvn(m,muy,Qtheta)
      print(score)
      return(score)
    }

    sig0 <- log(1)
    kappa0<--2
    tau0<-0.1


    starttime <- Sys.time()
    o1<-optim(par=c(sig0,kappa0,tau0),my_obj_func_3,control=list(maxit=50000))

    endtime <- Sys.time()
    times_score_rep<-difftime(endtime,starttime, units="secs")
    o_score_list_rep<-o1

    starttime <- Sys.time()
    o2<-optim(par=c(sig0,kappa0,tau0),my_log_obj_func,control=list(maxit=50000))
    endtime <- Sys.time()
    times_log_rep<-difftime(endtime,starttime, units="secs")
    o_log_list_rep<-o2

    score_ll_est <- my_obj_func_3(o2$par)

    starttime <- Sys.time()
    o3<-optim(par=c(sig0,kappa0,tau0),my_log_score_obj_func,control=list(maxit=50000))
    endtime <- Sys.time()
    times_log_score_rep<-difftime(endtime,starttime, units="secs")
    o_log_score_list_rep<-o3




  return(list(o1=o1,o2=o2,o3=o3))
}

repeated_score_norm_resp <- function(res,spde,n_mesh,n_rep,Q,n_outlier=0,outlier_val=NULL,sigma_val=1,A=NULL){
  narr_rep <- rep(n_mesh,n_rep)
  times_score_rep <- rep(0,n_rep)
  times_log_rep <- rep(0,n_rep)
  score_ll_est <- rep(0,n_rep)
  times_log_score_rep <- rep(0,n_rep)
  o_score_list_rep <- vector("list", n_rep)
  o_log_list_rep <- vector("list", n_rep)
  o_log_score_list_rep <- vector("list", n_rep)
  times_score_rep_2 <- rep(0,n_rep)
  o_score_list_rep_2 <- vector("list", n_rep)


  for(i_n in c(1:n_rep)){

    print(paste("Iteration",i_n))

    n_x <-  n_mesh #narr_rep[i_n]
    if(is.null(A)){
      A<-Diagonal(n_x)
    }
    n_y <- nrow(A)
    I<-Diagonal(n_y)
    mu<- rep(0,n_x)
    # I<-Diagonal(n)
    # A<-Diagonal(n)
    # m<-t(inla.qsample(n=1, Q = Q, mu=mu)) #observations of latent field
    # m <- m+rnorm(n=n,mean = 0,sd = sigma_val) #added noise
    n_sample <- 10
    m<-inla.qsample(n=n_sample, Q = Q, mu=mu) #observations of latent field
    m <- A%*%m+rnorm(n=n_y*n_sample,mean = 0,sd = sigma_val) #added noise
    m <- Matrix::t(m)
    #m[,1]<-4
    if(n_outlier>0){ #add outliers
      m[ matrix(c(1:n_outlier,sample(1:n,n_outlier)),ncol=2)]<-outlier_val #set random observation at each sample to 4.
    }


    my_obj_func_3_x <- function(par){
      print(par)
      theta <- par
      mu <- rep(0,n)
      #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
      Qx <- inla.spde.precision(spde, theta=theta)
      mux <- mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- loo_score_vectorised_eps(m,mux,Qx,sigma_val,A)
      print(score)
      return(score)
    }

    my_obj_func_3 <- function(par){
      print(par)

      sig <- exp(par[1])
      theta <- par[-1]
      mu <- rep(0,n_x)
      #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
      Qx <- inla.spde.precision(spde, theta=theta)
      Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))
      Qeps <- Qeps <- I/sig^2

      Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
      muy <- A%*%mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- loo_score_vectorised(m,muy,Qtheta)
      print(score)
      return(score)
    }

    # my_log_score_obj_func <- function(par){
    #   theta <- par
    #   mu <- rep(0,n)
    #   Qx <- inla.spde.precision(spde, theta=theta)
    #   mux <- mu
    #   score <- loo_log_score_eps(m,mux,Qx,sigma_val,A)
    #   print(score)
    #   return(score)
    # }

    my_log_score_obj_func <- function(par){
      print(par)

      sig <- exp(par[1])
      theta <- par[-1]
      mu <- rep(0,n_x)
      #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
      Qx <- inla.spde.precision(spde, theta=theta)
      Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))
      Qeps <- I/sig^2

      Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
      muy <- A%*%mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- loo_log_score(m,muy,Qtheta)
      print(score)
      return(score)
    }

    #
    # my_log_obj_func <- function(par){
    #   theta <- par
    #   mu <- rep(0,n)
    #   Qxy <- inla.spde.precision(spde, theta=theta)+I*sigma_val^2
    #   muxy <- mu
    #   if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
    #   return(-log_dmvn_eps(m,muxy,Qxy,sigma_val))
    # }


    my_log_obj_func <- function(par){
      print(par)
      sig <- exp(par[1])
      theta <- par[-1]
      mu <- rep(0,n_x)
      #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
      Qx <- inla.spde.precision(spde, theta=theta)
      Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))
      Qeps <- I/sig^2
      Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
      muy <- A%*%mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- -log_dmvn(m,muy,Qtheta)
      print(score)
      return(score)
    }



    o_score_list_rep[[i_n]]<-c(my_obj_func_3(res$o_sroot[[i_n]]$par),my_log_obj_func(res$o_sroot[[i_n]]$par),my_log_score_obj_func(res$o_sroot[[i_n]]$par))


    o_log_list_rep[[i_n]]<-c(my_obj_func_3(res$o_ll[[i_n]]$par),my_log_obj_func(res$o_ll[[i_n]]$par),my_log_score_obj_func(res$o_ll[[i_n]]$par))


    o_log_score_list_rep[[i_n]]<-c(my_obj_func_3(res$o_slog[[i_n]]$par),my_log_obj_func(res$o_slog[[i_n]]$par),my_log_score_obj_func(res$o_slog[[i_n]]$par))


  }

  return(list(o_sroot =o_score_list_rep ,
              o_ll =o_log_list_rep ,
              o_slog=o_log_score_list_rep))
}


repeated_inference_norm_resp <- function(spde,n_mesh,n_rep,Q,n_outlier=0,outlier_val=NULL,sigma_val=1,A=NULL,Atest=NULL,tnu=NULL,scoretypes=c("sroot","ll","slog")){
  narr_rep <- rep(n_mesh,n_rep)
  times_score_rep <- rep(0,n_rep)
  times_log_rep <- rep(0,n_rep)
  score_ll_est <- rep(0,n_rep)
  times_score_rep_crps <- rep(0,n_rep)
  times_score_rep_scrps <- rep(0,n_rep)
  times_log_score_rep <- rep(0,n_rep)
  o_score_list_rep <- vector("list", n_rep)
  o_log_list_rep <- vector("list", n_rep)
  o_log_score_list_rep <- vector("list", n_rep)
  o_score_list_rep_crps <- vector("list", n_rep)
  o_score_list_rep_scrps <- vector("list", n_rep)
  times_score_rep_2 <- rep(0,n_rep)
  o_score_list_rep_2 <- vector("list", n_rep)

  pred_score_o_sroot <- rep(0,n_rep)
  pred_score_o_ll <- rep(0,n_rep)

  for(i_n in c(1:n_rep)){

    print(paste("Iteration",i_n))

    n_x <-  n_mesh #narr_rep[i_n]
    if(is.null(A)){
      A<-Diagonal(n_x)
    }
    n_y <- nrow(A)
    I<-Diagonal(n_y)
    mu<- rep(0,n_x)
    # I<-Diagonal(n)
    # A<-Diagonal(n)
    # m<-t(inla.qsample(n=1, Q = Q, mu=mu)) #observations of latent field
    # m <- m+rnorm(n=n,mean = 0,sd = sigma_val) #added noise
    n_sample <- 10
    mfield<-inla.qsample(n=n_sample, Q = Q, mu=mu) #observations of latent field
    if(is.null(tnu)){
      m <- A%*%mfield+rnorm(n=n_y*n_sample,mean = 0,sd = sigma_val) #added noise
    }else if(tnu>0){
      m <- A%*%mfield+sigma_val*rt(n=n_y*n_sample,df=tnu)
    }else{
      m <- A%*%mfield
    }

    m <- Matrix::t(m)


    #m[,1]<-4
    if(n_outlier>0){ #add outliers
      m[ matrix(c(1:n_outlier,sample(1:n_y,n_outlier)),ncol=2)]<-abs(m[ matrix(c(1:n_outlier,sample(1:n_y,n_outlier)),ncol=2)])+outlier_val #set random observation at each sample to the absolute observed value plus 4.
      #m[ matrix(c(1:n_outlier,sample(1:n_y,n_outlier)),ncol=2)]<-outlier_val #set random observation at each sample to 4.
    }


    my_obj_func_3_x <- function(par){
      print(par)
      theta <- par
      mu <- rep(0,n)
      #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
      Qx <- inla.spde.precision(spde, theta=theta)
      mux <- mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- loo_score_vectorised_eps(m,mux,Qx,sigma_val,A)
      print(score)
      return(score)
    }

    my_obj_func_3 <- function(par,A,m,scoretype="sroot"){
      print(par)

      sig <- exp(par[1])
      theta <- par[-1]
      mu <- rep(0,n_x)
      #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
      Qx <- inla.spde.precision(spde, theta=theta)
      Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))
      Qeps <- Qeps <- I/sig^2

      Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
      muy <- A%*%mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- loo_score_vectorised(m,muy,Qtheta,scoretype)
      print(score)
      return(score)
    }

    # my_log_score_obj_func <- function(par){
    #   theta <- par
    #   mu <- rep(0,n)
    #   Qx <- inla.spde.precision(spde, theta=theta)
    #   mux <- mu
    #   score <- loo_log_score_eps(m,mux,Qx,sigma_val,A)
    #   print(score)
    #   return(score)
    # }

    my_log_score_obj_func <- function(par,A,m){
      print(par)

      sig <- exp(par[1])
      theta <- par[-1]
      mu <- rep(0,n_x)
      #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
      Qx <- inla.spde.precision(spde, theta=theta)
      Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))
      Qeps <- I/sig^2

      Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
      muy <- A%*%mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- loo_log_score(m,muy,Qtheta)
      print(score)
      return(score)
    }

    #
    # my_log_obj_func <- function(par){
    #   theta <- par
    #   mu <- rep(0,n)
    #   Qxy <- inla.spde.precision(spde, theta=theta)+I*sigma_val^2
    #   muxy <- mu
    #   if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
    #   return(-log_dmvn_eps(m,muxy,Qxy,sigma_val))
    # }


    my_log_obj_func <- function(par,A,m){
      print(par)
      sig <- exp(par[1])
      theta <- par[-1]
      mu <- rep(0,n_x)
      #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
      Qx <- inla.spde.precision(spde, theta=theta)
      Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))
      Qeps <- I/sig^2
      Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
      muy <- A%*%mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- -log_dmvn(m,muy,Qtheta)
      print(score)
      return(score)
    }

    sig0 <- log(1)
    kappa0<--2
    tau0<-0.1

    if("sroot"%in%scoretypes){
    starttime <- Sys.time()
    o1<-optim(par=c(sig0,kappa0,tau0),my_obj_func_3,A=A,m=m,scoretype="sroot",control=list(maxit=50000))

    endtime <- Sys.time()
    times_score_rep[i_n]<-difftime(endtime,starttime, units="secs")
    o_score_list_rep[[i_n]]<-o1
    }else{
      times_score_rep[i_n]<-NULL
      o_score_list_rep[[i_n]]<-NULL
    }

    if("ll"%in%scoretypes){
    starttime <- Sys.time()
    o2<-optim(par=c(sig0,kappa0,tau0),my_log_obj_func,A=A,m=m,control=list(maxit=50000))
    endtime <- Sys.time()
    times_log_rep[i_n]<-difftime(endtime,starttime, units="secs")
    o_log_list_rep[[i_n]]<-o2

    score_ll_est[[i_n]] <- my_obj_func_3(o2$par,A=A,m=m)
    }else{
      o_log_list_rep[[i_n]]<-NULL

      score_ll_est[[i_n]] <- NULL
    }

    if("slog"%in%scoretypes){
      starttime <- Sys.time()
    o3<-optim(par=c(sig0,kappa0,tau0),my_log_score_obj_func,A=A,m=m,control=list(maxit=50000))
    endtime <- Sys.time()
    times_log_score_rep[i_n]<-difftime(endtime,starttime, units="secs")
    o_log_score_list_rep[[i_n]]<-o3
    }else{
      times_log_score_rep[i_n]<- -1
      o_log_score_list_rep[[i_n]]<- -1
    }


    if(!is.null(Atest)){
      n_test <- nrow(Atest)
      mtest <- Atest%*%mfield+rnorm(n=n_test*n_sample,mean = 0,sd = sigma_val) #added noise
      mtest <- Matrix::t(mtest)

      pred_score_o_sroot[[i_n]] <- my_obj_func_3(o1$par,A=Atest,m=mtest)
      pred_score_o_ll[[i_n]] <- my_obj_func_3(o2$par,A=Atest,m=mtest)
    }

    if("crps"%in%scoretypes){
      starttime <- Sys.time()
      o4<-optim(par=c(sig0,kappa0,tau0),my_obj_func_3,A=A,m=m,scoretype="crps",control=list(maxit=50000))

      endtime <- Sys.time()
      times_score_rep_crps[i_n]<-difftime(endtime,starttime, units="secs")
      o_score_list_rep_crps[[i_n]]<-o4
    }else{
      times_score_rep_crps[i_n]<-NULL
      o_score_list_rep_crps[[i_n]]<-NULL
    }

    if("scrps"%in%scoretypes){
      starttime <- Sys.time()
      o5<-optim(par=c(sig0,kappa0,tau0),my_obj_func_3,A=A,m=m,scoretype="scrps",control=list(maxit=50000))

      endtime <- Sys.time()
      times_score_rep_scrps[i_n]<-difftime(endtime,starttime, units="secs")
      o_score_list_rep_scrps[[i_n]]<-o5
    }else{
      times_score_rep_scrps[i_n]<-NULL
      o_score_list_rep_scrps[[i_n]]<-NULL
    }

  }

  return(list(times_sroot=times_score_rep,
              times_ll =times_log_rep ,
              times_slog=times_log_score_rep,
              times_crps=times_score_rep_crps,
              times_scrps=times_score_rep_scrps,
              o_sroot =o_score_list_rep ,
              o_ll =o_log_list_rep ,
              o_slog=o_log_score_list_rep,
              o_crps =o_score_list_rep_crps ,
              o_scrps =o_score_list_rep_scrps ,
              score_ll=score_ll_est,
              pred_sroot=pred_score_o_sroot,
              pred_ll=pred_score_o_ll))
}



#' Title
#'
#' @param spde
#' @param n_mesh
#' @param n_rep
#' @param Q
#' @param n_outlier
#' @param outlier_val
#'
#' @return
#' @export
#'
#' @examples
repeated_inference <- function(spde,n_mesh,n_rep,Q,n_outlier=0,outlier_val=NULL,m=NULL){
  narr_rep <- rep(n_mesh,n_rep)
  times_score_rep <- rep(0,n_rep)
  times_score_rep_crps <- rep(0,n_rep)
  times_score_rep_scrps <- rep(0,n_rep)
  times_log_rep <- rep(0,n_rep)
  times_log_score_rep <- rep(0,n_rep)
  o_score_list_rep <- vector("list", n_rep)
  o_score_list_rep_crps <- vector("list", n_rep)
  o_score_list_rep_scrps <- vector("list", n_rep)
  o_log_list_rep <- vector("list", n_rep)
  o_log_score_list_rep <- vector("list", n_rep)
  times_score_rep_2 <- rep(0,n_rep)
  o_score_list_rep_2 <- vector("list", n_rep)


  for(i_n in c(1:n_rep)){

    print(paste("Iteration",i_n))
    n <- narr_rep[i_n]
    mu<- rep(0,n)
    #if(is.null(m)){
    m<-t(inla.qsample(n=10, Q = Q, mu=mu))
    #m[,1]<-4
    if(n_outlier>0){
      m[ matrix(c(1:n_outlier,sample(1:n,n_outlier)),ncol=2)]<-outlier_val #set random observation at each sample to 4.
    }
    #}


    my_obj_func_3 <- function(par,scoretype="sroot"){
      theta <- par
      mu <- rep(0,n)
      Qtheta <- inla.spde.precision(spde, theta=theta)
      return(loo_score_vectorised(m,mu,Qtheta,score=scoretype))
    }

    my_log_score_obj_func <- function(par){
      theta <- par
      mu <- rep(0,n)
      Qtheta <- inla.spde.precision(spde, theta=theta)
      return(loo_log_score(m,mu,Qtheta))
    }

    my_log_obj_func <- function(par){
      theta <- par
      mu <- rep(0,n)
      Qtheta <- inla.spde.precision(spde, theta=theta)
      return(-log_dmvn(m,mu,Qtheta))
    }

    kappa0<--2
    tau0<-1


    starttime <- Sys.time()
    o1<-optim(par=c(kappa0,tau0),my_obj_func_3,control=list(maxit=50000))
    endtime <- Sys.time()
    times_score_rep[i_n]<-difftime(endtime,starttime, units="secs")
    o_score_list_rep[[i_n]]<-o1

    starttime <- Sys.time()
    o2<-optim(par=c(kappa0,tau0),my_log_obj_func,control=list(maxit=50000))
    endtime <- Sys.time()
    times_log_rep[i_n]<-difftime(endtime,starttime, units="secs")
    o_log_list_rep[[i_n]]<-o2

    starttime <- Sys.time()
    o3<-optim(par=c(kappa0,tau0),my_log_score_obj_func,control=list(maxit=50000))
    endtime <- Sys.time()
    times_log_score_rep[i_n]<-difftime(endtime,starttime, units="secs")
    o_log_score_list_rep[[i_n]]<-o3

    starttime <- Sys.time()
    o4<-optim(par=c(kappa0,tau0),my_obj_func_3,scoretype="crps",control=list(maxit=50000))
    endtime <- Sys.time()
    times_score_rep_crps[i_n]<-difftime(endtime,starttime, units="secs")
    o_score_list_rep_crps[[i_n]]<-o4

    starttime <- Sys.time()
    o5<-optim(par=c(kappa0,tau0),my_obj_func_3,scoretype="scrps",control=list(maxit=50000))
    endtime <- Sys.time()
    times_score_rep_scrps[i_n]<-difftime(endtime,starttime, units="secs")
    o_score_list_rep_scrps[[i_n]]<-o5

  }

  return(list(times_sroot=times_score_rep,
              times_crps=times_score_rep_crps,
              times_scrps=times_score_rep_scrps,
              times_ll =times_log_rep ,
              times_slog=times_log_score_rep,
              o_sroot =o_score_list_rep ,
              o_crps =o_score_list_rep_crps ,
              o_scrps =o_score_list_rep_scrps ,
              o_ll =o_log_list_rep ,
              o_slog=o_log_score_list_rep))
}

#n_res <- 100 #set n_res to smaller value if we want to use part of results
#' @import ggplot2
#' @import dplyr
#'
#' @export
plot_results <- function(res,n_res=NULL,extended=TRUE){

  times_score_rep <- res$times_sroot
  times_log_rep <- res$times_ll
  times_log_score_rep<-res$times_slog
  o_score_list_rep<-res$o_sroot
  o_log_list_rep <- res$o_ll
  o_log_score_list_rep<- res$o_slog

  if(is.null(n_res)){
    n_res <- length(o_score_list_rep)
  }
  #extract estimated parameters
  score_par <- sapply(o_score_list_rep[1:n_res], function(o) o$par)
  log_par <- sapply(o_log_list_rep[1:n_res], function(o) o$par)
  log_score_par <- sapply(o_log_score_list_rep[1:n_res], function(o) o$par)

  #Check if all optimisations converged
  print(sum(sapply(o_score_list_rep[1:n_res], function(o) o$convergence)))
  print(sum(sapply(o_log_list_rep[1:n_res], function(o) o$convergence)))
  print(sum(sapply(o_log_score_list_rep[1:n_res], function(o) o$convergence)))

  if(extended){
  score_par_crps <- sapply(res$o_crps[1:n_res], function(o) o$par)
  score_par_scrps <- sapply(res$o_scrps[1:n_res], function(o) o$par)
  print(sum(sapply(res$o_crps[1:n_res], function(o) o$convergence)))
  print(sum(sapply(res$o_scrps[1:n_res], function(o) o$convergence)))
  #Join results in dataframe
  par_df <- data.frame(method=rep(c("Sroot","LL","Slog","CRPS","SCRPS"),each=n_res), par = Matrix::t(cbind(score_par,log_par,log_score_par,score_par_crps,score_par_scrps)), run.time=c(times_score_rep[1:n_res],times_log_rep[1:n_res],times_log_score_rep[1:n_res],res$times_crps[1:n_res],res$times_scrps[1:n_res]),i=rep(c(1:n_res),5))
  }else{

  #Join results in dataframe
  par_df <- data.frame(method=rep(c("Sroot","LL","Slog"),each=n_res), par = Matrix::t(cbind(score_par,log_par,log_score_par)), run.time=c(times_score_rep[1:n_res],times_log_rep[1:n_res],times_log_score_rep[1:n_res]),i=rep(c(1:n_res),3))
  }
  par_df_mean <- par_df %>%
    group_by(method) %>%
    summarise_all(.funs = c(mean="mean"))
  p.scatter.1 <- ggplot(par_df,aes(x=par.1,y=par.2,color=method,shape=method))+geom_point(alpha=0.75)+scale_color_brewer(palette="Dark2")+annotate("point",x=params_true[1],y=params_true[2],col="black",shape=8)

  p.hist.p1<-ggplot(par_df,aes(x=par.1,fill=method,color=method))+geom_histogram(alpha=0.5, position="identity")+geom_vline(xintercept = params_true[1])+
    geom_vline(data = par_df_mean,aes(xintercept=par.1_mean,color=method,linetype=method))+
    scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")
  p.hist.p2<-ggplot(par_df,aes(x=par.2,fill=method,color=method))+geom_histogram(alpha=0.5, position="identity")+geom_vline(xintercept = params_true[2])+
    geom_vline(data = par_df_mean,aes(xintercept=par.2_mean,color=method,linetype=method))+
    scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")

  p.box.p1 <- ggplot(par_df,aes(x=method, y=par.1,fill=method,color=method))+geom_boxplot(alpha=0.5, position="identity")+geom_hline(yintercept = params_true[1])+
    scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")
  p.box.p2 <- ggplot(par_df,aes(x=method, y=par.2,fill=method,color=method))+geom_boxplot(alpha=0.5, position="identity")+geom_hline(yintercept = params_true[2])+
    scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")

  p.time <- ggplot(par_df,aes(x=i,y=run.time,color=method,shape=method))+geom_point(alpha=0.75)+ scale_color_brewer(palette = "Dark2")
  p.time.hist <- ggplot(par_df,aes(x=run.time,color=method,fill=method))+geom_histogram(alpha=0.5,position="identity")+ scale_color_brewer(palette = "Dark2")+ scale_fill_brewer(palette = "Dark2")


  if(nrow(score_par)>2){
  p.hist.p3<-ggplot(par_df,aes(x=par.3,fill=method,color=method))+geom_histogram(alpha=0.5, position="identity")+geom_vline(xintercept = params_true[3])+
    geom_vline(data = par_df_mean,aes(xintercept=par.3_mean,color=method,linetype=method))+
    scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")
  p.box.p3 <- ggplot(par_df,aes(x=method, y=par.3,fill=method,color=method))+geom_boxplot(alpha=0.5, position="identity")+geom_hline(yintercept = params_true[3])+
    scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")
  p.scatter.2 <- ggplot(par_df,aes(x=par.2,y=par.3,color=method,shape=method))+geom_point(alpha=0.75)+scale_color_brewer(palette="Dark2")+annotate("point",x=params_true[2],y=params_true[3],col="black",shape=8)
  p.scatter.3 <- ggplot(par_df,aes(x=par.1,y=par.3,color=method,shape=method))+geom_point(alpha=0.75)+scale_color_brewer(palette="Dark2")+annotate("point",x=params_true[1],y=params_true[3],col="black",shape=8)
  return(list(p.scatter=p.scatter.1,p.scatter2=p.scatter.2,p.scatter3=p.scatter.3,p.hist.p1=p.hist.p1,p.hist.p2=p.hist.p2,p.hist.p3=p.hist.p3,p.time=p.time,p.time.hist=p.time.hist,p.box.p1=p.box.p1,p.box.p2=p.box.p2,p.box.p3=p.box.p3,df=par_df))
}else{
  return(list(p.scatter=p.scatter.1,p.hist.p1=p.hist.p1,p.hist.p2=p.hist.p2,p.time=p.time,p.time.hist=p.time.hist,p.box.p1=p.box.p1,p.box.p2=p.box.p2,df=par_df))
}

}

#' @import cowplot ggplot2
#'
#' @export
plot_grid_2 <- function(p1,p2){
  prow <- plot_grid(
    p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B"),
    hjust = -1,
    nrow = 1
  )

  # extract a legend that is laid out horizontally
  legend_b <- get_legend(
    p2 +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )

  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  return(plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1)))
}

#' @import cowplot ggplot2
#'
#' @export
plot_grid_3 <- function(p1,p2,p3){
  prow <- plot_grid(
    p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    p3 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B", "C"),
    hjust = -1,
    nrow = 1
  )

  # extract a legend that is laid out horizontally
  legend_b <- get_legend(
    p2 +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )

  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  return(plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1)))
}



plot_grid_3_nolegend <- function(p1,p2,p3){
  prow <- plot_grid(
    p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    p3 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B", "C"),
    hjust = -1,
    nrow = 1
  )


  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  return(plot_grid(prow, ncol = 1, rel_heights = c(1, .1)))
}
# plot(score_par[1,],log_par[1,])
# abline(a=0,b=1)
# plot(score_par[2,],log_par[2,])
# abline(a=0,b=1)


#' @import cowplot ggplot2
#'
#' @export
plot_grid_4 <- function(p1,p2,p3,p4){
  prow <- plot_grid(
    p4 + theme(legend.position="none"),
    p5 + theme(legend.position="none"),
    p6 + theme(legend.position="none"),
    p7 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B", "C", "D"),
    hjust = -1,
    nrow = 2
  )

  # extract a legend that is laid out horizontally
  legend_b <- get_legend(
    p4 +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )

  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  return(plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .05)))
}

plot_grid_2_nolegend <- function(p1,p2){
  prow <- plot_grid(
    p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B"),
    hjust = -1,
    nrow = 1
  )


  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  return(plot_grid(prow, ncol = 1, rel_heights = c(1, .1)))
}
