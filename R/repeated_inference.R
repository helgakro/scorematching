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
    Qx <- Matrix::solve(A%*%Matrix::solve(Qx,Matrix::t(A)))

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
    Qx <- Matrix::solve(A%*%Matrix::solve(Qx,Matrix::t(A)))

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
    Qx <- Matrix::solve(A%*%Matrix::solve(Qx,Matrix::t(A)))
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

    Qtheta <- Qeps-Qeps%*%A%*%Matrix::solve(Qx+Matrix::t(A)%*%Qeps%*%A,Matrix::t(A)%*%Qeps)
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

    Qtheta <- Qeps-Qeps%*%A%*%Matrix::solve(Qx+Matrix::t(A)%*%Qeps%*%A,Matrix::t(A)%*%Qeps)
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
    Qtheta <- Qeps-Qeps%*%A%*%Matrix::solve(Qx+Matrix::t(A)%*%Qeps%*%A,Matrix::t(A)%*%Qeps)
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



# TODO commented out
#
# repeated_score_norm_resp <- function(res,spde,n_mesh,n_rep,Q,n_outlier=0,outlier_val=NULL,sigma_val=1,A=NULL){
#   narr_rep <- rep(n_mesh,n_rep)
#   times_score_rep <- rep(0,n_rep)
#   times_log_rep <- rep(0,n_rep)
#   score_ll_est <- rep(0,n_rep)
#   times_log_score_rep <- rep(0,n_rep)
#   o_score_list_rep <- vector("list", n_rep)
#   o_log_list_rep <- vector("list", n_rep)
#   o_log_score_list_rep <- vector("list", n_rep)
#   times_score_rep_2 <- rep(0,n_rep)
#   o_score_list_rep_2 <- vector("list", n_rep)
#
#
#   for(i_n in c(1:n_rep)){
#
#     print(paste("Iteration",i_n))
#
#     n_x <-  n_mesh #narr_rep[i_n]
#     if(is.null(A)){
#       A<-Diagonal(n_x)
#     }
#     n_y <- nrow(A)
#     I<-Diagonal(n_y)
#     mu<- rep(0,n_x)
#     # I<-Diagonal(n)
#     # A<-Diagonal(n)
#     # m<-t(inla.qsample(n=1, Q = Q, mu=mu)) #observations of latent field
#     # m <- m+rnorm(n=n,mean = 0,sd = sigma_val) #added noise
#     n_sample <- 10
#     m<-inla.qsample(n=n_sample, Q = Q, mu=mu) #observations of latent field
#     m <- A%*%m+rnorm(n=n_y*n_sample,mean = 0,sd = sigma_val) #added noise
#     m <- Matrix::t(m)
#     #m[,1]<-4
#     if(n_outlier>0){ #add outliers
#       m[ matrix(c(1:n_outlier,sample(1:n,n_outlier)),ncol=2)]<-outlier_val #set random observation at each sample to 4.
#     }
#
#
#     my_obj_func_3_x <- function(par){
#       print(par)
#       theta <- par
#       mu <- rep(0,n)
#       #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
#       Qx <- inla.spde.precision(spde, theta=theta)
#       mux <- mu
#       #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
#       score <- loo_score_vectorised_eps(m,mux,Qx,sigma_val,A)
#       print(score)
#       return(score)
#     }
#
#     my_obj_func_3 <- function(par){
#       print(par)
#
#       sig <- exp(par[1])
#       theta <- par[-1]
#       mu <- rep(0,n_x)
#       #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
#       Qx <- inla.spde.precision(spde, theta=theta)
#       Qx <- Matrix::solve(A%*%Matrix::solve(Qx,Matrix::t(A)))
#       Qeps <- Qeps <- I/sig^2
#
#       Qtheta <- Qx-Qx%*%Matrix::solve(Qx+Qeps,Qx)
#       muy <- A%*%mu
#       #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
#       score <- loo_score_vectorised(m,muy,Qtheta)
#       print(score)
#       return(score)
#     }
#
#     # my_log_score_obj_func <- function(par){
#     #   theta <- par
#     #   mu <- rep(0,n)
#     #   Qx <- inla.spde.precision(spde, theta=theta)
#     #   mux <- mu
#     #   score <- loo_log_score_eps(m,mux,Qx,sigma_val,A)
#     #   print(score)
#     #   return(score)
#     # }
#
#     my_log_score_obj_func <- function(par){
#       print(par)
#
#       sig <- exp(par[1])
#       theta <- par[-1]
#       mu <- rep(0,n_x)
#       #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
#       Qx <- inla.spde.precision(spde, theta=theta)
#       Qx <- Matrix::solve(A%*%Matrix::solve(Qx,Matrix::t(A)))
#       Qeps <- I/sig^2
#
#       Qtheta <- Qx-Qx%*%Matrix::solve(Qx+Qeps,Qx)
#       muy <- A%*%mu
#       #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
#       score <- loo_log_score(m,muy,Qtheta)
#       print(score)
#       return(score)
#     }
#
#     #
#     # my_log_obj_func <- function(par){
#     #   theta <- par
#     #   mu <- rep(0,n)
#     #   Qxy <- inla.spde.precision(spde, theta=theta)+I*sigma_val^2
#     #   muxy <- mu
#     #   if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
#     #   return(-log_dmvn_eps(m,muxy,Qxy,sigma_val))
#     # }
#
#
#     my_log_obj_func <- function(par){
#       print(par)
#       sig <- exp(par[1])
#       theta <- par[-1]
#       mu <- rep(0,n_x)
#       #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
#       Qx <- inla.spde.precision(spde, theta=theta)
#       Qx <- Matrix::solve(A%*%Matrix::solve(Qx,Matrix::t(A)))
#       Qeps <- I/sig^2
#       Qtheta <- Qx-Qx%*%Matrix::solve(Qx+Qeps,Qx)
#       muy <- A%*%mu
#       #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
#       score <- -log_dmvn(m,muy,Qtheta)
#       print(score)
#       return(score)
#     }
#
#
#
#     o_score_list_rep[[i_n]]<-c(my_obj_func_3(res$o_sroot[[i_n]]$par),my_log_obj_func(res$o_sroot[[i_n]]$par),my_log_score_obj_func(res$o_sroot[[i_n]]$par))
#
#
#     o_log_list_rep[[i_n]]<-c(my_obj_func_3(res$o_ll[[i_n]]$par),my_log_obj_func(res$o_ll[[i_n]]$par),my_log_score_obj_func(res$o_ll[[i_n]]$par))
#
#
#     o_log_score_list_rep[[i_n]]<-c(my_obj_func_3(res$o_slog[[i_n]]$par),my_log_obj_func(res$o_slog[[i_n]]$par),my_log_score_obj_func(res$o_slog[[i_n]]$par))
#
#
#   }
#
#   return(list(o_sroot =o_score_list_rep ,
#               o_ll =o_log_list_rep ,
#               o_slog=o_log_score_list_rep))
# }


repeated_inference_norm_resp <- function(spde,n_mesh,n_rep,Q,n_outlier=0,outlier_val=NULL,sigma_val=1,A=NULL,Atest=NULL,tnu=NULL,scoretypes=c("sroot","ll","slog")){
  narr_rep <- rep(n_mesh,n_rep)
  times_score_rep <- rep(0,n_rep)
  times_log_rep <- rep(0,n_rep)
  score_ll_est <- rep(0,n_rep)
  times_score_rep_crps <- rep(0,n_rep)
  times_score_rep_scrps <- rep(0,n_rep)
  times_score_rep_rcrps <- rep(0,n_rep)
  times_log_score_rep <- rep(0,n_rep)
  o_score_list_rep <- vector("list", n_rep)
  o_log_list_rep <- vector("list", n_rep)
  o_log_score_list_rep <- vector("list", n_rep)
  o_score_list_rep_crps <- vector("list", n_rep)
  o_score_list_rep_scrps <- vector("list", n_rep)
  o_score_list_rep_rcrps <- vector("list", n_rep)
  times_score_rep_2 <- rep(0,n_rep)
  o_score_list_rep_2 <- vector("list", n_rep)

  pred_score_o_sroot <- rep(0,n_rep)
  pred_score_o_ll <- rep(0,n_rep)

  rmse_sroot <- rep(-1,n_rep)
  rmse_ll <- rep(-1,n_rep)
  rmse_slog <- rep(-1,n_rep)
  rmse_crps <- rep(-1,n_rep)
  rmse_scrps <- rep(-1,n_rep)
  rmse_rcrps <- rep(-1,n_rep)

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
      Qtheta <- Matrix::solve(A%*%Matrix::solve(Qx,Matrix::t(A))+I*sig^2)
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
      Qtheta <- Matrix::solve(A%*%Matrix::solve(Qx,Matrix::t(A))+I*sig^2)
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
      Qtheta <- Matrix::solve(A%*%Matrix::solve(Qx,Matrix::t(A))+I*sig^2)
      muy <- A%*%mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      score <- -log_dmvn(m,muy,Qtheta)
      print(score)
      return(score)
    }


    my_rmse <- function(par,A,m){
      sig <- exp(par[1])
      theta <- par[-1]
      mu <- rep(0,n_x)
      Qx <- inla.spde.precision(spde, theta=theta)
      Qtheta <- Matrix::solve(A%*%Matrix::solve(Qx,Matrix::t(A))+I*sig^2)
      muy <- A%*%mu
      score <- get_rmse(m,muy,Qtheta)
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
      times_score_rep[i_n]<--1
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
      o_log_list_rep[[i_n]]<--1

      score_ll_est[[i_n]] <- -1
    }

    if("slog"%in%scoretypes){
      starttime <- Sys.time()
    o3<-optim(par=c(sig0,kappa0,tau0),my_log_score_obj_func,A=A,m=m,control=list(maxit=50000))
    endtime <- Sys.time()
    times_log_score_rep[i_n]<-difftime(endtime,starttime, units="secs")
    o_log_score_list_rep[[i_n]]<-o3
    }else{
      times_log_score_rep[i_n]<- -1
      o_log_score_list_rep[[i_n]]<- NULL
    }




    if("crps"%in%scoretypes){
      starttime <- Sys.time()
      o4<-optim(par=c(sig0,kappa0,tau0),my_obj_func_3,A=A,m=m,scoretype="crps",control=list(maxit=50000))

      endtime <- Sys.time()
      times_score_rep_crps[i_n]<-difftime(endtime,starttime, units="secs")
      o_score_list_rep_crps[[i_n]]<-o4
    }else{
      times_score_rep_crps[i_n]<--1
      o_score_list_rep_crps[[i_n]]<-NULL
    }

    if("scrps"%in%scoretypes){
      starttime <- Sys.time()
      o5<-optim(par=c(sig0,kappa0,tau0),my_obj_func_3,A=A,m=m,scoretype="scrps",control=list(maxit=50000))

      endtime <- Sys.time()
      times_score_rep_scrps[i_n]<-difftime(endtime,starttime, units="secs")
      o_score_list_rep_scrps[[i_n]]<-o5
    }else{
      times_score_rep_scrps[i_n]<--1
      o_score_list_rep_scrps[[i_n]]<-NULL
    }

    if("rcrps"%in%scoretypes){
      starttime <- Sys.time()
      o6<-optim(par=c(sig0,kappa0,tau0),my_obj_func_3,A=A,m=m,scoretype="rcrps",control=list(maxit=50000))

      endtime <- Sys.time()
      times_score_rep_rcrps[i_n]<-difftime(endtime,starttime, units="secs")
      o_score_list_rep_rcrps[[i_n]]<-o6
    }else{
      times_score_rep_rcrps[i_n]<--1
      o_score_list_rep_rcrps[[i_n]]<-NULL
    }

    if(!is.null(Atest)){
      n_test <- nrow(Atest)
      mtest <- Atest%*%mfield+rnorm(n=n_test*n_sample,mean = 0,sd = sigma_val) #added noise
      mtest <- Matrix::t(mtest)
      if("sroot"%in%scoretypes){
        pred_score_o_sroot[[i_n]] <- my_obj_func_3(o1$par,A=Atest,m=mtest)
        rmse_sroot[[i_n]] <- my_rmse(o1$par,A=Atest,m=mtest)
      }
      if("ll"%in%scoretypes){
        pred_score_o_ll[[i_n]] <- my_obj_func_3(o2$par,A=Atest,m=mtest)
        rmse_ll[[i_n]] <- my_rmse(o2$par,A=Atest,m=mtest)
      }
      if("slog"%in%scoretypes){
        rmse_slog[[i_n]] <- my_rmse(o3$par,A=Atest,m=mtest)
      }
      if("crps"%in%scoretypes){
        rmse_crps[[i_n]] <- my_rmse(o4$par,A=Atest,m=mtest)
      }
      if("scrps"%in%scoretypes){
        rmse_scrps[[i_n]] <- my_rmse(o5$par,A=Atest,m=mtest)
      }
      if("rcrps"%in%scoretypes){
        rmse_rcrps[[i_n]] <- my_rmse(o6$par,A=Atest,m=mtest)
      }
    }

  }

  return(list(times_sroot=times_score_rep,
              times_ll =times_log_rep ,
              times_slog=times_log_score_rep,
              times_crps=times_score_rep_crps,
              times_scrps=times_score_rep_scrps,
              times_rcrps=times_score_rep_rcrps,
              o_sroot =o_score_list_rep ,
              o_ll =o_log_list_rep ,
              o_slog=o_log_score_list_rep,
              o_crps =o_score_list_rep_crps ,
              o_scrps =o_score_list_rep_scrps ,
              o_rcrps =o_score_list_rep_rcrps ,
              score_ll=score_ll_est,
              pred_sroot=pred_score_o_sroot,
              pred_ll=pred_score_o_ll,
              rmse_sroot=rmse_sroot,
              rmse_ll=rmse_ll,
              rmse_slog=rmse_slog,
              rmse_crps=rmse_crps,
              rmse_scrps=rmse_scrps,
              rmse_rcrps=rmse_rcrps))
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
    m<-t(inla.qsample(n=1, Q = Q, mu=mu)) #10
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
plot_results <- function(res,n_res=NULL,extended=TRUE,transf=FALSE){

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
  if(transf){
    score_par <- exp(score_par)
    log_par <- exp(log_par)
    log_score_par <- exp(log_score_par)
  }
  #Check if all optimisations converged
  print(sum(sapply(o_score_list_rep[1:n_res], function(o) o$convergence)))
  print(sum(sapply(o_log_list_rep[1:n_res], function(o) o$convergence)))
  print(sum(sapply(o_log_score_list_rep[1:n_res], function(o) o$convergence)))

  if(extended){
  score_par_crps <- sapply(res$o_crps[1:n_res], function(o) o$par)
  score_par_scrps <- sapply(res$o_scrps[1:n_res], function(o) o$par)
  score_par_rcrps <- sapply(res$o_rcrps[1:n_res], function(o) o$par)
  if(transf){
    score_par_crps <- exp(score_par_crps)
    score_par_scrps <- exp(score_par_scrps)
    score_par_rcrps <- exp(score_par_rcrps)
  }
  print(sum(sapply(res$o_crps[1:n_res], function(o) o$convergence)))
  print(sum(sapply(res$o_scrps[1:n_res], function(o) o$convergence)))
  print(sum(sapply(res$o_rcrps[1:n_res], function(o) o$convergence)))
  #Join results in dataframe
  par_df <- data.frame(method=rep(c("Sroot","LL","Slog","CRPS","SCRPS","rCRPS"),each=n_res), par = Matrix::t(cbind(score_par,log_par,log_score_par,score_par_crps,score_par_scrps,score_par_rcrps)), run.time=c(times_score_rep[1:n_res],times_log_rep[1:n_res],times_log_score_rep[1:n_res],res$times_crps[1:n_res],res$times_scrps[1:n_res],res$times_rcrps[1:n_res]),i=rep(c(1:n_res),6))
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


time_comparison <- function(nlist,nlist2=NULL,Qsparse=TRUE){
  n_rep <- length(nlist)
  times_score <- rep(0,n_rep)
  times_ll <- rep(0,n_rep)

  times_df <- data.frame()

  for(i_n in c(1:n_rep)){

    print(paste("Iteration",i_n))
    h<-1
    kappa <- 1
    tau <- 2

    getQ_regular <- function(kappa,tau,n_dim,h=1){
      C<-1/h*Matrix::Diagonal(n_dim,1)
      C[abs(row(C) - col(C)) == 1] <- -1/(3*h)
      G = 1/h*Matrix::Diagonal(n_dim, 2)
      G[abs(row(G) - col(G)) == 1] <- -1/h

      K <- kappa^2*C+G
      Q <- tau^2*K%*%Matrix::solve(C,K)
      Q[abs(Q)<0.01]<-0
      Q<-Matrix::Matrix(Q,sparse = TRUE)

      return(Q)
    }

    getQ_dense <- function(kappa,tau,n_dim,h=1){
      C<-1/h*Matrix::Diagonal(n_dim,1)
      C[abs(row(C) - col(C)) == 1] <- -1/(3*h)
      G = 1/h*Matrix::Diagonal(n_dim, 2)
      G[abs(row(G) - col(G)) == 1] <- -1/h

      K <- kappa^2*C+G
      Q <- tau^2*K%*%Matrix::solve(C,K)
      #Q[abs(Q)<0.01]<-0
      #Q<-Matrix::Matrix(Q,sparse = TRUE)

      return(Q)
    }

    getQ_regular_multi <- function(kappa,tau,n_dim_1,n_dim_2,h=1){
      Q<-kronecker(getQ_regular(kappa,tau,n_dim_1,h),getQ_regular(kappa,tau,n_dim_2,h))
      return(Q)
    }


    if(is.null(nlist2)){
      if(Qsparse){
        Q<-getQ_regular(kappa,tau,nlist[i_n],h) #getQ_regular(kappa,tau,n_dim,h)
      }else{
        Q<-getQ_dense(kappa,tau,nlist[i_n],h)
      }
    }else{
      Q<-getQ_regular_multi(kappa,tau,nlist[i_n],nlist2[i_n],h)
    }

    n_dim <- nrow(Q)
    print(n_dim)
    mu<- rep(0,n_dim)
    m<-t(inla.qsample(n=1, Q = Q, mu=mu))

    my_obj_func_3 <- function(par,scoretype="sroot"){
      theta <- par
      kappa <- exp(par[1])
      tau <- exp(par[2])
      mu <- rep(0,n_dim)
      Qtheta <- getQ_regular(kappa,tau,n_dim,h)
      return(loo_score_vectorised(m,mu,Qtheta,score=scoretype))
    }

    my_log_obj_func <- function(par){
      theta <- par
      kappa <- exp(par[1])
      tau <- exp(par[2])
      mu <- rep(0,n_dim)
      Qtheta <- getQ_regular(kappa,tau,n_dim,h)
      return(-log_dmvn(m,mu,Qtheta))
    }



    mbtest<-microbenchmark::microbenchmark(
    loo_score_vectorised(m,mu,Q,score="sroot"),
    log_dmvn(m,mu,Q),
    times=100)

    # mbtest<-microbenchmark::microbenchmark(
    #   my_obj_func_3(c(kappa0,tau0)),
    #   my_log_obj_func(c(kappa0,tau0)),
    #   times=1
    # )

    times_df <- rbind(times_df,data.frame(summary(mbtest,unit="ms"),n=n_dim))
  }
  return(times_df)
}
