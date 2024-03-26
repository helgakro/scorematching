repeated_inference_norm_resp <- function(n_mesh,n_rep,Q,n_outlier=0,outlier_val=NULL,sigma_val=1){
  narr_rep <- rep(n_mesh,n_rep)
  times_score_rep <- rep(0,n_rep)
  times_log_rep <- rep(0,n_rep)
  times_log_score_rep <- rep(0,n_rep)
  o_score_list_rep <- vector("list", n_rep)
  o_log_list_rep <- vector("list", n_rep)
  o_log_score_list_rep <- vector("list", n_rep)
  times_score_rep_2 <- rep(0,n_rep)
  o_score_list_rep_2 <- vector("list", n_rep)


  for(i_n in c(1:n_rep)){

    print(paste("Iteration",i_n))
    n <- narr_rep[i_n]
    mu<- rep(0,n)
    I<-Diagonal(n)
    m<-t(inla.qsample(n=1, Q = Q, mu=mu)) #observations of latent field
    m <- m+rnorm(n=n,mean = 0,sd = sigma_val) #added noise
    #m[,1]<-4
    if(n_outlier>0){ #add outliers
      m[ matrix(c(1:n_outlier,sample(1:n,n_outlier)),ncol=2)]<-outlier_val #set random observation at each sample to 4.
    }


    my_obj_func_3 <- function(par){
      theta <- par
      mu <- rep(0,n)
      Qxy <- inla.spde.precision(spde, theta=theta)+I*sigma_val^2
      mux <- mu
      #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      return(loo_score_vectorised_eps(m,mux,Qxy,sigma_val))
    }

    my_log_score_obj_func <- function(par){
      theta <- par
      mu <- rep(0,n)
      Qxy <- inla.spde.precision(spde, theta=theta)+I*sigma_val^2
      muxy <- mu
      if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      return(loo_log_score_eps(m,muxy,Qxy,sigma_val))
    }

    my_log_obj_func <- function(par){
      theta <- par
      mu <- rep(0,n)
      Qxy <- inla.spde.precision(spde, theta=theta)+I*sigma_val^2
      muxy <- mu
      if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
      return(-log_dmvn_eps(m,muxy,Qxy,sigma_val))
    }

    kappa0<-1
    tau0<-0.1


    starttime <- Sys.time()
    o1<-optim(par=c(kappa0,tau0),my_obj_func_3,control=list(maxit=50000))
    endtime <- Sys.time()
    times_score_rep[i_n]<-difftime(endtime,starttime, units="secs")
    o_score_list_rep[[i_n]]<-o1

    starttime <- Sys.time()
    #o2<-optim(par=c(kappa0,tau0),my_log_obj_func,control=list(maxit=50000))
    endtime <- Sys.time()
    times_log_rep[i_n]<-difftime(endtime,starttime, units="secs")
    #o_log_list_rep[[i_n]]<-o2

    starttime <- Sys.time()
    #o3<-optim(par=c(kappa0,tau0),my_log_score_obj_func,control=list(maxit=50000))
    endtime <- Sys.time()
    times_log_score_rep[i_n]<-difftime(endtime,starttime, units="secs")
    #o_log_score_list_rep[[i_n]]<-o3


  }

  return(list(times_sroot=times_score_rep,
              times_ll =times_log_rep ,
              times_slog=times_log_score_rep,
              o_sroot =o_score_list_rep ,
              o_ll =o_log_list_rep ,
              o_slog=o_log_score_list_rep))
}



repeated_inference <- function(n_mesh,n_rep,Q,n_outlier=0,outlier_val=NULL){
  narr_rep <- rep(n_mesh,n_rep)
  times_score_rep <- rep(0,n_rep)
  times_log_rep <- rep(0,n_rep)
  times_log_score_rep <- rep(0,n_rep)
  o_score_list_rep <- vector("list", n_rep)
  o_log_list_rep <- vector("list", n_rep)
  o_log_score_list_rep <- vector("list", n_rep)
  times_score_rep_2 <- rep(0,n_rep)
  o_score_list_rep_2 <- vector("list", n_rep)


  for(i_n in c(1:n_rep)){

    print(paste("Iteration",i_n))
    n <- narr_rep[i_n]
    mu<- rep(0,n)
    m<-t(inla.qsample(n=10, Q = Q, mu=mu))
    #m[,1]<-4
    if(n_outlier>0){
      m[ matrix(c(1:n_outlier,sample(1:n,n_outlier)),ncol=2)]<-outlier_val #set random observation at each sample to 4.
    }


    my_obj_func_3 <- function(par){
      theta <- par
      mu <- rep(0,n)
      return(loo_score_vectorised(m,mu,inla.spde.precision(spde, theta=theta)))
    }

    my_log_score_obj_func <- function(par){
      theta <- par
      mu <- rep(0,n)
      return(loo_log_score(m,mu,inla.spde.precision(spde, theta=theta)))
    }

    my_log_obj_func <- function(par){
      theta <- par
      mu <- rep(0,n)
      return(-log_dmvn(m,mu,inla.spde.precision(spde, theta=theta)))
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


  }

  return(list(times_sroot=times_score_rep,
              times_ll =times_log_rep ,
              times_slog=times_log_score_rep,
              o_sroot =o_score_list_rep ,
              o_ll =o_log_list_rep ,
              o_slog=o_log_score_list_rep))
}

#n_res <- 100 #set n_res to smaller value if we want to use part of results

plot_results <- function(res,n_res=NULL){

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

  #Join results in dataframe
  par_df <- data.frame(method=rep(c("Sroot","LL","Slog"),each=n_res), mu = c(score_par[1,],log_par[1,],log_score_par[1,]), rho = c(score_par[2,],log_par[2,],log_score_par[2,]), run.time=c(times_score_rep[1:n_res],times_log_rep[1:n_res],times_log_score_rep[1:n_res]),i=rep(c(1:n_res),3))
  par_df_mean <- par_df %>%
    group_by(method) %>%
    summarize(meanmu = mean(mu),meanrho = mean(rho))
  p.scatter <- ggplot(par_df,aes(x=mu,y=rho,color=method,shape=method))+geom_point(alpha=0.75)+scale_color_brewer(palette="Dark2")+annotate("point",x=params_true[1],y=params_true[2],col="black",shape=8)

  p.hist.mu<-ggplot(par_df,aes(x=mu,fill=method,color=method))+geom_histogram(alpha=0.5, position="identity")+geom_vline(xintercept = params_true[1])+
    geom_vline(data = par_df_mean,aes(xintercept=meanmu,color=method,linetype=method))+
    scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")
  p.hist.rho<-ggplot(par_df,aes(x=rho,fill=method,color=method))+geom_histogram(alpha=0.5, position="identity")+geom_vline(xintercept = params_true[2])+
    geom_vline(data = par_df_mean,aes(xintercept=meanrho,color=method,linetype=method))+
    scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")

  p.time <- ggplot(par_df,aes(x=i,y=run.time,color=method,shape=method))+geom_point(alpha=0.75)+ scale_color_brewer(palette = "Dark2")
  p.time.hist <- ggplot(par_df,aes(x=run.time,color=method,fill=method))+geom_histogram(alpha=0.5,position="identity")+ scale_color_brewer(palette = "Dark2")+ scale_fill_brewer(palette = "Dark2")


  return(list(p.scatter=p.scatter,p.hist.mu=p.hist.mu,p.hist.rho=p.hist.rho,p.time=p.time,p.time.hist=p.time.hist))

}


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
# plot(score_par[1,],log_par[1,])
# abline(a=0,b=1)
# plot(score_par[2,],log_par[2,])
# abline(a=0,b=1)

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
