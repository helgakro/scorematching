#Time comparison
# M1,M2,M3 in 2 dimensions
# LL and Sroot
# n=100,1000,10000



# M1
#

get_M1 <- function(n,rho=0.5){
  print(n)
  Q <- get_prec_mat_sparse(rho,n)
  n_dim <- nrow(Q)
  mu<- rep(0,n_dim)
  m<-t(inla.qsample(n=1, Q = Q, mu=mu))#bayesSurv::rMVNorm(n = 100, mean=mu,Q=Q,param="canonical")

  return(list(m=m,mu=mu,Q=Q,n_dim=n_dim))
}



# M2

#Generate m, mu, Q
#call benchmark_scoretime
get_M2_regular <- function(n){
  # if(Qsparse){
  #   Q<-getQ_regular(kappa,tau,n,h) #getQ_regular(kappa,tau,n_dim,h)
  # }else{
  #   Q<-getQ_dense(kappa,tau,n,h)
  # }
  # n_dim <- nrow(Q)
  sim_loc = matrix(runif(2*n), ncol = 2, byrow = T)
  mesh_sim = inla.mesh.2d(loc = sim_loc, max.edge=c(0.5, 1))
  plot(mesh_sim)
  points(sim_loc,col='red',pch=19)

  spde = inla.spde2.matern(mesh_sim, alpha = 2)
  params_true=spde$param.inla$theta.initial
  Q = inla.spde.precision(spde, theta=spde$param.inla$theta.initial)
  n_dim <- nrow(Q)
  mu<- rep(0,n_dim)
  m<-t(inla.qsample(n=1, Q = Q, mu=mu))
  return(list(m=m,mu=mu,Q=Q,n_dim=n_dim,spde=spde))
}



# M3

get_M3 <- function(n,sig=0.5){
  sim_loc = matrix(runif(2*n), ncol = 2, byrow = T)
  mesh_sim = inla.mesh.2d(loc = sim_loc, max.edge=c(0.5, 1))
  plot(mesh_sim)
  points(sim_loc,col='red',pch=19)

  spde = inla.spde2.matern(mesh_sim, alpha = 2)
  params_true=spde$param.inla$theta.initial
  Q = inla.spde.precision(spde, theta=spde$param.inla$theta.initial)
  n_x <- nrow(Q)
  A <- inla.spde.make.A(mesh=mesh_sim,loc=sim_loc)
  n_y <- nrow(A)
  I<-Diagonal(n_y)
  mu<- rep(0,n_x)
  # I<-Diagonal(n)
  # A<-Diagonal(n)
  # m<-t(inla.qsample(n=1, Q = Q, mu=mu)) #observations of latent field
  # m <- m+rnorm(n=n,mean = 0,sd = sigma_val) #added noise
  n_sample <- 10
  mfield<-inla.qsample(n=n_sample, Q = Q, mu=mu) #observations of latent field
  m <- A%*%mfield+rnorm(n=n_y*n_sample,mean = 0,sd = sigma_val) #added noise
  m <- Matrix::t(m)

  Qeps <- I/sig^2
  Qtheta <- Qeps-(Qeps%*%A)%*%Matrix::solve(Q+Matrix::t(A)%*%Qeps%*%A,Matrix::t(A)%*%Qeps)
  n_dim <- nrow(Qtheta)
  mu=A%*%mu
  return(list(m=m,mu=mu,Q=Qtheta,n_dim=n_dim))
}





# time comparison

time_comparison <- function(nlist,modeltype,nlist2=NULL,Qsparse=TRUE){
  n_rep <- length(nlist)
  times_df <- data.frame()

  for(i_n in c(1:n_rep)){

    if(modeltype == "m1"){
      modelsample<-get_M1(nlist[i_n])
    }
    if(modeltype=="m2"){
      modelsample<-get_M2(nlist[i_n])
    }
    if(modeltype=="m3"){
      modelsample<-get_M3(nlist[i_n])
    }


    mbtest<-microbenchmark::microbenchmark(
      loo_score_vectorised(modelsample$m,modelsample$mu,modelsample$Q,score="sroot"),
      #loo_score_vectorised(modelsample$m,modelsample$mu,modelsample$Q,score="crps"),
      #loo_score_vectorised(modelsample$m,modelsample$mu,modelsample$Q,score="scrps"),
      #loo_score_vectorised(modelsample$m,modelsample$mu,modelsample$Q,score="rcrps"),
      #loo_log_score(modelsample$m,modelsample$mu,modelsample$Q),
      log_dmvn(modelsample$m,modelsample$mu,modelsample$Q),
      times=1000)

    times_df <- rbind(times_df,data.frame(summary(mbtest,unit="ms"),n=modelsample$n_dim,modeltype=modeltype))
  }
  return(times_df)
}


run_time_comparison <- function(nlist){

t_df1<-time_comparison(nlist,"m2")
t_df2<-time_comparison(nlist,"m1")
t_df3<-time_comparison(nlist,"m3")

t_dfall<-rbind(t_df1,t_df2,t_df3)
levels(t_dfall$expr)<-c("Sr","LL")
t_dfall%>%ggplot(aes(x=n,y=mean,shape=modeltype,color=expr))+geom_point()+
  geom_smooth(data=subset(t_dfall,log(n)>5),method = "lm", se = FALSE, aes(fill=expr,linetype=expr),alpha=0.3,size=0.5)+
  scale_fill_manual(values=c("#1B9E77", "#D95F02","#1B9E77", "#D95F02","#1B9E77", "#D95F02"))+
  scale_color_manual(values=c("#1B9E77", "#D95F02","#1B9E77", "#D95F02","#1B9E77", "#D95F02"))+
#  scale_shape_manual(values=c(17,17,19,19,15,15))+
#  scale_linetype_manual(values=c("dotted","dotted","solid","solid","dashed","dashed"))+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +ylab("time (ms)")

tmptitle <- paste("run_time_comparison_",format(Sys.time(), "%Y%m%d_%H%M%S"),sep = "")
ggsave(paste(tmptitle,".pdf",sep=""),dpi = 1200,width = 12,height = 8,units = 'cm')

save(t_df1,t_df2,t_df3,file=paste(tmptitle,".Rda",sep=""))

try(
subset(t_dfall,log(n)>5) %>%
  group_by(expr,modeltype) %>%
  do({
    mod = lm(log(mean) ~ log(n), data = .)
    # print(confint(mod))
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
)
return(t_dfall)
}


