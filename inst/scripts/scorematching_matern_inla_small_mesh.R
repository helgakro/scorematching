library(INLA)
library(grid)
library(cowplot)

sim_loc = matrix(c(0,0,1,1, 0, 1, 1, 0), nrow = 4, byrow = T)
mesh_sim = inla.mesh.2d(loc = sim_loc, max.edge=c(0.5, 1))
plot(mesh_sim)
points(sim_loc,col='red',pch=19)
mesh_sim$n

spde = inla.spde2.matern(mesh_sim, alpha = 2)
params_true=spde$param.inla$theta.initial
print(params_true)
Q = inla.spde.precision(spde, theta=spde$param.inla$theta.initial)



################################# normal response model once
spde$param.inla$theta.initial
Q = inla.spde.precision(spde, theta=spde$param.inla$theta.initial)
sigma_val <- 0.001
n <- mesh_sim$n
n_sample <- 100
mu<- rep(0,n)
I<-Diagonal(n)
A<-Diagonal(n)
m<-t(inla.qsample(n=n_sample, Q = Q, mu=mu)) #observations of latent field
m <- m+rnorm(n=n*n_sample,mean = 0,sd = sigma_val) #added noise
#m[,1]<-4
# if(n_outlier>0){ #add outliers
#   m[ matrix(c(1:n_outlier,sample(1:n,n_outlier)),ncol=2)]<-outlier_val #set random observation at each sample to 4.
# }




my_obj_func_new_x <- function(par){
  print(par)
  theta <- par
  mu <- rep(0,n)
  #Qxy <- inla.spde.precision(spde, theta=theta)+I*sigma_val^2
  Qx <- inla.spde.precision(spde, theta=theta)
  mux <- mu
  #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
  #return(loo_score_vectorised_eps(m,mux,Qx,sigma_val,A))
  return(loo_log_score_eps(m,mux,Qx,sigma_val,A))
}

my_obj_func_new_y <- function(par){
  print(par)
  theta <- par
  mu <- rep(0,n)
  #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
  Qx <- inla.spde.precision(spde, theta=theta)
  Qeps <- Qeps <- I/sigma_val^2
  Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
  mux <- mu
  #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
  score <- loo_score_vectorised(m,mux,Qtheta)
  print(score)
  return(score)
}

my_log_obj_func_y <- function(par){
  print(par)
  theta <- par
  mu <- rep(0,n)
  #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
  Qx <- inla.spde.precision(spde, theta=theta)
  Qeps <- I/sigma_val^2
  Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
  mux <- mu
  #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
  score <- -log_dmvn(m,mux,Qtheta)
  print(score)
  return(score)
}

kappa0<-1
tau0<-0.1


microbenchmark::microbenchmark(o1<-optim(par=c(kappa0,tau0),my_obj_func_new_x,control=list(maxit=50000)),
o2<-optim(par=c(kappa0,tau0),my_obj_func_new_y,control=list(maxit=50000)),times=1)

o2<-optim(par=c(kappa0,tau0),my_log_obj_func_y,control=list(maxit=50000))

microbenchmark::microbenchmark(my_obj_func_new_x(c(kappa0,tau0)),my_obj_func_new_y(c(kappa0,tau0)),times=5)


loo_score_vectorised_eps(m,mu,Q,sigma_val,A)
loo_log_score_eps(m,mu,Q,sigma_val,A)

o1$value
o1$par
spde$param.inla$theta.initial
my_obj_func_new(spde$param.inla$theta.initial)

my_log_obj_func_y(spde$param.inla$theta.initial)

##################### Test S-M

A <- Diagonal(nrow(Q))
Qeps <- Diagonal(nrow(Q))/sigma_val^2
Qxy <- Q + t(A)%*%Qeps%*%A
invQxy <- solve(Qxy)
Qeps <- Qeps[-1,-1]
invQe <- solve(Qeps)

a1<-inv_sherman_morrison(invQxy,A[1,]/sigma_val,-A[1,]/sigma_val)
a2<-solve(Q + t(A[-1,])%*%Qeps%*%A[-1,])
min(a1-a2)


##################################
n_rep<-100
res_no_outliers_nresp_small <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,sigma_val=0.001)
n_res=n_rep
score_par <- sapply(res_no_outliers_nresp_small$o_sroot[1:n_res], function(o) o$par)
log_score_par <- sapply(res_no_outliers_nresp_small$o_slog[1:n_res], function(o) o$par)
print(sum(sapply(res_no_outliers_nresp_small$o_sroot[1:n_res], function(o) o$convergence)))

par_df <- data.frame(method=rep(c("Sroot","Slog"),each=n_res), par = t(cbind(score_par,log_score_par)), run.time=c(res_no_outliers_nresp_small$times_sroot[1:n_res],res_no_outliers_nresp_small$times_slog[1:n_res]),i=rep(c(1:n_res),2))
par_df_mean <- par_df %>%
  group_by(method) %>%
  summarize(meanmu = mean(mu),meanrho = mean(rho))

par_df <- data.frame(method=rep(c("Sroot","Slog"),each=n_res), mu = c(score_par[1,],log_score_par[1,]), rho = c(score_par[2,],log_score_par[2,]), run.time=c(res_no_outliers_nresp_small$times_sroot[1:n_res],res_no_outliers_nresp_small$times_slog[1:n_res]),i=rep(c(1:n_res),2))
par_df_mean <- par_df %>%
  group_by(method) %>%
  summarize(meanmu = mean(mu),meanrho = mean(rho))
p.scatter <- ggplot(par_df,aes(x=mu,y=rho,color=method,shape=method))+geom_point(alpha=0.75)+scale_color_brewer(palette="Dark2")+annotate("point",x=params_true[1],y=params_true[2],col="black",shape=8)
p.scatter
op.hist.mu<-ggplot(par_df,aes(x=mu,fill=method,color=method))+geom_histogram(alpha=0.5, position="identity")+geom_vline(xintercept = params_true[1])+
  geom_vline(data = par_df_mean,aes(xintercept=meanmu,color=method,linetype=method))+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")
p.hist.rho<-ggplot(par_df,aes(x=rho,fill=method,color=method))+geom_histogram(alpha=0.5, position="identity")+geom_vline(xintercept = params_true[2])+
  geom_vline(data = par_df_mean,aes(xintercept=meanrho,color=method,linetype=method))+
  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")

p.time <- ggplot(par_df,aes(x=i,y=run.time,color=method,shape=method))+geom_point(alpha=0.75)+ scale_color_brewer(palette = "Dark2")
p.time.hist <- ggplot(par_df,aes(x=run.time,color=method,fill=method))+geom_histogram(alpha=0.5,position="identity")+ scale_color_brewer(palette = "Dark2")+ scale_fill_brewer(palette = "Dark2")



res_no_outliers_nresp_small_2 <- repeated_inference_norm_resp(spde,mesh_sim$n,5,Q,sigma_val=0.01)
p.res_no_outliers_nresp_small<- plot_results(res_no_outliers_nresp_small)
p.res_no_outliers_nresp_small$p.scatter
p.res_no_outliers_nresp_small$p.hist.mu
p.res_no_outliers_nresp_small$p.hist.rho
p.res_no_outliers_nresp_small$p.time
p.res_no_outliers_nresp_small$p.time.hist

