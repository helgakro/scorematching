library(INLA)
library(grid)
library(cowplot)

set.seed(111101)

sim_loc = matrix(runif(100,min=0,max=5), ncol = 2, byrow = T)
test_loc = matrix(runif(100,min=0,max=5), ncol = 2, byrow = T)
mesh_sim = inla.mesh.2d(loc = matrix(c(0,0,5,5, 0, 5, 5, 0), nrow = 4, byrow = T), max.edge=c(0.5, 1))
plot(mesh_sim)
points(sim_loc,col='red',pch=19)
mesh_sim$n
A <- inla.spde.make.A(mesh=mesh_sim,loc=sim_loc)
Atest <- inla.spde.make.A(mesh=mesh_sim,loc=test_loc)

spde = inla.spde2.matern(mesh_sim, alpha = 2)
sigma_val <- 0.1
params_true=spde$param.inla$theta.initial
params_true<-c(log(sigma_val),params_true)
print(params_true)
Q = inla.spde.precision(spde, theta=spde$param.inla$theta.initial)


library(ggplot2)
library(inlabru)
Amap.df <- data.frame(lon=c(sim_loc[,1],test_loc[,1]),lat=c(sim_loc[,2],test_loc[,2]),type=rep(c("Train","Test"),each=50))
p.A.map<-ggplot(Amap.df)+gg(mesh_sim,edge.color = "gray",
                                    int.color = "black",
                                    ext.color = "black")+geom_point(aes(x=lon,y=lat,shape=type),size=2)+scale_color_viridis_d()+xlab("longitude")+ylab("latitude")
p.A.map
ggsave('maternA_map_all.pdf',p.A.map,dpi = 1200,width = 12,height = 10,units = 'cm')

########### rep optim ##############
n_rep<-500#100
res_no_outliers_nresp <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,sigma_val=sigma_val,A=A) #no outliers 0.0002
p.res_no_outliers_nresp <- plot_results(res_no_outliers_nresp)
p.res_no_outliers_nresp$p.scatter
p.res_no_outliers_nresp$p.hist.p1
p.res_no_outliers_nresp$p.hist.p2
p.res_no_outliers_nresp$p.hist.p3
p.res_no_outliers_nresp$p.time
p.res_no_outliers_nresp$p.time.hist


p.par.hist <- plot_grid_3(p.res_no_outliers_nresp$p.hist.p1+xlab(TeX("$\\log(\\sigma)$")),
                          p.res_no_outliers_nresp$p.hist.p2+xlab(TeX("$\\log(\\kappa)$")),
                          p.res_no_outliers_nresp$p.hist.p3+xlab(TeX("$\\log(\\tau)$")))
ggsave('Ainlaparhist_sigma01.pdf',p.par.hist,dpi = 1200,width = 18,height = 8,units = 'cm')
ggsave('Ainlaparscatter_sigma01.pdf',p.res_no_outliers_nresp$p.scatter+xlab(TeX("$\\log(\\kappa)$"))+ylab(TeX("$\\log(\\tau)$")),dpi = 1200,width = 18,height = 8,units = 'cm')
ggsave('Ainlatimehist_sigma01.pdf',p.res_no_outliers_nresp$p.time.hist,dpi = 1200,width = 18,height = 8,units = 'cm')
################ score of repeated other

scores_no_outliers_nresp <- repeated_score_norm_resp(res_no_outliers_nresp,spde,mesh_sim$n,n_rep,Q,sigma_val=sigma_val,A=A) #no outliers 0.0002


score_res_nresp_df <- data.frame(val=c(sapply(scores_no_outliers_nresp$o_sroot,function(x) x[1]),sapply(scores_no_outliers_nresp$o_ll,function(x) x[1])),type=rep(c("sroot","ll"),each=n_rep))
p.score.hist1<-ggplot(score_res_nresp_df,aes(x=val,color=type,fill=type))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")

score_res_nresp_df_2 <- data.frame(val.sroot=sapply(scores_no_outliers_nresp$o_sroot,function(x) x[1]),val.ll=sapply(scores_no_outliers_nresp$o_ll,function(x) x[1]))
p.score.hist2<-ggplot(score_res_nresp_df_2,aes(x=val.ll-val.sroot,color="#1B9E77",fill="#1B9E77"))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+theme(legend.position = "none")

hist(sapply(res_no_outliers_nresp$o_sroot, function(o) o$value))
hist(sapply(scores_no_outliers_nresp$o_sroot,function(x) x[1]))

hist(sapply(res_no_outliers_nresp$o_ll, function(o) o$value))
hist(sapply(scores_no_outliers_nresp$o_ll,function(x) x[1]))

hist(sapply(scores_no_outliers_nresp$o_sroot,function(x) x[1])-sapply(scores_no_outliers_nresp$o_ll,function(x) x[1]))

################ score of LL est #####################



score_sroot <- sapply(res_no_outliers_nresp$o_sroot[1:n_rep], function(o) o$value)
score_ll <- res_no_outliers_nresp$score_ll

score_res_nresp_df_3 <- data.frame(val.sroot=score_sroot,val.ll=score_ll)
p.score.hist3<-ggplot(score_res_nresp_df_3,aes(x=val.ll-val.sroot,color="#1B9E77",fill="#1B9E77"))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+theme(legend.position = "none")

plot(score_sroot,score_ll)
abline(a=0,b=1)


plot(score_ll-score_sroot)
abline(a=0,b=0)


hist(score_ll-score_sroot)


p.score.hist <- plot_grid_2_nolegend(p.score.hist3,p.score.hist2)
p.score.hist
ggsave('Ainlascorehist_sigma01.pdf',p.score.hist,dpi = 1200,width = 18,height = 8,units = 'cm')
#####################################

score_par <- sapply(res_no_outliers_nresp$o_sroot[1:n_res], function(o) o$par)
log_score_par <- sapply(res_no_outliers_nresp$o_slog[1:n_res], function(o) o$par)
print(sum(sapply(res_no_outliers_nresp$o_sroot[1:n_res], function(o) o$convergence)))
par_df <- data.frame(method=rep(c("Sroot","Slog"),each=n_res), mu = c(score_par[1,],log_score_par[1,]), rho = c(score_par[2,],log_score_par[2,]), run.time=c(res_no_outliers_nresp$times_sroot[1:n_res],res_no_outliers_nresp$times_slog[1:n_res]),i=rep(c(1:n_res),2))
par_df_mean <- par_df %>%
  group_by(method) %>%
  summarize(meanmu = mean(mu),meanrho = mean(rho))
p.scatter <- ggplot(par_df,aes(x=mu,y=rho,color=method,shape=method))+geom_point(alpha=0.75)+scale_color_brewer(palette="Dark2")+annotate("point",x=params_true[1],y=params_true[2],col="black",shape=8)
p.scatter
p.hist.mu<-ggplot(par_df,aes(x=mu,fill=method,color=method))+geom_histogram(alpha=0.5, position="identity")+geom_vline(xintercept = params_true[1])+
  geom_vline(data = par_df_mean,aes(xintercept=meanmu,color=method,linetype=method))+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")
p.hist.rho<-ggplot(par_df,aes(x=rho,fill=method,color=method))+geom_histogram(alpha=0.5, position="identity")+geom_vline(xintercept = params_true[2])+
  geom_vline(data = par_df_mean,aes(xintercept=meanrho,color=method,linetype=method))+
  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")

p.time <- ggplot(par_df,aes(x=i,y=run.time,color=method,shape=method))+geom_point(alpha=0.75)+ scale_color_brewer(palette = "Dark2")
p.time.hist <- ggplot(par_df,aes(x=run.time,color=method,fill=method))+geom_histogram(alpha=0.5,position="identity")+ scale_color_brewer(palette = "Dark2")+ scale_fill_brewer(palette = "Dark2")


########### optim once ############


n_x <- mesh_sim$n

if(is.null(A)){
  A<-Diagonal(n_x)
}
n_y <- nrow(A)
I<-Diagonal(n_y)

n_sample <- 10
mu<- rep(0,n_x)

m<-inla.qsample(n=n_sample, Q = Q, mu=mu) #observations of latent field
m <- A%*%m+rnorm(n=n_y*n_sample,mean = 0,sd = sigma_val) #added noise
m <- Matrix::t(m)
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
  mu <- rep(0,n_x)
  #Qxy <- inla.spde.precision(spde, theta=theta) +I*sigma_val^2
  Qx <- inla.spde.precision(spde, theta=theta)
  Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A))
  Qeps  <- I/sigma_val^2
  Qtheta <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx
  muy <- A%*%mu

  #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
  score <- loo_score_vectorised(m,muy,Qtheta)
  print(score)
  return(score)
}

kappa0<-1
tau0<-0.1


o1<-optim(par=c(kappa0,tau0),my_obj_func_new_y,control=list(maxit=50000))


my_obj_func_new_y(spde$param.inla$theta.initial)

loo_score_vectorised_eps(m,mu,Q,sigma_val,A)
loo_log_score_eps(m,mu,Q,sigma_val,A)

o1$value
o1$par
spde$param.inla$theta.initial
my_obj_func_new(spde$param.inla$theta.initial)



##################### TEST set ############################

res_traintest_nresp <- repeated_inference_norm_resp(spde,mesh_sim$n,100,Q,sigma_val=sigma_val,A=A,Atest=Atest)

hist((res_traintest_nresp$pred_ll-res_traintest_nresp$pred_sroot)/res_traintest_nresp$pred_sroot)
