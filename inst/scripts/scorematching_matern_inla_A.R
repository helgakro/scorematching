library(INLA)
library(grid)
library(cowplot)
library(latex2exp)

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
p.score.hist2<-ggplot(score_res_nresp_df_2,aes(x=(val.ll-val.sroot)/val.sroot,color="#1B9E77",fill="#1B9E77"))+geom_histogram(alpha=0.5, position="identity")+
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
p.score.hist3<-ggplot(score_res_nresp_df_3,aes(x=(val.ll-val.sroot)/val.sroot,color="#1B9E77",fill="#1B9E77"))+geom_histogram(alpha=0.5, position="identity")+
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


score_res_nresp_df_4 <- data.frame(val.sroot=res_traintest_nresp$pred_sroot,val.ll=res_traintest_nresp$pred_ll)
p.score.hist4<-ggplot(score_res_nresp_df_4,aes(x=(val.ll-val.sroot)/val.sroot,color="#1B9E77",fill="#1B9E77"))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+theme(legend.position = "none")


p.score.hist.multi<-plot_grid_3_nolegend(p.score.hist3+xlim(c(-0.015,0.015)),p.score.hist4+xlim(c(-0.015,0.015)),p.score.hist2+xlim(c(-0.015,0.015)))
p.score.hist.multi
ggsave('Ainlascorehist_sigma01_traintest.pdf',p.score.hist.multi,dpi = 1200,width = 18,height = 8,units = 'cm')



####################### t-response ############################

n_rep<-500#100
res_no_outliers_tresp <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,sigma_val=sigma_val,A=A,tnu=2,scoretypes=c("sroot","ll","slog","crps","scrps")) #no outliers 0.0002
p.res_no_outliers_tresp <- plot_results(res_no_outliers_tresp)
p.res_no_outliers_tresp$p.scatter
p.res_no_outliers_tresp$p.hist.p1
p.res_no_outliers_tresp$p.hist.p2
p.res_no_outliers_tresp$p.hist.p3
p.res_no_outliers_tresp$p.time
p.res_no_outliers_tresp$p.time.hist


p.par.hist.t <- plot_grid_3(p.res_no_outliers_tresp$p.hist.p1+xlab(TeX("$\\log(\\sigma)$")),
                          p.res_no_outliers_tresp$p.hist.p2+xlab(TeX("$\\log(\\kappa)$")),
                          p.res_no_outliers_tresp$p.hist.p3+xlab(TeX("$\\log(\\tau)$")))
ggsave('Ainlaparhist_t_500_nu2.pdf',p.par.hist.t,dpi = 1200,width = 18,height = 8,units = 'cm')
ggsave('Ainlaparscatter_t_500_nu2.pdf',p.res_no_outliers_tresp$p.scatter+xlab(TeX("$\\log(\\kappa)$"))+ylab(TeX("$\\log(\\tau)$")),dpi = 1200,width = 18,height = 8,units = 'cm')
ggsave('Ainlatimehist_t_500_nu2.pdf',p.res_no_outliers_tresp$p.time.hist,dpi = 1200,width = 18,height = 8,units = 'cm')


tscore_sroot <- sapply(res_no_outliers_tresp$o_sroot[1:n_rep], function(o) o$value)
tscore_ll <- res_no_outliers_tresp$score_ll

score_res_tresp_df_3 <- data.frame(val.sroot=tscore_sroot,val.ll=tscore_ll)
p.t.score.hist3<-ggplot(score_res_tresp_df_3,aes(x=(val.ll-val.sroot)/val.sroot,color="#1B9E77",fill="#1B9E77"))+geom_histogram(alpha=0.5, position="identity",binwidth = 0.005)+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+theme(legend.position = "none")
p.t.score.hist3
ggsave('Ainlascorehist_t_500_nu2_binw.pdf',p.t.score.hist3,dpi = 1200,width = 18,height = 8,units = 'cm')
############# full response
res_no_outliers_fresp <- repeated_inference_norm_resp(spde,mesh_sim$n,100,Q,sigma_val=sigma_val,A=A,tnu=0) #no outliers 0.0002
p.res_no_outliers_fresp <- plot_results(res_no_outliers_fresp)
p.res_no_outliers_fresp$p.scatter
p.res_no_outliers_fresp$p.hist.p1
p.res_no_outliers_fresp$p.hist.p2
p.res_no_outliers_fresp$p.hist.p3
p.res_no_outliers_fresp$p.time
p.res_no_outliers_fresp$p.time.hist

fscore_sroot <- sapply(res_no_outliers_fresp$o_sroot[1:100], function(o) o$value)
fscore_ll <- res_no_outliers_fresp$score_ll

score_res_fresp_df_3 <- data.frame(val.sroot=fscore_sroot,val.ll=fscore_ll)
p.f.score.hist3<-ggplot(score_res_fresp_df_3,aes(x=(val.ll-val.sroot)/val.sroot,color="#1B9E77",fill="#1B9E77"))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+theme(legend.position = "none")
p.f.score.hist3





########### rep optim outlier ##############

set.seed(20240514)
n_rep<-100#100
res_10_outliers_nresp <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 10, outlier_val = 5,sigma_val=0.5,A=A,Atest=Atest,scoretypes=c("sroot","ll","slog","crps","scrps") ) #no outliers 0.0002
params_true=spde$param.inla$theta.initial
params_true<-c(log(0.5),params_true)
p.res_10_outliers_nresp <- plot_results(res_10_outliers_nresp)
p.res_10_outliers_nresp$p.scatter
p.res_10_outliers_nresp$p.scatter2
p.res_10_outliers_nresp$p.scatter3
p.res_10_outliers_nresp$p.hist.p1
p.res_10_outliers_nresp$p.hist.p2
p.res_10_outliers_nresp$p.hist.p3
p.res_10_outliers_nresp$p.time
p.res_10_outliers_nresp$p.time.hist

hist(as.vector(A%*%inla.qsample(n=100, Q = Q, mu=rep(0,nrow(Q))))) #observations of latent field
hist(as.vector(A%*%inla.qsample(n=100, Q = Q, mu=rep(0,nrow(Q))))+rnorm(n=nrow(A)*100,mean = 0,sd = sigma_val)) #added noise


p.res_10_outliers_nresp$p.box.p1
p.res_10_outliers_nresp$p.box.p2
p.res_10_outliers_nresp$p.box.p3
p.par.box_10_outliers <- plot_grid_3(p.res_10_outliers_nresp$p.box.p1+ylab(TeX("$\\log(\\sigma)$")),
                          p.res_10_outliers_nresp$p.box.p2+ylab(TeX("$\\log(\\kappa)$")),
                          p.res_10_outliers_nresp$p.box.p3+ylab(TeX("$\\log(\\tau)$")))
ggsave('Ainlaparbox_sigma05_10outliers.pdf',p.par.box_10_outliers,dpi = 1200,width = 18,height = 8,units = 'cm')


data.frame(val.sroot=res_10_outliers_nresp$pred_sroot,val.ll=res_10_outliers_nresp$pred_ll,outlier="yes")%>%ggplot(aes(x=(val.ll-val.sroot)/val.sroot,color=outlier,fill=outlier))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")

#----------------
set.seed(20240514)
n_rep<-100#100
res_0_outliers_nresp <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 0, outlier_val = 5,sigma_val=0.5,A=A,Atest=Atest,scoretypes=c("sroot","ll","slog","crps","scrps") ) #no outliers 0.0002
params_true=spde$param.inla$theta.initial
params_true<-c(log(0.5),params_true)
p.res_0_outliers_nresp <- plot_results(res_0_outliers_nresp)
p.res_0_outliers_nresp$p.scatter
p.res_0_outliers_nresp$p.scatter2
p.res_0_outliers_nresp$p.scatter3
p.res_0_outliers_nresp$p.hist.p1
p.res_0_outliers_nresp$p.hist.p2
p.res_0_outliers_nresp$p.hist.p3
p.res_0_outliers_nresp$p.time
p.res_0_outliers_nresp$p.time.hist

p.res_0_outliers_nresp$p.box.p1
p.res_0_outliers_nresp$p.box.p2
p.res_0_outliers_nresp$p.box.p3

hist(as.vector(A%*%inla.qsample(n=100, Q = Q, mu=rep(0,nrow(Q))))) #observations of latent field
hist(as.vector(A%*%inla.qsample(n=100, Q = Q, mu=rep(0,nrow(Q))))+rnorm(n=nrow(A)*100,mean = 0,sd = sigma_val)) #added noise
p.par.box_0_outliers <- plot_grid_3(p.res_0_outliers_nresp$p.box.p1+ylab(TeX("$\\log(\\sigma)$")),
                                     p.res_0_outliers_nresp$p.box.p2+ylab(TeX("$\\log(\\kappa)$")),
                                     p.res_0_outliers_nresp$p.box.p3+ylab(TeX("$\\log(\\tau)$")))
ggsave('Ainlaparbox_sigma05_0outliers.pdf',p.par.box_0_outliers,dpi = 1200,width = 18,height = 8,units = 'cm')




#---- plot boxplots -------------

df0<-p.res_0_outliers_nresp_xy$df
df10<-p.res_10_outliers_nresp_xy$df

df0$outlier <- "no"
df10$outlier <- "yes"
df.all <- rbind(df0,df10)

df.all$method <- factor(df.all$method, levels=c('LL', 'Slog', 'SCRPS', 'Sroot','CRPS'))

p.par.box_outliers<-plot_grid_3(ggplot(df.all,aes(x=outlier, y=par.1,fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = params_true[1])+
  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\log(\\sigma)$")),
ggplot(df.all,aes(x=outlier, y=par.2,fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = params_true[2])+
  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\log(\\kappa)$")),
ggplot(df.all,aes(x=outlier, y=par.3,fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = params_true[3])+
  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\log(\\tau)$")))
ggsave('Ainlaparbox_sigma05_outliers_val5_nonstationary_xy.pdf',p.par.box_outliers,dpi = 1200,width = 18,height = 8,units = 'cm')




score_res_nresp_outliers_df <- rbind(data.frame(val.sroot=res_0_outliers_nresp_xy$pred_sroot,val.ll=res_0_outliers_nresp_xy$pred_ll,outlier="no"),data.frame(val.sroot=res_10_outliers_nresp_xy$pred_sroot,val.ll=res_10_outliers_nresp_xy$pred_ll,outlier="yes"))
p.outliers.score.test<-ggplot(score_res_nresp_outliers_df,aes(x=(val.ll-val.sroot)/val.sroot,color=outlier,fill=outlier))+geom_histogram(alpha=0.5, position="identity",binwidth = 0.003)+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")
ggsave('Ainlascoretest_sigma05_outliers_val5_nonstationary_xy.pdf',p.outliers.score.test,dpi = 1200,width = 18,height = 8,units = 'cm')

################33 Non-stationarity

n_rep<-100#100
res_0_outliers_nresp <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 0, outlier_val = 5,sigma_val=0.5,A=t(t(A)*mesh_sim$loc[,1]),Atest=t(t(Atest)*mesh_sim$loc[,1]),scoretypes=c("sroot","ll","slog","crps","scrps") ) #no outliers 0.0002
p.res_0_outliers_nresp <- plot_results(res_0_outliers_nresp)

p.res_0_outliers_nresp$p.box.p1
p.res_0_outliers_nresp$p.box.p2
p.res_0_outliers_nresp$p.box.p3


n_rep<-100#100
res_10_outliers_nresp <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 10, outlier_val = 5,sigma_val=0.5,A=t(t(A)*mesh_sim$loc[,1]),Atest=t(t(Atest)*mesh_sim$loc[,1]),scoretypes=c("sroot","ll","slog","crps","scrps") ) #no outliers 0.0002
#res_10_outliers_nresp <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 10, outlier_val = 5,sigma_val=0.5,A=t(t(A)*mesh_sim$loc[,1]),Atest=t(t(Atest)*mesh_sim$loc[,1]),scoretypes=c("sroot","ll","slog","crps","scrps") ) #no outliers 0.0002

p.res_10_outliers_nresp <- plot_results(res_10_outliers_nresp)

p.res_10_outliers_nresp$p.box.p1
p.res_10_outliers_nresp$p.box.p2
p.res_10_outliers_nresp$p.box.p3


###xy spatial

n_rep<-100#100
res_10_outliers_nresp_xy <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 10, outlier_val = 5,sigma_val=0.5,A=t(t(A)*sqrt(abs(mesh_sim$loc[,1]*mesh_sim$loc[,2]))),Atest=t(t(Atest)*mesh_sim$loc[,1]),scoretypes=c("sroot","ll","slog","crps","scrps") ) #no outliers 0.0002

p.res_10_outliers_nresp_xy <- plot_results(res_10_outliers_nresp_xy)

p.res_10_outliers_nresp_xy$p.box.p1
p.res_10_outliers_nresp_xy$p.box.p2
p.res_10_outliers_nresp_xy$p.box.p3


n_rep<-100#100
res_0_outliers_nresp_xy <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 0, outlier_val = 5,sigma_val=0.5,A=t(t(A)*sqrt(abs(mesh_sim$loc[,1]*mesh_sim$loc[,2]))),Atest=t(t(Atest)*mesh_sim$loc[,1]),scoretypes=c("sroot","ll","slog","crps","scrps") ) #no outliers 0.0002

p.res_0_outliers_nresp_xy <- plot_results(res_0_outliers_nresp_xy)

p.res_0_outliers_nresp_xy$p.box.p1
p.res_0_outliers_nresp_xy$p.box.p2
p.res_0_outliers_nresp_xy$p.box.p3



hist(as.vector(A%*%inla.qsample(n=1, Q = Q, mu=rep(0,nrow(Q)))+rnorm(n=nrow(A)*1,mean = 0,sd = 0.5)))
library(ggplot2)
library(inlabru)


Amap_ns.df <- rbind(data.frame(lon=sim_loc[,1],lat=sim_loc[,2],rain=as.vector((A*(mesh_sim$loc[,1]))%*%testfield+(sim_loc[,1])*rnorm(n=nrow(A)*1,mean = 0,sd = 0.5)),type="Train"),
                    data.frame(lon=test_loc[,1],lat=test_loc[,2],rain=as.vector((Atest*(mesh_sim$loc[,1]))%*%testfield+(test_loc[,1])*rnorm(n=nrow(A)*1,mean = 0,sd = 0.5)),type="Test"))


Amap_ns.df <- rbind(data.frame(lon=mesh_sim$loc[,1],lat=mesh_sim$loc[,2],rain=as.vector(inla.qsample(n=1, Q = Q, mu=rep(0,nrow(Q)))),type="Train"))

Amap_ns.df <- rbind(data.frame(lon=mesh_sim$loc[,1],lat=mesh_sim$loc[,2],rain=as.vector((mesh_sim$loc[,1])*inla.qsample(n=1, Q = Q, mu=rep(0,nrow(Q)))),type="Train"))



Amap_ns.df <- rbind(data.frame(lon=sim_loc[,1],lat=sim_loc[,2],rain=as.vector(((sim_loc[,1]*A)%*%inla.qsample(n=1, Q = Q, mu=rep(0,nrow(Q))))),type="Train"),
                    data.frame(lon=test_loc[,1],lat=test_loc[,2],rain=as.vector(((test_loc[,1]*Atest)%*%inla.qsample(n=1, Q = Q, mu=rep(0,nrow(Q))))),type="Test"))
Amap_ns.df <- rbind(data.frame(lon=sim_loc[,1],lat=sim_loc[,2],rain=as.vector((sim_loc[,1])*(A%*%inla.qsample(n=1, Q = Q, mu=rep(0,nrow(Q))))),type="Train"),
                    data.frame(lon=test_loc[,1],lat=test_loc[,2],rain=as.vector((test_loc[,1])*(Atest%*%inla.qsample(n=1, Q = Q, mu=rep(0,nrow(Q))))),type="Test"))


testfield <- inla.qsample(n=1, Q = Q, mu=rep(0,nrow(Q)))
map_ns.df <- data.frame(lon=mesh_sim$loc[,1],lat=mesh_sim$loc[,2],
                        field=as.vector(testfield),
                        scaledfieldx = mesh_sim$loc[,1]*as.vector(testfield),
                        scaledfield = sqrt(abs(mesh_sim$loc[,1]*mesh_sim$loc[,2]))*as.vector(testfield))
nugget <- rnorm(n=nrow(Q),mean = 0,sd = 0.5)
scalednuggetx <- abs(mesh_sim$loc[,1])*nugget
scalednugget <- sqrt(abs(mesh_sim$loc[,1]*mesh_sim$loc[,2]))*nugget
ggplot(map_ns.df)+gg(mesh_sim,edge.color = "gray",
                      int.color = "black",
                      ext.color = "black")+geom_point(aes(x=lon,y=lat,color=field),size=2)+xlab("longitude")+
  ylab("latitude")+scale_color_distiller(palette = "Spectral",limits=c(-max(abs(map_ns.df$field)),max(abs(map_ns.df$field))))
ggplot(map_ns.df)+gg(mesh_sim,edge.color = "gray",
                     int.color = "black",
                     ext.color = "black")+geom_point(aes(x=lon,y=lat,color=field+nugget),size=2)+xlab("longitude")+
  ylab("latitude")+scale_color_distiller(palette = "Spectral",limits=c(-max(abs(map_ns.df$field)),max(abs(map_ns.df$field))))
ggplot(map_ns.df)+gg(mesh_sim,edge.color = "gray",
                     int.color = "black",
                     ext.color = "black")+geom_point(aes(x=lon,y=lat,color=field+scalednugget),size=2)+xlab("longitude")+
  ylab("latitude")+scale_color_distiller(palette = "Spectral",limits=c(-max(abs(map_ns.df$field)),max(abs(map_ns.df$field))))

#x-directional non-stationary field
ggplot(map_ns.df)+gg(mesh_sim,edge.color = "gray",
                     int.color = "black",
                     ext.color = "black")+geom_point(aes(x=lon,y=lat,color=scaledfieldx),size=2)+xlab("longitude")+
  ylab("latitude")+scale_color_distiller(palette = "Spectral",limits=c(-max(abs(map_ns.df$scaledfield)),max(abs(map_ns.df$scaledfield))))
ggplot(map_ns.df)+gg(mesh_sim,edge.color = "gray",
                     int.color = "black",
                     ext.color = "black")+geom_point(aes(x=lon,y=lat,color=scaledfieldx+nugget),size=2)+xlab("longitude")+
  ylab("latitude")+scale_color_distiller(palette = "Spectral",limits=c(-max(abs(map_ns.df$scaledfield+nugget)),max(abs(map_ns.df$scaledfield+nugget))))
ggplot(map_ns.df)+gg(mesh_sim,edge.color = "gray",
                     int.color = "black",
                     ext.color = "black")+geom_point(aes(x=lon,y=lat,color=scaledfieldx+scalednuggetx),size=2)+xlab("longitude")+
  ylab("latitude")+scale_color_distiller(palette = "Spectral",limits=c(-max(abs(map_ns.df$scaledfield+scalednugget)),max(abs(map_ns.df$scaledfield+scalednugget))))


#two-directional non-stationary field
ggplot(map_ns.df)+gg(mesh_sim,edge.color = "gray",
                     int.color = "black",
                     ext.color = "black")+geom_point(aes(x=lon,y=lat,color=scaledfield),size=2)+xlab("longitude")+
  ylab("latitude")+scale_color_distiller(palette = "Spectral",limits=c(-max(abs(map_ns.df$scaledfield)),max(abs(map_ns.df$scaledfield))))
ggplot(map_ns.df)+gg(mesh_sim,edge.color = "gray",
                     int.color = "black",
                     ext.color = "black")+geom_point(aes(x=lon,y=lat,color=scaledfield+nugget),size=2)+xlab("longitude")+
  ylab("latitude")+scale_color_distiller(palette = "Spectral",limits=c(-max(abs(map_ns.df$scaledfield+nugget)),max(abs(map_ns.df$scaledfield+nugget))))
ggplot(map_ns.df)+gg(mesh_sim,edge.color = "gray",
                     int.color = "black",
                     ext.color = "black")+geom_point(aes(x=lon,y=lat,color=scaledfield+scalednugget),size=2)+xlab("longitude")+
  ylab("latitude")+scale_color_distiller(palette = "Spectral",limits=c(-max(abs(map_ns.df$scaledfield+scalednugget)),max(abs(map_ns.df$scaledfield+scalednugget))))






Amap_ns.df <- rbind(data.frame(lon=sim_loc[,1],lat=sim_loc[,2],rain=as.vector((A*(5*mesh_sim$loc[,1]))%*%testfield+(sim_loc[,1])*rnorm(n=nrow(A)*1,mean = 0,sd = 0.5)),type="Train"),
                    data.frame(lon=test_loc[,1],lat=test_loc[,2],rain=as.vector((Atest*(5*mesh_sim$loc[,1]))%*%testfield+(test_loc[,1])*rnorm(n=nrow(A)*1,mean = 0,sd = 0.5)),type="Test"))

Amap_ns.df <- rbind(data.frame(lon=sim_loc[,1],lat=sim_loc[,2],rain=as.vector((A*sqrt(abs(mesh_sim$loc[,1]*mesh_sim$loc[,2])))%*%testfield),type="Train"),
                    data.frame(lon=test_loc[,1],lat=test_loc[,2],rain=as.vector((Atest*sqrt(abs(mesh_sim$loc[,1]*mesh_sim$loc[,2])))%*%testfield),type="Test"))

Amap_ns.df <- rbind(data.frame(lon=sim_loc[,1],lat=sim_loc[,2],rain=as.vector(A%*%(sqrt(abs(mesh_sim$loc[,1]*mesh_sim$loc[,2]))*testfield)),type="Train"),
                    data.frame(lon=test_loc[,1],lat=test_loc[,2],rain=as.vector(Atest%*%(sqrt(abs(mesh_sim$loc[,1]*mesh_sim$loc[,2]))*testfield)),type="Test"))

Amap_ns.df <- rbind(data.frame(lon=sim_loc[,1],lat=sim_loc[,2],rain=as.vector((t(t(A)*sqrt(abs(mesh_sim$loc[,1]*mesh_sim$loc[,2]))))%*%testfield),type="Train"),
                    data.frame(lon=test_loc[,1],lat=test_loc[,2],rain=as.vector(t(t(Atest)*sqrt(abs(mesh_sim$loc[,1]*mesh_sim$loc[,2])))%*%testfield),type="Test"))

p.A.map.ns<-ggplot(Amap_ns.df)+gg(mesh_sim,edge.color = "gray",
                            int.color = "black",
                            ext.color = "black")+geom_point(aes(x=lon,y=lat,shape=type,color=rain),size=2)+xlab("longitude")+
  ylab("latitude")+scale_color_distiller(palette = "Spectral",limits=c(-max(abs(Amap_ns.df$rain)),max(abs(Amap_ns.df$rain))))
p.A.map.ns

m <- Matrix::t(m)
















############################ Different outlier


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

outlierval<- seq(5,15,5)
n_rep<-100#100
predictionscore.ll.10 <- rep(0,length(outlierval))
predictionscore.sroot.10 <- rep(0,length(outlierval))
predictionscore.ll.0 <- rep(0,length(outlierval))
predictionscore.sroot.0 <- rep(0,length(outlierval))
predictionscore.diff.0 <- rep(0,length(outlierval))
predictionscore.diff.10 <- rep(0,length(outlierval))
for(i in c(2:length(outlierval))){
  res_10 <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 10, outlier_val = outlierval[i],sigma_val=0.5,A=t(t(A)*sqrt(mesh_sim$loc[,1]^2+mesh_sim$loc[,2]^2)),Atest=t(t(Atest)*mesh_sim$loc[,1]),scoretypes=c("sroot","ll","slog","crps","scrps") ) #no outliers 0.0002

  res_0 <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 0, outlier_val = outlierval[i],sigma_val=0.5,A=t(t(A)*sqrt(mesh_sim$loc[,1]^2+mesh_sim$loc[,2]^2)),Atest=t(t(Atest)*mesh_sim$loc[,1]),scoretypes=c("sroot","ll","slog","crps","scrps") ) #no outliers 0.0002
  predictionscore.ll.10[i] <- mean(res_10$pred_ll)
  predictionscore.sroot.10[i] <- mean(res_10$pred_sroot)
  predictionscore.diff.10[i] <- mean((res_10$pred_ll-res_10$pred_sroot)/res_10$pred_sroot)
  predictionscore.ll.0[i] <- mean(res_0$pred_ll)
  predictionscore.sroot.0[i] <- mean(res_0$pred_sroot)
  predictionscore.diff.0[i] <- mean((res_0$pred_ll-res_0$pred_sroot)/res_0$pred_sroot)

}


plot(c(0,outlierval),exp(c(predictionscore.diff.0[1],predictionscore.diff.10)),xlab="outlier",
     ylab="exp(score diff)")
abline(a=1,b=(exp(predictionscore.diff.10[6])-exp(predictionscore.diff.0[1]))/150)


plot(c(0,outlierval),c(predictionscore.diff.0[1],predictionscore.diff.10),xlab="outlier",
     ylab="score diff")

plot(c(0,outlierval),-c(predictionscore.ll.0[1],predictionscore.ll.10))
lines(c(0,outlierval),-c(predictionscore.sroot.0[1],predictionscore.sroot.10))

plot(c(0,outlierval),-c(predictionscore.sroot.0[1],predictionscore.sroot.10))
lines(c(0,outlierval),-c(predictionscore.ll.0[1],predictionscore.ll.10))


res_10 <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 10, outlier_val = 150 ,sigma_val=0.5,A=t(t(A)*sqrt(mesh_sim$loc[,1]^2+mesh_sim$loc[,2]^2)),Atest=t(t(Atest)*mesh_sim$loc[,1]),scoretypes=c("sroot","ll") ) #no outliers 0.0002

predictionscore.ll.10 <- c(predictionscore.ll.10,mean(na.omit(res_10$pred_ll)))
predictionscore.sroot.10 <- c(predictionscore.sroot.10,mean(res_10$pred_sroot))
predictionscore.diff.10 <- c(predictionscore.diff.10,mean(na.omit((res_10$pred_ll-res_10$pred_sroot)/res_10$pred_sroot)))
outlierval <- c(outlierval,150)


hist(sapply(res_10$o_ll,function(x) x$par)[3,],breaks = 100)
hist(sapply(res_10$o_sroot,function(x) x$par)[3,],breaks=100)

hist((res_10$pred_ll-res_10$pred_sroot)/res_10$pred_sroot)


p.res_10 <- plot_results(res_10)

score_par <- sapply(res_10$o_sroot, function(o) o$par)
log_par <- sapply(res_10$o_ll, function(o) o$par)

#Check if all optimisations converged
print(sum(sapply(res_10$o_sroot, function(o) o$convergence)))
print(sum(sapply(res_10$o_ll, function(o) o$convergence)))

  #Join results in dataframe
  par_df <- data.frame(method=rep(c("Sroot","LL"),each=ncol(score_par)), par = Matrix::t(cbind(score_par,log_par)))

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


  p.hist.p3<-ggplot(par_df,aes(x=par.3,fill=method,color=method))+geom_histogram(alpha=0.5, position="identity")+geom_vline(xintercept = params_true[3])+
    geom_vline(data = par_df_mean,aes(xintercept=par.3_mean,color=method,linetype=method))+
    scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")
  p.box.p3 <- ggplot(par_df,aes(x=method, y=par.3,fill=method,color=method))+geom_boxplot(alpha=0.5, position="identity")+geom_hline(yintercept = params_true[3])+
    scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")
  p.scatter.2 <- ggplot(par_df,aes(x=par.2,y=par.3,color=method,shape=method))+geom_point(alpha=0.75)+scale_color_brewer(palette="Dark2")+annotate("point",x=params_true[2],y=params_true[3],col="black",shape=8)
  p.scatter.3 <- ggplot(par_df,aes(x=par.1,y=par.3,color=method,shape=method))+geom_point(alpha=0.75)+scale_color_brewer(palette="Dark2")+annotate("point",x=params_true[1],y=params_true[3],col="black",shape=8)

p.box.p1
p.box.p2
p.box.p3



########### rcrps
res_10_rcrps <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 10, outlier_val = 10 ,sigma_val=0.5,A=t(t(A)*sqrt(mesh_sim$loc[,1]^2+mesh_sim$loc[,2]^2)),Atest=t(t(Atest)*mesh_sim$loc[,1]),scoretypes=c("sroot","ll","rcrps") ) #no outliers 0.0002
