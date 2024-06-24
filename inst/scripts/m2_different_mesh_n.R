library(INLA)
library(grid)
library(cowplot)
library(latex2exp)

set.seed(20240621)

reslist <- vector("list", 4)
nresmesh <- rep(0,4)
reslist1d <- vector("list", 4)
nresmesh1d <- rep(0,4)
nlist <- c(100,200,500,1000)#,1500,2000)

for (i in c(1:length(nlist))){
  n<- nlist[i]
sim_loc = matrix(runif(2*n), ncol = 2, byrow = T)
mesh_sim = inla.mesh.2d(loc = sim_loc, max.edge=c(0.5, 1))
mesh_sim1d = inla.mesh.1d(loc = sim_loc[,1], max.edge=c(0.5, 1))
plot(mesh_sim)
points(sim_loc,col='red',pch=19)
nresmesh[i]<-mesh_sim$n
nresmesh1d[i]<-mesh_sim1d$n

spde = inla.spde2.matern(mesh_sim, alpha = 2)
spde1d = inla.spde2.matern(mesh_sim1d, alpha = 2)
params_true=spde$param.inla$theta.initial
#params_true<-c(log(0.1),params_true)
print(params_true)
Q = inla.spde.precision(spde, theta=spde$param.inla$theta.initial)
Q1d = inla.spde.precision(spde1d, theta=spde$param.inla$theta.initial)


# n_mesh= mesh_sim$n
n_rep=100
res_no_outliers <- repeated_inference(spde,mesh_sim$n,n_rep,Q) #no outliers
reslist[[i]]<-res_no_outliers
res_no_outliers1d <- repeated_inference(spde1d,mesh_sim1d$n,n_rep,Q1d) #no outliers
reslist1d[[i]]<-res_no_outliers1d
}

#save(reslist, file="m2_different_mesh_reslist.Rda")

rbind(data.frame(n = rep(nlist,each=n_rep), t=as.vector(sapply(reslist,function(x) x$times_sroot)),type="Sroot"),
      data.frame(n = rep(nlist,each=n_rep), t=as.vector(sapply(reslist,function(x) x$times_ll)),type="LL"))%>%ggplot(aes(x=as.factor(log(n)),y=log(t),color=type))+geom_boxplot()


rbind(data.frame(n = rep(nlist,each=n_rep), t=as.vector(sapply(reslist,function(x) x$times_sroot)),type="Sroot"),
      data.frame(n = rep(nlist,each=n_rep), t=as.vector(sapply(reslist,function(x) x$times_ll)),type="LL"))%>%ggplot(aes(x=log(n),y=log(t),color=type))+geom_boxplot()

rbind(data.frame(n = rep(nlist,each=n_rep), t=as.vector(sapply(reslist,function(x) x$times_sroot)),type="Sroot"),
      data.frame(n = rep(nlist,each=n_rep), t=as.vector(sapply(reslist,function(x) x$times_ll)),type="LL"))%>%ggplot(aes(x=log(n),y=log(t),color=type))+geom_boxplot(position="identity")


rbind(data.frame(n = rep(nlist,each=n_rep), t=as.vector(sapply(reslist,function(x) x$times_sroot)),type="Sroot"),
      data.frame(n = rep(nlist,each=n_rep), t=as.vector(sapply(reslist,function(x) x$times_ll)),type="LL"))%>%ggplot(aes(x=log(n),y=log(t),color=type))+geom_point()


rbind(data.frame(n = rep(nlist,each=n_rep), t=as.vector(sapply(reslist,function(x) x$times_sroot)),type="Sroot"),
      data.frame(n = rep(nlist,each=n_rep), t=as.vector(sapply(reslist,function(x) x$times_ll)),type="LL"))%>%ggplot(aes(x=log(n),y=log(t),color=type))+
  geom_point()+
  geom_line(data=predframe)+
  geom_ribbon(data=predframe,aes(ymin=lwr,ymax=upr),alpha=0.3)



library(scales)
rbind(data.frame(n = rep(nresmesh,each=n_rep), t=as.vector(t(sapply(reslist,function(x) x$times_sroot))),type="Sroot, 2d"),
      data.frame(n = rep(nresmesh,each=n_rep), t=as.vector(t(sapply(reslist,function(x) x$times_ll))),type="LL, 2d"),
      data.frame(n = rep(nresmesh1d,each=n_rep), t=as.vector(t(sapply(reslist1d,function(x) x$times_sroot))),type="Sroot, 1d"),
      data.frame(n = rep(nresmesh1d,each=n_rep), t=as.vector(t(sapply(reslist1d,function(x) x$times_ll))),type="LL, 1d")) %>%
  group_by(n,type) %>%
  summarize(min_t = min(t), max_t = max(t), mean_t=mean(t)) %>% ggplot(aes(x=n,y=mean_t,color=type,fill=type,shape=type))+geom_point()+geom_line()+geom_ribbon(aes(ymin=min_t,ymax=max_t),alpha=0.3,linetype="dashed")+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +ylab("time (sek)")
ggsave('m2_different_mesh_time.pdf',dpi = 1200,width = 10,height = 8,units = 'cm')


  # + coord_trans(x="log", y="log")

subset(rbind(data.frame(n = rep(nlist,each=n_rep), t=as.vector(sapply(reslist,function(x) x$times_sroot)),type="Sroot"),
             data.frame(n = rep(nlist,each=n_rep), t=as.vector(sapply(reslist,function(x) x$times_ll)),type="LL")) %>%
         group_by(n,type) %>%
         summarize(min_t = min(t), max_t = max(t), mean_t=mean(t)),log(n)>5) %>% group_by(type)%>%
  do({
    mod = lm(log(mean_t) ~ log(n), data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })







########################## just call score #####################
library(INLA)
library(grid)
library(cowplot)
library(latex2exp)

set.seed(20240621)

reslist <- vector("list", 21)
nresmesh <- rep(0,21)
reslist1d <- vector("list", 21)
nresmesh1d <- rep(0,21)
reslist3d <- vector("list", 21)
nresmesh3d <- rep(0,21)
nlist <- c(100,200,500,1000,1500,2000)

nlist <- exp(seq(2,8,0.7))
times_df <- data.frame()

for (i in c(1:length(nlist))){
  n<- nlist[i]
  sim_loc = matrix(runif(2*n), ncol = 2, byrow = T)
  mesh_sim = inla.mesh.2d(loc = sim_loc, max.edge=c(0.5, 1))
  mesh_sim1d = inla.mesh.1d(loc = sim_loc[,1], max.edge=c(0.5, 1))
  plot(mesh_sim)
  points(sim_loc,col='red',pch=19)
#  nresmesh[i]<-n#mesh_sim$n
#  nresmesh1d[i]<-n#mesh_sim1d$n
#  nresmesh3d[i]<-n#mesh_sim1d$n

  spde = inla.spde2.matern(mesh_sim, alpha = 2)
  spde1d = inla.spde2.matern(mesh_sim1d, alpha = 2)
  params_true=spde$param.inla$theta.initial
  #params_true<-c(log(0.1),params_true)
  print(params_true)
  Q = inla.spde.precision(spde, theta=spde$param.inla$theta.initial)
  Q1d = inla.spde.precision(spde1d, theta=spde$param.inla$theta.initial)
  # Q3d <- kronecker(Q,Q1d)


  # n_mesh= mesh_sim$n
  # n_rep=100
  # res_no_outliers <- repeated_inference(spde,mesh_sim$n,n_rep,Q) #no outliers
  # reslist[[i]]<-res_no_outliers
  # res_no_outliers1d <- repeated_inference(spde1d,mesh_sim1d$n,n_rep,Q1d) #no outliers
  # reslist1d[[i]]<-res_no_outliers1d
  if(mesh_sim$n<exp(9.3)){
    #Q = inla.spde.precision(spde, theta=spde$param.inla$theta.initial)
  mu<- rep(0,mesh_sim$n)
  #if(is.null(m)){
  m<-t(inla.qsample(n=1, Q = Q, mu=mu)) #10
  mbtest<-microbenchmark::microbenchmark(
    loo_score_vectorised(m,mu,Q,score="sroot"),
    log_dmvn(m,mu,Q),
    times=10)
  times_df <- rbind(times_df,data.frame(summary(mbtest,unit="ms"),n=mesh_sim$n,dim="2d"))
  }
  if(mesh_sim1d$n<exp(9.3)){
    #Q1d = inla.spde.precision(spde1d, theta=spde$param.inla$theta.initial)
  mu1d<- rep(0,mesh_sim1d$n)
  #if(is.null(m)){
  m1d<-t(inla.qsample(n=1, Q = Q1d, mu=mu1d)) #10
  mbtest1d<-microbenchmark::microbenchmark(
    loo_score_vectorised(m1d,mu1d,Q1d,score="sroot"),
    log_dmvn(m1d,mu1d,Q1d),
    times=10)
  times_df <- rbind(times_df,data.frame(summary(mbtest1d,unit="ms"),n=mesh_sim1d$n,dim="1d"))
  }
  if(mesh_sim$n*mesh_sim1d$n<exp(9.3)){
    Q3d <- kronecker(Q,Q1d)
  mu3d<- rep(0,mesh_sim$n*mesh_sim1d$n)
  #if(is.null(m)){
  m3d<-t(inla.qsample(n=1, Q = Q3d, mu=mu3d)) #10
  mbtest3d<-microbenchmark::microbenchmark(
    loo_score_vectorised(m3d,mu3d,Q3d,score="sroot"),
    log_dmvn(m3d,mu3d,Q3d),
    times=10)
  times_df <- rbind(times_df,data.frame(summary(mbtest3d,unit="ms"),n=mesh_sim$n*mesh_sim1d$n,dim="3d"))
  }


  # mbtest<-microbenchmark::microbenchmark(
  #   my_obj_func_3(c(kappa0,tau0)),
  #   my_log_obj_func(c(kappa0,tau0)),
  #   times=1
  # )

  #times_df <- rbind(times_df,data.frame(summary(mbtest,unit="ms"),n=mesh_sim$n,dim="2d"),data.frame(summary(mbtest1d,unit="ms"),n=mesh_sim1d$n,dim="1d"),data.frame(summary(mbtest3d,unit="ms"),n=mesh_sim$n*mesh_sim1d$n,dim="3d"))
}

levels(times_df$expr)<-c("2d, Sroot","2d, LL","1d, Sroot","1d, LL","3d, Sroot","3d, LL")
times_df$expr<-factor(times_df$expr,c("1d, Sroot","1d, LL","2d, Sroot","2d, LL","3d, Sroot","3d, LL"))

save(times_df, file="m2_different_mesh_reslist_diffdim.Rda")
ggplot(data=times_df,aes(x=log(n),y=log(mean),color=expr,shape=expr))+geom_point()+
  geom_smooth(data=subset(times_df,log(n)>4),method = "lm", se = FALSE, aes(fill=expr,linetype=expr),alpha=0.3)+scale_fill_manual(values=c("#95D840FF","#2D708EFF","#95D840FF","#2D708EFF","#95D840FF","#2D708EFF","#3CBB75FF","#39568CFF"))+scale_color_manual(values=c("#95D840FF","#2D708EFF","#95D840FF","#2D708EFF","#95D840FF","#2D708EFF","#3CBB75FF","#39568CFF"))+scale_shape_manual(values=c(19,19,17,17,15,15))+scale_linetype_manual(values=c("solid","dotted","solid","dotted","solid","dotted"))+
  ylab("log(mean runtime)")


ggplot(data=times_df,aes(x=n,y=mean,color=expr,shape=expr))+geom_point()+
  geom_smooth(data=subset(times_df,log(n)>4),method = "lm", se = FALSE, aes(fill=expr,linetype=expr),alpha=0.3,size=0.5)+scale_fill_manual(values=c("#95D840FF","#2D708EFF","#95D840FF","#2D708EFF","#95D840FF","#2D708EFF","#3CBB75FF","#39568CFF"))+scale_color_manual(values=c("#95D840FF","#2D708EFF","#95D840FF","#2D708EFF","#95D840FF","#2D708EFF","#3CBB75FF","#39568CFF"))+scale_shape_manual(values=c(19,19,17,17,15,15))+scale_linetype_manual(values=c("solid","solid","dotted","dotted","dashed","dashed"))+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +ylab("time (ms)")

ggplot(data=times_df,aes(x=n,y=mean,color=expr,shape=expr))+geom_point()+
  geom_smooth(data=subset(times_df,log(n)>4),method = "lm", se = FALSE, aes(fill=expr,linetype=expr),alpha=0.3,size=0.5)+scale_fill_manual(values=c("#1B9E77", "#D95F02","#1B9E77", "#D95F02","#1B9E77", "#D95F02"))+scale_color_manual(values=c("#1B9E77", "#D95F02","#1B9E77", "#D95F02","#1B9E77", "#D95F02"))+scale_shape_manual(values=c(17,17,19,19,15,15))+scale_linetype_manual(values=c("dotted","dotted","solid","solid","dashed","dashed"))+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +ylab("time (ms)")
ggsave('m2_different_mesh_time_diffdim.pdf',dpi = 1200,width = 12,height = 8,units = 'cm')

subset(times_df,log(n)>5) %>%
  group_by(expr) %>%
  do({
    mod = lm(log(mean) ~ log(n), data = .)
    # print(confint(mod))
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
