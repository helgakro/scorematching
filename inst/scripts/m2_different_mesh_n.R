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
