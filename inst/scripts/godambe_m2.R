#godambe numerical derivatives

# set.seed(20240629)
# n<-100
# sim_loc = matrix(runif(2*n), ncol = 2, byrow = T)
# mesh_sim = inla.mesh.2d(loc = sim_loc, max.edge=c(0.5, 1))
# plot(mesh_sim)
# points(sim_loc,col='red',pch=19)
# spde = inla.spde2.matern(mesh_sim, alpha = 2)


params_true=spde$param.inla$theta.initial
godambe_res_df<-data.frame()
for (p1 in seq(-4,-3,0.2)){
  for(p2 in seq(1,4,0.2)){
    EgradLOOS2crps <- Matrix(0,nrow=2,ncol=2)
    EgradLOOScrps <- Matrix(0,nrow=2,ncol=2)
    EhessLOOScrps <- Matrix(0,nrow=2,ncol=2)
    EgradLOOS2scrps <- Matrix(0,nrow=2,ncol=2)
    EgradLOOSscrps <- Matrix(0,nrow=2,ncol=2)
    EhessLOOSscrps <- Matrix(0,nrow=2,ncol=2)
    EgradLOOS2rcrps <- Matrix(0,nrow=2,ncol=2)
    EgradLOOSrcrps <- Matrix(0,nrow=2,ncol=2)
    EhessLOOSrcrps <- Matrix(0,nrow=2,ncol=2)
    EgradLOOS2 <- Matrix(0,nrow=2,ncol=2)
    EgradLOOS <- Matrix(0,nrow=2,ncol=2)
    EhessLOOS <- Matrix(0,nrow=2,ncol=2)
    EgradLL2 <- Matrix(0,nrow=2,ncol=2)
    EgradLL <- Matrix(0,nrow=2,ncol=2)
    EhessLL <- Matrix(0,nrow=2,ncol=2)

    params_true = c(p1,p2)
Q = inla.spde.precision(spde, theta=params_true)
n_dim <- nrow(Q)
mu<- rep(0,n_dim)

n_rep <- 500

for(i in c(1:n_rep)){ #repeat 100 times
  print(i)
m<-t(inla.qsample(n=1, Q = Q, mu=mu))

LOOscore <- function(par,scoretype="sroot"){
  theta <- par
  mu <- rep(0,n_dim)
  Qtheta <- inla.spde.precision(spde, theta=theta)
  return(loo_score_vectorised(m,mu,Qtheta,score=scoretype))
}

LL <- function(par){
  theta <- par
  mu <- rep(0,n_dim)
  Qtheta <- inla.spde.precision(spde, theta=theta)
  return(-log_dmvn(m,mu,Qtheta))
}

LOOSlog <- function(par,A,m){
  # print(par)
  theta <- par
  mu <- rep(0,n_dim)
  Qtheta <- inla.spde.precision(spde, theta=theta)
  score <- loo_log_score(m,muy,Qtheta)
  return(score)
}

gradLOOS<-pracma::grad(LOOscore,params_true)
hessLOOS<-pracma::hessian(LOOscore,params_true)
gradLOOScrps<-pracma::grad(LOOscore,params_true,scoretype="crps")
hessLOOScrps<-pracma::hessian(LOOscore,params_true,scoretype="crps")
gradLOOSscrps<-pracma::grad(LOOscore,params_true,scoretype="scrps")
hessLOOSscrps<-pracma::hessian(LOOscore,params_true,scoretype="scrps")
gradLOOSrcrps<-pracma::grad(LOOscore,params_true,scoretype="rcrps")
hessLOOSrcrps<-pracma::hessian(LOOscore,params_true,scoretype="rcrps")

gradLOOSlog<-pracma::grad(LOOSlog,params_true)
hessLOOSlog<-pracma::hessian(LOOSlog,params_true)

gradLL<-pracma::grad(LL,params_true)
hessLL<-pracma::hessian(LL,params_true)



EgradLOOS <-EgradLOOS + gradLOOS
EgradLOOS2 <-EgradLOOS2 + gradLOOS%*%t(gradLOOS)
EhessLOOS <-EhessLOOS + hessLOOS

EgradLOOScrps <-EgradLOOScrps + gradLOOScrps
EgradLOOS2crps <-EgradLOOS2crps + gradLOOScrps%*%t(gradLOOScrps)
EhessLOOScrps <-EhessLOOScrps + hessLOOScrps

EgradLOOSscrps <-EgradLOOSscrps + gradLOOSscrps
EgradLOOS2scrps <-EgradLOOS2scrps + gradLOOSscrps%*%t(gradLOOSscrps)
EhessLOOSscrps <-EhessLOOSscrps + hessLOOSscrps

EgradLOOSrcrps <-EgradLOOSrcrps + gradLOOSrcrps
EgradLOOS2rcrps <-EgradLOOS2rcrps + gradLOOSrcrps%*%t(gradLOOSrcrps)
EhessLOOSrcrps <-EhessLOOSrcrps + hessLOOSrcrps

EgradLOOSlog <-EgradLOOSlog + gradLOOSlog
EgradLOOS2log <-EgradLOOS2log + gradLOOSlog%*%t(gradLOOSlog)
EhessLOOSlog <-EhessLOOSlog + hessLOOSlog

EgradLL <-EgradLL + gradLL
EgradLL2 <-EgradLL2 + gradLL%*%t(gradLL)
EhessLL <-EhessLL + hessLL

}

EgradLOOScrps <- EgradLOOScrps/n_rep
EgradLOOSscrps <- EgradLOOSscrps/n_rep
EgradLOOSrcrps <- EgradLOOSrcrps/n_rep
EgradLOOS <- EgradLOOS/n_rep
EgradLOOSlog <- EgradLOOSlog/n_rep
EgradLL <- EgradLL/n_rep

EgradLOOS2crps <- EgradLOOS2crps/n_rep
EgradLOOS2scrps <- EgradLOOS2scrps/n_rep
EgradLOOS2rcrps <- EgradLOOS2rcrps/n_rep
EgradLOOS2 <- EgradLOOS2/n_rep
EgradLOOS2log <- EgradLOOS2log/n_rep
EgradLL2 <- EgradLL2/n_rep


EhessLOOScrps <- EhessLOOScrps/n_rep
EhessLOOSscrps <- EhessLOOSscrps/n_rep
EhessLOOSrcrps <- EhessLOOSrcrps/n_rep
EhessLOOS <- EhessLOOS/n_rep
EhessLOOSlog <- EhessLOOSlog/n_rep
EhessLL <- EhessLL/n_rep

VLOOScrps <- solve(EhessLOOScrps,EgradLOOS2crps)%*%t(solve(EhessLOOScrps))
VLOOSscrps <- solve(EhessLOOSscrps,EgradLOOS2scrps)%*%t(solve(EhessLOOSscrps))
VLOOSrcrps <- solve(EhessLOOSrcrps,EgradLOOS2rcrps)%*%t(solve(EhessLOOSrcrps))
VLOOS <- solve(EhessLOOS,EgradLOOS2)%*%t(solve(EhessLOOS))
VLOOSlos <- solve(EhessLOOSlog,EgradLOOS2log)%*%t(solve(EhessLOOSlog))
VLL<-solve(EhessLL,EgradLL2)%*%t(solve(EhessLL))


diag(VLOOScrps)
diag(VLOOSscrps)
diag(VLOOSrcrps)
diag(VLOOS)
diag(VLOOSlog)
diag(VLL)


tmptitle <- paste("godambe_m2",format(Sys.time(), "%Y%m%d_%H%M%S"),sep = "")
save(VLOOScrps,VLOOSscrps,VLOOSrcrps,VLOOS,VLOOSlog,VLL,params_true,file=paste(tmptitle,".Rda",sep=""))
godambe_res_df<-rbind(godambe_res_df,data.frame(god.sd1 = sqrt(c(diag(VLOOScrps)[1],diag(VLOOSscrps)[1],diag(VLOOSrcrps)[1],diag(VLOOS)[1],diag(VLL)[1])),p1=p1,p2=p2,
god.sd2 = sqrt(c(diag(VLOOScrps)[2],diag(VLOOSscrps)[2],diag(VLOOSrcrps)[2],diag(VLOOS)[2],diag(VLL)[2])),type = c("crps","scrps","rcrps","sroot","slog","ll")))
  }
}

m2_test <- repeated_inference(spde,mesh_sim$n,100,Q,n_outlier=0,outlier_val=NULL,m=NULL)
par_sroot <- sapply(m2_test$o_sroot,function(o) o$par)
par_ll <- sapply(m2_test$o_ll,function(o) o$par)
par_crps <- sapply(m2_test$o_crps,function(o) o$par)
par_scrps <- sapply(m2_test$o_scrps,function(o) o$par)
par_rcrps <- sapply(m2_test$o_rcrps,function(o) o$par)
save(m2_test,file=paste("m2test",tmptitle,".Rda",sep=""))
godambe_df <- data.frame(sd.par1=c(sd(par_crps[1,]),sd(par_scrps[1,]),sd(par_rcrps[1,]),sd(par_sroot[1,]),sd(par_ll[1,])),
                         sd.par2=c(sd(par_crps[2,]),sd(par_scrps[2,]),sd(par_rcrps[2,]),sd(par_sroot[2,]),sd(par_ll[2,])),
                         god.sd1 = sqrt(c(diag(VLOOScrps)[1],diag(VLOOSscrps)[1],diag(VLOOSrcrps)[1],diag(VLOOS)[1],diag(VLL)[1])),
                         god.sd2 = sqrt(c(diag(VLOOScrps)[2],diag(VLOOSscrps)[2],diag(VLOOSrcrps)[2],diag(VLOOS)[2],diag(VLL)[2])),
                         type = c("crps","scrps","rcrps","sroot","ll"))
godambe_df$n=100

m2_test200 <- repeated_inference(spde,mesh_sim$n,200,Q,n_outlier=0,outlier_val=NULL,m=NULL)
par_sroot <- sapply(m2_test200$o_sroot,function(o) o$par)
par_ll <- sapply(m2_test200$o_ll,function(o) o$par)
par_crps <- sapply(m2_test200$o_crps,function(o) o$par)
par_scrps <- sapply(m2_test200$o_scrps,function(o) o$par)
par_rcrps <- sapply(m2_test200$o_rcrps,function(o) o$par)
save(m2_test200,file=paste("m2test200",tmptitle,".Rda",sep=""))
godambe_df <- rbind(godambe_df,data.frame(sd.par1=c(sd(par_crps[1,]),sd(par_scrps[1,]),sd(par_rcrps[1,]),sd(par_sroot[1,]),sd(par_ll[1,])),
                         sd.par2=c(sd(par_crps[2,]),sd(par_scrps[2,]),sd(par_rcrps[2,]),sd(par_sroot[2,]),sd(par_ll[2,])),
                         god.sd1 = sqrt(c(diag(VLOOScrps)[1],diag(VLOOSscrps)[1],diag(VLOOSrcrps)[1],diag(VLOOS)[1],diag(VLL)[1])),
                         god.sd2 = sqrt(c(diag(VLOOScrps)[2],diag(VLOOSscrps)[2],diag(VLOOSrcrps)[2],diag(VLOOS)[2],diag(VLL)[2])),
                         type = c("crps","scrps","rcrps","sroot","ll"),n=200))

m2_test10 <- repeated_inference(spde,mesh_sim$n,10,Q,n_outlier=0,outlier_val=NULL,m=NULL)
par_sroot <- sapply(m2_test10$o_sroot,function(o) o$par)
par_ll <- sapply(m2_test10$o_ll,function(o) o$par)
par_crps <- sapply(m2_test10$o_crps,function(o) o$par)
par_scrps <- sapply(m2_test10$o_scrps,function(o) o$par)
par_rcrps <- sapply(m2_test10$o_rcrps,function(o) o$par)
save(m2_test200,file=paste("m2test10",tmptitle,".Rda",sep=""))
godambe_df <- rbind(godambe_df,data.frame(sd.par1=c(sd(par_crps[1,]),sd(par_scrps[1,]),sd(par_rcrps[1,]),sd(par_sroot[1,]),sd(par_ll[1,])),
                                          sd.par2=c(sd(par_crps[2,]),sd(par_scrps[2,]),sd(par_rcrps[2,]),sd(par_sroot[2,]),sd(par_ll[2,])),
                                          god.sd1 = sqrt(c(diag(VLOOScrps)[1],diag(VLOOSscrps)[1],diag(VLOOSrcrps)[1],diag(VLOOS)[1],diag(VLL)[1])),
                                          god.sd2 = sqrt(c(diag(VLOOScrps)[2],diag(VLOOSscrps)[2],diag(VLOOSrcrps)[2],diag(VLOOS)[2],diag(VLL)[2])),
                                          type = c("crps","scrps","rcrps","sroot","ll"),n=10))


godambe_df %>% tidyr::gather(partype,sd,-type,-n)%>%ggplot(aes(x=partype,y=sd,color=type,shape=as.factor(n)))+geom_point()



godambe_res_df%>%ggplot(aes(x=p1,y=god.sd1,color=type,shape=as.factor(p2)))+geom_point()+geom_line()
subset(godambe_res_df,p1==-3.6)%>%ggplot(aes(x=p2,y=god.sd2,color=type,shape=as.factor(p1)))+geom_point()+geom_line()
godambe_res_df$type<-as.factor(godambe_res_df$type)
levels(godambe_res_df$type)<-c("CRPS","LL","rCRPS","SCRPS","Sr")
godambe_res_df$type<-factor(godambe_res_df$type,c("LL","SCRPS","Sr","CRPS","rCRPS"))
p.g.22<-subset(godambe_res_df,p1==-3.6&p2>2)%>%ggplot(aes(x=p2,y=god.sd2,color=type,linetype=type))+geom_line()+scale_color_viridis_d()+xlab(TeX("$\\log(\\tau)$"))+ylab(TeX("sd $\\log(\\tau)$"))
p.g.12<-subset(godambe_res_df,p2==2.6)%>%ggplot(aes(x=p1,y=god.sd2,color=type,linetype=type))+geom_line()+scale_color_viridis_d()+xlab(TeX("$\\log(\\kappa)$"))+ylab(TeX("sd $\\log(\\tau)$"))
p.g.21<-subset(godambe_res_df,p1==-3.6&p2>2)%>%ggplot(aes(x=p2,y=god.sd1,color=type,linetype=type))+geom_line()+scale_color_viridis_d()+xlab(TeX("$\\log(\\tau)$"))+ylab(TeX("sd $\\log(\\kappa)$"))
p.g.11<-subset(godambe_res_df,p2==2.6)%>%ggplot(aes(x=p1,y=god.sd1,color=type,linetype=type))+geom_line()+scale_color_viridis_d()+xlab(TeX("$\\log(\\kappa)$"))+ylab(TeX("sd $\\log(\\kappa)$"))

p.g.22<-subset(godambe_res_df,p1==-3.6&p2>2)%>%ggplot(aes(x=p2,y=god.sd2,color=type,linetype=type))+geom_line()+scale_color_brewer(palette="Dark2")+xlab(TeX("$\\log(\\tau)$"))+ylab(TeX("sd $\\log(\\tau)$"))
p.g.12<-subset(godambe_res_df,p2==2.6)%>%ggplot(aes(x=p1,y=god.sd2,color=type,linetype=type))+geom_line()+scale_color_brewer(palette="Dark2")+xlab(TeX("$\\log(\\kappa)$"))+ylab(TeX("sd $\\log(\\tau)$"))
p.g.21<-subset(godambe_res_df,p1==-3.6&p2>2)%>%ggplot(aes(x=p2,y=god.sd1,color=type,linetype=type))+geom_line()+scale_color_brewer(palette="Dark2")+xlab(TeX("$\\log(\\tau)$"))+ylab(TeX("sd $\\log(\\kappa)$"))
p.g.11<-subset(godambe_res_df,p2==2.6)%>%ggplot(aes(x=p1,y=god.sd1,color=type,linetype=type))+geom_line()+scale_color_brewer(palette="Dark2")+xlab(TeX("$\\log(\\kappa)$"))+ylab(TeX("sd $\\log(\\kappa)$"))


plot_grid_4(p.g.11,p.g.12,p.g.21,p.g.22)
tmptitle <- paste("godambe_res_m2",format(Sys.time(), "%Y%m%d_%H%M%S"),sep = "")
ggsave(paste(tmptitle,"_lines.pdf",sep=""),dpi = 1200,width = 12,height = 12,units = 'cm')


ggsave("godambe_res_m220240630_103127_lines_dark2.pdf",dpi = 1200,width = 12,height = 12,units = 'cm')
subset(godambe_res_df,type=="LL")%>%ggplot(aes(x=p1,y=god.sd1,color=as.factor(p2)))+geom_point()+geom_line()+scale_color_viridis_d()+xlab(TeX("$\\log(\\tau)$"))+ylab(TeX("sd $\\log(\\tau)$"))
save(godambe_res_df,file=paste(tmptitle,".Rda",sep=""))
