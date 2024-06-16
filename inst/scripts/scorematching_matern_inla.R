library(INLA)
library(grid)
library(cowplot)
library(latex2exp)

set.seed(20240509)
sim_loc = matrix(c(0,0,5,5, 0, 5, 5, 0), nrow = 4, byrow = T)
mesh_sim = inla.mesh.2d(loc = sim_loc, max.edge=c(0.5, 1))
plot(mesh_sim)
points(sim_loc,col='red',pch=19)
mesh_sim$n

spde = inla.spde2.matern(mesh_sim, alpha = 2)
params_true=spde$param.inla$theta.initial
#params_true<-c(log(0.1),params_true)
print(params_true)
Q = inla.spde.precision(spde, theta=spde$param.inla$theta.initial)


# n_mesh= mesh_sim$n
n_rep=100
res_no_outliers <- repeated_inference(spde,mesh_sim$n,n_rep,Q) #no outliers
res_50_outliers <- repeated_inference(spde,mesh_sim$n,n_rep,Q,5,4) #5 outliers (50%)
res_100_outliers <- repeated_inference(spde,mesh_sim$n,n_rep,Q,10,4) #10 outliers (100%)
res_25_outliers <- repeated_inference(spde,mesh_sim$n,n_rep,Q,25,4) #5 outliers (25%)
res_75_outliers <- repeated_inference(spde,mesh_sim$n,n_rep,Q,75,4) #5 outliers (75%)

p.res_no_outliers <- plot_results(res_no_outliers,extended = FALSE)
p.res_5_outliers <- plot_results(res_50_outliers,extended = FALSE)
p.res_10_outliers <- plot_results(res_100_outliers,extended = FALSE)


########################### Save results #######################

ggsave(paste("GMRF_scatter_est_no_5_10","_",n_rep,"rep",".pdf",sep=""),plot_grid_3(p.res_no_outliers$p.scatter,p.res_5_outliers$p.scatter,p.res_10_outliers$p.scatter), dpi=1200,width=18,height=8,unit="cm")
ggsave("GMRF_p1_est_no_5_10.pdf",plot_grid_3(p.res_no_outliers$p.hist.p1,p.res_5_outliers$p.hist.p1,p.res_10_outliers$p.hist.p1), dpi=1200,width=18,height=8,unit="cm")
ggsave("GMRF_p2_est_no_5_10.pdf",plot_grid_3(p.res_no_outliers$p.hist.p2,p.res_5_outliers$p.hist.p2,p.res_10_outliers$p.hist.p2), dpi=1200,width=18,height=8,unit="cm")
ggsave("GMRF_time_est_no_5_10.pdf",plot_grid_3(p.res_no_outliers$p.time,p.res_5_outliers$p.time,p.res_10_outliers$p.time), dpi=1200,width=18,height=8,unit="cm")
ggsave("GMRF_time_hist_est_no_5_10.pdf",plot_grid_3(p.res_no_outliers$p.time.hist+xlim(0,2)+ylim(0,61),p.res_5_outliers$p.time.hist+xlim(0,2)+ylim(0,61),p.res_10_outliers$p.time.hist+xlim(0,2)+ylim(0,61)), dpi=1200,width=18,height=8,unit="cm")


################# normal response model
n_rep<-100
res_no_outliers_nresp <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,sigma_val=0.1) #no outliers 0.0002
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
ggsave('inlaparhist_sigma01.pdf',p.par.hist,dpi = 1200,width = 18,height = 8,units = 'cm')
ggsave('inlaparscatter_sigma01.pdf',p.res_no_outliers_nresp$p.scatter+xlab(TeX("$\\log(\\kappa)$"))+ylab(TeX("$\\log(\\tau)$")),dpi = 1200,width = 18,height = 8,units = 'cm')
ggsave('inlatimehist_sigma01.pdf',p.res_no_outliers_nresp$p.time.hist,dpi = 1200,width = 18,height = 8,units = 'cm')

score_par <- sapply(res_no_outliers_nresp$o_sroot[1:n_res], function(o) o$par)
print(sum(sapply(res_no_outliers_nresp$o_sroot[1:n_res], function(o) o$convergence)))
par_df <- data.frame(method=rep(c("Sroot"),each=n_res), mu = c(score_par[1,]), rho = c(score_par[2,]), run.time=c(res_no_outliers_nresp$times_sroot[1:n_res]),i=rep(c(1:n_res),1))
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



################################# normal response model once
spde$param.inla$theta.initial
Q = inla.spde.precision(spde, theta=spde$param.inla$theta.initial)
sigma_val <- 0.001
n <- mesh_sim$n
mu<- rep(0,n)
I<-Diagonal(n)
A<-Diagonal(n)
m<-t(inla.qsample(n=100, Q = Q, mu=mu)) #observations of latent field
m <- m+rnorm(n=n,mean = 0,sd = sigma_val) #added noise
#m[,1]<-4
# if(n_outlier>0){ #add outliers
#   m[ matrix(c(1:n_outlier,sample(1:n,n_outlier)),ncol=2)]<-outlier_val #set random observation at each sample to 4.
# }




my_obj_func_new <- function(par){
  print(par)
  theta <- par
  mu <- rep(0,n)
  #Qxy <- inla.spde.precision(spde, theta=theta)+I*sigma_val^2
  Qx <- inla.spde.precision(spde, theta=theta)
  mux <- mu
  #if(sigma_val>0) muxy<- muxy+solve(Qxy,(I/sigma_val^2)%*%(t(m)-I%*%mu))
  return(loo_score_vectorised_eps(m,mux,Qx,sigma_val,A))
}

kappa0<-1
tau0<-0.1


o1<-optim(par=c(kappa0,tau0),my_obj_func_new,control=list(maxit=50000))



spde$param.inla$theta.initial
my_obj_func_new(spde$param.inla$theta.initial)

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

Q
A=Diagonal(nrow(Q))
Qe = 0.5^2*Diagonal(nrow(Q))

Ahat <- function(A,i){
  Ahat <- A
  Ahat[i,]<-0
  return(Ahat)
}

Qxy <- Q + t(A)%*%Qe%*%A
invQxy <- solve(Qxy)

invQe <- solve(Qe)
invQe <- solve(Qe[-1,-1]) #n-1 dimensions
tic()
muxyi <- sapply(c(1:nrow(Q)),function(i) mu + inv_sherman_morrison(invQxy,0.5*A[i,],-0.5*A[i,])%*%t(A[-i,])%*%invQe[-i,-i]%*%(A[-i,]%*%t(m)-A[-i,]%*%mu))
params<-sapply(c(1:nrow(Q)),function(i) c(as.numeric(A[i,]%*%muxyi[[i]]),as.numeric(A[i,]%*%inv_sherman_morrison(invQxy,0.5*A[i,],-0.5*A[i,])%*%A[i,]+sigma_val^2)))
toc()

tic()
  invQxyi<- sapply(c(1:nrow(Q)),function(i) inv_sherman_morrison(invQxy,0.5*A[i,],-0.5*A[i,]))
  muxyi <- sapply(c(1:nrow(Q)),function(i) mu + invQxyi[[i]]%*%t(A[-i,])%*%invQe[-i,-i]%*%(A[-i,]%*%t(m)-A[-i,]%*%mu))
  params<-sapply(c(1:nrow(Q)),function(i) c(as.numeric(A[i,]%*%muxyi[[i]]),as.numeric(A[i,]%*%invQxyi[[i]]%*%A[i,]+sigma_val^2)))
toc()


tic()
for (i in c(1:nrow(Q))){
  invQxyi<- inv_sherman_morrison(invQxy,0.5*A[i,],-0.5*A[i,])
}
toc()

(invQ-(invQ%*%(0.5*A[i,])%*%t(-0.5*A[i,])%*%invQ)/as.numeric(1+t(-0.5*A[i,])%*%invQ%*%(0.5*A[i,])))

tic()
bla<-apply(A,1,function(x) (invQxy-(invQxy%*%(0.5*x)%*%t(-0.5*x)%*%invQxy)/as.numeric(1+t(-0.5*x)%*%invQxy%*%(0.5*x))))
toc()

tic()
for (i in c(1:nrow(Q))){
  invQxyi<- inv_sherman_morrison_I(invQxy,i)
}
toc()

tic()
  invQxyi<- sapply(c(1:nrow(Q)),function(i)inv_sherman_morrison_I(invQxy,i,sigma_val))
toc()


tic()
for (i in c(1:nrow(Q))){
invQxyi<- inv_sherman_morrison(invQxy,0.5*A[i,],-0.5*A[i,])

#max(abs(invQxyi-invQxyhat))

# muxyi <- mu + invQxyi%*%t(Ahat(A,i))%*%invQe%*%(Ahat(A,i)%*%t(m)-Ahat(A,i)%*%mu)
muxyi <- mu + invQxyi%*%t(A[-i,])%*%invQe[-i,-i]%*%(A[-i,]%*%t(m)-A[-i,]%*%mu)

params<-c(as.numeric(A[i,]%*%muxyi),as.numeric(A[i,]%*%invQxyi%*%A[i,]+sigma_val^2))
}
toc()

invQxy
(invQxy%*%A[i,]%*%t(A[i,])%*%invQxy)/as.numeric(1+t(A[i,])%*%invQxy%*%A[i,])
# Ahat=A
# Ahat[2,]=0
# Ai = A[2,]
#
#
# t(A)%*%Qe%*%A-t(Ahat)%*%Qe%*%Ahat
# Ai%*%t(Ai)/2
#
#
# inla.qinv(Q) #compute diagonals of the inverse of the precision matrix




###################### parallelise

library(parallel)
library(MASS)
numCores <- detectCores()
numCores
system.time(
  results <- lapply(c(1:nrow(A)),function(i) (invQxy-(invQxy%*%(0.5*A[i,])%*%t(-0.5*A[i,])%*%invQxy)/as.numeric(1+t(-0.5*A[i,])%*%invQxy%*%(0.5*A[i,]))))
)
system.time(
  results <- mclapply(c(1:nrow(A)),function(i) (invQxy-(invQxy%*%(0.5*A[i,])%*%t(-0.5*A[i,])%*%invQxy)/as.numeric(1+t(-0.5*A[i,])%*%invQxy%*%(0.5*A[i,]))))
)



#################################33
#
#
#

Q
A=Diagonal(nrow(Q))
Qe = 0.5^2*Diagonal(nrow(Q))

Ahat <- function(A,i){
  Ahat <- A
  Ahat[i,]<-0
  return(Ahat)
}

Qxy <- Q + t(A)%*%Qe%*%A
tic()
invQxy <- solve(Qxy)
toc()

invQe <- solve(Qe)
invQe <- solve(Qe[-1,-1, drop = FALSE]) #n-1 dimensions
tic()
resA <- sapply(c(1:nrow(Q)),function(i)  inv_sherman_morrison(invQxy,0.5*A[i,, drop = FALSE],-0.5*A[i,, drop = FALSE]))
toc()

tic()
resB <- sapply(c(1:nrow(Q)),function(i)  solve(Qxy-0.5^2*(t(A[i,, drop = FALSE])%*%A[i,, drop = FALSE])))
toc()

tic()
muxyi <- sapply(c(1:nrow(Q)),function(i)  solve(Qxy))
toc()

tic()
muxyi <- sapply(c(1:nrow(Q)),function(i)  inv_sherman_morrison_I(invQxy,i,0.5))
toc()

tic()
resC<-sapply(c(1:nrow(Q)),function(i) invQxy%*%(Matrix::t(A[i,, drop = FALSE])%*%A[i,, drop = FALSE])%*%invQxy)
toc()

tic()
resC<-sapply(c(1:nrow(Q)),function(i) invQxy%*%(A[i,]%*%Matrix::t(A[i,]))%*%invQxy)
toc()

tic()
resD<-sapply(c(1:nrow(Q)),function(i) (invQxy%*%Matrix::t(A[i,, drop = FALSE]))%*%(A[i,, drop = FALSE]%*%invQxy))
toc()

tic()
resE<-sapply(c(1:nrow(Q)),function(i) invQxy%*%Matrix::t(A[i,, drop = FALSE])%*%A[i,, drop = FALSE]%*%invQxy)
toc()

tic()
muxyi <- sapply(c(1:nrow(Q)),function(i)  invQxy+invQxy[i,i])
toc()

tic()
muxyi <- sapply(c(1:nrow(Q)),function(i)  inv_sherman_morrison(invQxy,t(0.5*A[i,]),-0.5*A[i,]))
toc()


params<-sapply(c(1:nrow(Q)),function(i) c(as.numeric(A[i,]%*%muxyi[[i]]),as.numeric(A[i,]%*%inv_sherman_morrison(invQxy,0.5*A[i,],-0.5*A[i,])%*%A[i,]+sigma_val^2)))
toc()


i<-1
(invQxy-(invQxy%*%(0.5*A[i,])%*%t(-0.5*A[i,])%*%invQxy)/as.numeric(1+t(-0.5*A[i,])%*%invQxy%*%(0.5*A[i,])))
(invQxy-(invQxy%*%t(0.5*A[i,, drop = FALSE])%*%(-0.5*A[i,, drop = FALSE])%*%invQxy)/as.numeric(1+t(-0.5*A[i,,drop=FALSE])%*%invQxy%*%(0.5*A[i,])))
