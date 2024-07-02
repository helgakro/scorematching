library(devtools)
load_all()
library(INLA)

load("inst/data/taperexample/anom1962.RData")

#idxneusa <- (loc[,1]>(-120)&loc[,1]<(-100)&loc[,2]>35)#(loc[,1]>-80 & loc[,2]>40) #choose norhteast usa
idxneusa <- z #sample(c(1:length(z)),1000)#randomly chosen subset
loc<-loc[idxneusa,]
z<-z[idxneusa]
z<-z-mean(z)#normalise for neusa
z<-t(z)
mesh_sim = inla.mesh.2d(loc = loc, max.edge=c(5, 10))
plot(mesh_sim)
points(loc,col='red',pch=19)
mesh_sim$n
A <- inla.spde.make.A(mesh=mesh_sim,loc=loc)

spde = inla.spde2.matern(mesh_sim, alpha = 2)
sigma_val <- 0.001
Q = inla.spde.precision(spde, theta=spde$param.inla$theta.initial)

library(ggplot2)
library(inlabru)
neus.df <- data.frame(lon=loc[,1],lat=loc[,2],rain=t(z))
p.anomalies.map<-ggplot(neus.df)+gg(mesh_sim,edge.color = "gray",
                   int.color = "black",
                   ext.color = "black")+geom_point(aes(x=lon,y=lat,color=rain),size=0.9)+scale_color_viridis_c()+xlab("longitude")+ylab("latitude")
p.anomalies.map
#ggsave('anomalies_map_all.pdf',p.anomalies.map,dpi = 1200,width = 12,height = 10,units = 'cm')
################################## timing inverses
# n_y <- nrow(A)
# I<-Diagonal(n_y)
# Qx <- Q
# Qeps <- I/0.5^2
#
# microbenchmark::microbenchmark(
# Qx <- solve(A%*%solve(Qx)%*%Matrix::t(A)), times=1
# )
#
# microbenchmark::microbenchmark(
#   Qtheta1 <- Qx-Qx%*%solve(Qx+Qeps)%*%Qx,
#   #method 2
#   Qtheta2 <- Qeps-Qeps%*%solve(Qx+Qeps)%*%Qeps,
#   #method 3
#   Qtheta3 <- Qeps-Qeps%*%A%*%solve(Q+Matrix::t(A)%*%Qeps%*%A)%*%Matrix::t(A)%*%Qeps, times = 1)



########### rep optim ##############
# microbenchmark::microbenchmark(
# res_no_outliers_nresp <- inference_norm_resp(z,spde,mesh_sim$n,Q,sigma_val=sigma_val,A=A), #no outliers 0.0002,
# res_no_outliers_nresp_old <- inference_norm_resp_old(z,spde,mesh_sim$n,Q,sigma_val=sigma_val,A=A), #no outliers 0.0002
# times=1)

res_no_outliers_nresp$o1
res_no_outliers_nresp$o2
res_no_outliers_nresp$o3

res_no_outliers_nresp$t1
res_no_outliers_nresp$t2
res_no_outliers_nresp$t3

res_no_outliers_nresp_old$o1
res_no_outliers_nresp_old$o2
res_no_outliers_nresp_old$o3


# > res_no_outliers_nresp$o1
# $par
# [1] -1.2380598 -1.1209509  0.5442361
# $value
# [1] 0.6570577
# $counts
# function gradient
# 78       NA
# $convergence
# [1] 0
# $message
# NULL
#
# > res_no_outliers_nresp$o2
# $par
# [1] -0.6313840 -1.4512003  0.2900854
# $value
# [1] 750.7099
# $counts
# function gradient
# 160       NA
# $convergence
# [1] 0
# $message
# NULL
#
# > res_no_outliers_nresp$o3
# $par
# [1] -1.1697293 -1.2403658  0.6511273
# $value
# [1] 0.4921993
# $counts
# function gradient
# 140       NA
# $convergence
# [1] 0
# $message
# NULL

res_no_outliers_nresp <- inference_norm_resp(z,spde,mesh_sim$n,Q,A=A,ll=FALSE,slog=FALSE,sroot=FALSE,crps=FALSE,scrps=TRUE) #no outliers 0.0002,


idxoutlier <- z<3
locoutlier<-loc[idxoutlier,]
zoutlier<-z[idxoutlier]
zoutlier<-zoutlier-mean(zoutlier)#normalise for neusa
zoutlier <- t(zoutlier)
Aoutlier <- inla.spde.make.A(mesh=mesh_sim,loc=locoutlier)

res_outliers_nresp <- inference_norm_resp(zoutlier,spde,mesh_sim$n,Q,A=Aoutlier) #no outliers 0.0002,


res_no_outliers_nresp$o1
par.res <- exp(res_no_outliers_nresp$o1$par)
par.res
nufield <- 2-2/2
sfield = gamma(nufield)/(par.res[3]^2*par.res[2]^(2*nufield)*(4*pi)^(2/2)*gamma(2))
sfield
range = sqrt(8)/par.res[1]
range
res_outliers_nresp$o1


res_no_outliers_nresp$o2
res_outliers_nresp$o2


res_no_outliers_nresp$o3
res_outliers_nresp$o3

res_no_outliers_nresp$score_ll
res_outliers_nresp$score_ll
############# direct observation ######################
res_no_outliers_fresp <- inference_fix_resp(z,spde,mesh_sim$n,Q,A=A)
res_outliers_fresp <- inference_fix_resp(zoutlier,spde,mesh_sim$n,Q,A=Aoutlier) #no outliers 0.0002,

res_no_outliers_fresp$o1
res_outliers_fresp$o1




####### SCRPS

res_no_outliers_nresp$o5
par.res <- exp(res_no_outliers_nresp$o5$par)
par.res
nufield <- 2-2/2
sfield = gamma(nufield)/(par.res[3]^2*par.res[2]^(2*nufield)*(4*pi)^(2/2)*gamma(2))
sfield
range = sqrt(8)/par.res[2]
range
par.res[1]

