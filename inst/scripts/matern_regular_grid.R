library(INLA)

n_dim <- 1000#50000
h<-1
kappa <- 1
tau <- 2

getQ_regular <- function(kappa,tau,n_dim,h=1){
  C<-1/h*Matrix::Diagonal(n_dim,1)
  C[abs(row(C) - col(C)) == 1] <- -1/(3*h)
  G = 1/h*Matrix::Diagonal(n_dim, 2)
  G[abs(row(G) - col(G)) == 1] <- -1/h

  K <- kappa^2*C+G
  Q <- tau^2*K%*%solve(C)%*%K
  Q[abs(Q)<0.01]<-0
  Q<-Matrix::Matrix(Q,sparse = TRUE)

  image(Q)
  return(Q)
}

getQ_regular_multi <- function(kappa,tau,n_dim_1,n_dim_2,h=1){
  Q<-kronecker(getQ_regular(kappa,tau,n_dim_1,h),getQ_regular(kappa,tau,n_dim_2,h))
  return(Q)
}

#print(paste("Iteration",i_n))
n <- n_dim
mu<- rep(0,n)
#if(is.null(m)){
Q<-getQ_regular(kappa,tau,n_dim,h) #getQ_regular(kappa,tau,n_dim,h)
m<-t(inla.qsample(n=10, Q = Q, mu=mu))
#m[,1]<-4
#if(n_outlier>0){
#  m[ matrix(c(1:n_outlier,sample(1:n,n_outlier)),ncol=2)]<-outlier_val #set random observation at each sample to 4.
#}
#}


my_obj_func_3 <- function(par,scoretype="sroot"){
  theta <- par
  kappa <- exp(par[1])
  tau <- exp(par[2])
  mu <- rep(0,n)
  Qtheta <- getQ_regular(kappa,tau,n_dim,h)
  return(loo_score_vectorised(m,mu,Qtheta,score=scoretype))
}

my_log_score_obj_func <- function(par){
  theta <- par
  mu <- rep(0,n)
  kappa <- exp(par[1])
  tau <- exp(par[2])
  Qtheta <- getQ_regular(kappa,tau,n_dim,h)
  return(loo_log_score(m,mu,Qtheta))
}

my_log_obj_func <- function(par){
  theta <- par
  kappa <- exp(par[1])
  tau <- exp(par[2])
  mu <- rep(0,n)
  Qtheta <- getQ_regular(kappa,tau,n_dim,h)
  return(-log_dmvn(m,mu,Qtheta))
}

# kappa0<--2
# tau0<-1

kappa0<--1
tau0<-1

microbenchmark::microbenchmark(
  my_obj_func_3(c(kappa0,tau0)),
  times=10
)


microbenchmark::microbenchmark(
  my_obj_func_3(c(kappa0,tau0)),
  my_log_obj_func(c(kappa0,tau0)),
  times=10
)


starttime <- Sys.time()
o1<-optim(par=c(kappa0,tau0),my_obj_func_3,control=list(maxit=50000))
endtime <- Sys.time()
#times_score_rep[i_n]<-difftime(endtime,starttime, units="secs")
print(difftime(endtime,starttime, units="secs"))
#o_score_list_rep[[i_n]]<-o1

starttime <- Sys.time()
o2<-optim(par=c(kappa0,tau0),my_log_obj_func,control=list(maxit=50000))
endtime <- Sys.time()
#times_log_rep[i_n]<-difftime(endtime,starttime, units="secs")
print(difftime(endtime,starttime, units="secs"))
#o_log_list_rep[[i_n]]<-o2

starttime <- Sys.time()
o3<-optim(par=c(kappa0,tau0),my_log_score_obj_func,control=list(maxit=50000))
endtime <- Sys.time()
#times_log_score_rep[i_n]<-difftime(endtime,starttime, units="secs")
print(difftime(endtime,starttime, units="secs"))
#o_log_score_list_rep[[i_n]]<-o3

starttime <- Sys.time()
o4<-optim(par=c(kappa0,tau0),my_obj_func_3,scoretype="crps",control=list(maxit=50000))
endtime <- Sys.time()
#times_score_rep_crps[i_n]<-difftime(endtime,starttime, units="secs")
print(difftime(endtime,starttime, units="secs"))
#o_score_list_rep_crps[[i_n]]<-o4

starttime <- Sys.time()
o5<-optim(par=c(kappa0,tau0),my_obj_func_3,scoretype="scrps",control=list(maxit=50000))
endtime <- Sys.time()
#times_score_rep_scrps[i_n]<-difftime(endtime,starttime, units="secs")
print(difftime(endtime,starttime, units="secs"))
#o_score_list_rep_scrps[[i_n]]<-o5




#######plot#############

plot(c(10,100,500,1000,2000,3000,4000,5000,6000,10000),(c(0.010,0.016,0.25,0.84,3.73,10.62,21.58,41.33,64.04,NaN))/c(10,100,500,1000,2000,3000,4000,5000,6000,10000))
lines(c(10,100,500,1000,2000,3000,4000,5000,6000,10000),(c(0.009,0.0014,0.23,0.62,2.03,4.43,6.02,9.17,11.74,31.21))/c(10,100,500,1000,2000,3000,4000,5000,6000,10000))

time.sroot <- c(0.009,0.0014,0.23,0.62,2.03,4.43,6.02,9.17,11.74,31.21)
time.ll <- c(0.010,0.016,0.25,0.84,3.73,10.62,21.58,41.33,64.04,NaN)
n.dim.list <- c(10,100,500,1000,2000,3000,4000,5000,6000,10000)
time.df <- data.frame(n=rep(n.dim.list,times=2),t=c(time.sroot,time.ll), type = rep(c("Sroot","LL"),each=length(n.dim.list)))

log(time.sroot[5]/time.sroot[9])/log(n.dim.list[5]/n.dim.list[9]) #1.58
log(time.ll[5]/time.ll[9])/log(n.dim.list[5]/n.dim.list[9]) #2.23

time.df%>%ggplot(aes(x=n,y=t,color=type))+geom_point()
time.df%>%ggplot(aes(x=log(n),y=log(t),color=type,shape=type))+geom_point()+
  geom_smooth(data=subset(time.df,log(n)>5),method = "lm", se = TRUE)+scale_color_brewer(palette="Dark2")
ggsave('matern_regular_grid.pdf',dpi = 1200,width = 12,height = 10,units = 'cm')

subset(time.df,log(n)>5)%>%ggplot(aes(x=log(n),y=log(t),color=type,shape=type))+geom_point()+
  geom_smooth(method = "lm", se = TRUE)+scale_color_brewer(palette="Dark2")
ggsave('matern_regular_grid_cut.pdf',dpi = 1200,width = 12,height = 10,units = 'cm')


plot(log(c(10,100,500,1000,2000,3000,4000,5000,6000,10000)),log(c(0.010,0.016,0.25,0.84,3.73,10.62,21.58,41.33,64.04,NaN)),col="blue")
points(log(c(10,100,500,1000,2000,3000,4000,5000,6000,10000)),log(c(0.009,0.0014,0.23,0.62,2.03,4.43,6.02,9.17,11.74,31.21)),col="red")
#lines(c(3:9),c(3:9)-8)
lines(c(3:10),1.58*(c(3:10))-11.25,col="red")
lines(c(3:10),2.58*(c(3:10))-18.25,col="blue")
axis(side = 1, labels="log(n)")
axis(side = 2, labels="log(time)")#lines(c(3:9),1.5*(c(3:9)-8)+1)
#lines(c(3:9),2*(c(3:9)-8)+2)
#lines(c(3:9),2.4*(c(3:9)-8)+2.5)


plot((c(0.009,0.0014,0.23,0.62,2.03,4.43,6.02,9.17,11.74))/(c(0.010,0.016,0.25,0.84,3.73,10.62,21.58,41.33,64.04)))



############################## Space time

n_dim_1 <- 100
n_dim_2 <- 100
Q<-kronecker(getQ_regular(kappa,tau,n_dim_1,h),getQ_regular(kappa,tau,n_dim_2,h)) #getQ_regular(kappa,tau,n_dim,h)
mu <- rep(0,nrow(Q))
m<-t(inla.qsample(n=1, Q = Q, mu=mu))

nnzero(Q)
sum(spam::bandwidth(spam::as.spam.dgCMatrix(Q))*100*100)
image(Q[1:200,1:1000])

microbenchmark::microbenchmark(
  loo_score_vectorised(m,mu,Q,score="sroot"),
  times=100)

microbenchmark::microbenchmark(
  loo_score_vectorised(m,mu,Q,score="sroot"),
  loo_log_score(m,mu,Q),
  -log_dmvn(m,mu,Q),
  times=100)
