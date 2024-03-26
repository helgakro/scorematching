library(INLA)
library(grid)
library(cowplot)

sim_loc = matrix(c(0,0,5,5, 0, 5, 5, 0), nrow = 4, byrow = T)
mesh_sim = inla.mesh.2d(loc = sim_loc, max.edge=c(0.5, 1))
plot(mesh_sim)
points(sim_loc,col='red',pch=19)
mesh_sim$n

spde = inla.spde2.matern(mesh_sim, alpha = 2)
params_true=spde$param.inla$theta.initial
print(params_true)
Q = inla.spde.precision(spde, theta=spde$param.inla$theta.initial)


# n_mesh= mesh_sim$n
n_rep=100
res_no_outliers <- repeated_inference(mesh_sim$n,n_rep,Q) #no outliers
res_50_outliers <- repeated_inference(mesh_sim$n,n_rep,Q,50,4) #5 outliers (50%)
res_100_outliers <- repeated_inference(mesh_sim$n,n_rep,Q,100,4) #10 outliers (100%)
res_25_outliers <- repeated_inference(mesh_sim$n,n_rep,Q,25,4) #5 outliers (25%)
res_25_outliers <- repeated_inference(mesh_sim$n,n_rep,Q,75,4) #5 outliers (75%)

p.res_no_outliers <- plot_results(res_no_outliers)
p.res_5_outliers <- plot_results(res_5_outliers)
p.res_10_outliers <- plot_results(res_10_outliers)


########################### Save results #######################

ggsave(paste("GMRF_scatter_est_no_5_10","_",n_rep,"rep",".pdf",sep=""),plot_grid_3(p.res_no_outliers$p.scatter,p.res_5_outliers$p.scatter,p.res_10_outliers$p.scatter), dpi=1200,width=18,height=8,unit="cm")
ggsave("GMRF_mu_est_no_5_10.pdf",plot_grid_3(p.res_no_outliers$p.hist.mu,p.res_5_outliers$p.hist.mu,p.res_10_outliers$p.hist.mu), dpi=1200,width=18,height=8,unit="cm")
ggsave("GMRF_rho_est_no_5_10.pdf",plot_grid_3(p.res_no_outliers$p.hist.rho,p.res_5_outliers$p.hist.rho,p.res_10_outliers$p.hist.rho), dpi=1200,width=18,height=8,unit="cm")
ggsave("GMRF_time_est_no_5_10.pdf",plot_grid_3(p.res_no_outliers$p.time,p.res_5_outliers$p.time,p.res_10_outliers$p.time), dpi=1200,width=18,height=8,unit="cm")
ggsave("GMRF_time_hist_est_no_5_10.pdf",plot_grid_3(p.res_no_outliers$p.time.hist+xlim(0,2)+ylim(0,61),p.res_5_outliers$p.time.hist+xlim(0,2)+ylim(0,61),p.res_10_outliers$p.time.hist+xlim(0,2)+ylim(0,61)), dpi=1200,width=18,height=8,unit="cm")


################# normal response model

res_no_outliers_nresp <- repeated_inference_norm_resp(mesh_sim$n,n_rep,Q,sigma_val=0.2) #no outliers
p.res_no_outliers_nresp <- plot_results(res_no_outliers_nresp)
p.res_no_outliers_nresp$p.scatter
p.res_no_outliers_nresp$p.hist.mu
p.res_no_outliers_nresp$p.hist.rho
p.res_no_outliers_nresp$p.time
p.res_no_outliers_nresp$p.time.hist

score_par <- sapply(res_no_outliers_nresp$o_sroot[1:n_res], function(o) o$par)
print(sum(sapply(res_no_outliers_nresp$o_sroot[1:n_res], function(o) o$convergence)))
par_df <- data.frame(method=rep(c("Sroot"),each=n_res), mu = c(score_par[1,]), rho = c(score_par[2,]), run.time=c(times_score_rep[1:n_res]),i=rep(c(1:n_res),1))
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

