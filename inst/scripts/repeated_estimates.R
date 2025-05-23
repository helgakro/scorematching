library(INLA)
library(grid)
library(cowplot)
library(latex2exp)

set.seed(111101)
# Generate mesh
# Generate points
# Generate A
# Generate test points
# Define parameters
tmptitle <- paste("repeated_estimates_",format(Sys.time(), "%Y%m%d_%H%M%S"),sep = "")
npoints <- 100
sim_loc = matrix(runif(2*npoints,min=0,max=5), ncol = 2, byrow = T)
test_loc = matrix(runif(2*npoints,min=0,max=5), ncol = 2, byrow = T)
mesh_sim = inla.mesh.2d(loc = rbind(matrix(c(0,0,5,5, 0, 5, 5, 0), nrow = 4, byrow = T)), max.edge=c(0.5, 1))
plot(mesh_sim)
points(sim_loc,col='red',pch=19)
mesh_sim$n
A <- inla.spde.make.A(mesh=mesh_sim,loc=sim_loc)
Atest <- inla.spde.make.A(mesh=mesh_sim,loc=test_loc)

spde = inla.spde2.matern(mesh_sim, alpha = 2)
sigma_val <- 0.5 #0.05 #0.1
params_true=spde$param.inla$theta.initial #-1.8301722  0.5646601
#params_true=log(c(1.25,sqrt(gamma(1)/(2^2*1.25^2*4*pi*gamma(2)))))
# params_true=log(c(sqrt(2),sqrt(gamma(1)/(sqrt(2)^2*1^2*4*pi*gamma(2)))))
#params_true=log(c(sqrt(2),1.75))
print(params_true)

# plot mesh
library(ggplot2)
library(inlabru)
Amap.df <- data.frame(lon=c(sim_loc[,1],test_loc[,1]),lat=c(sim_loc[,2],test_loc[,2]),type=rep(c("Train","Test"),each=50))
p.A.map<-ggplot(Amap.df)+gg(mesh_sim,edge.color = "gray",
                            int.color = "black",
                            ext.color = "black")+geom_point(aes(x=lon,y=lat,shape=type),size=2)+scale_color_viridis_d()+xlab("longitude")+ylab("latitude")+ guides(shape="none")
p.A.map
#ggsave('maternA_map_all.pdf',p.A.map,dpi = 1200,width = 12,height = 10,units = 'cm')
ggsave(paste(tmptitle,"_map.pdf",sep=""),p.A.map,dpi = 1200,width = 12,height = 10,units = 'cm')

#range
sqrt(8)/exp(params_true[2])

#sigma
sqrt(gamma(1)/(exp(params_true[2])^2*exp(params_true[1])^2*(4*pi)*gamma(2)))


#params_true<-log(c(sigma_val,0.5,1.75))#params_true<-c(log(sigma_val),params_true)

# M2

Q = inla.spde.precision(spde, theta=params_true)
n_rep=300 #100
res_no_outliers <- repeated_inference(spde,mesh_sim$n,n_rep,Q) #no outliers
res_5_outliers <- repeated_inference(spde,mesh_sim$n,n_rep,Q,5,5) #5 outliers (50%)
res_10_outliers <- repeated_inference(spde,mesh_sim$n,n_rep,Q,10,5) #10 outliers (100%)

save(res_no_outliers,res_5_outliers,res_10_outliers,file=paste(tmptitle,"resM2",".Rda",sep=""))

#params_true<-log(c(sigma_val,0.5,1.75))#params_true<-c(log(sigma_val),params_true)
#params_true<-params_true[c(2,3)]
p.res_no_outliers <- plot_results(res_no_outliers,extended = TRUE)
p.res_5_outliers <- plot_results(res_5_outliers,extended = TRUE)
p.res_10_outliers <- plot_results(res_10_outliers,extended = TRUE)

ggsave(paste(tmptitle,"GMRF_time_hist_est_no_rep100.pdf",sep=""),p.res_no_outliers$p.time.hist, dpi=1200,width=18,height=8,unit="cm")
ggsave(paste(tmptitle,paste("GMRF_scatter_est_no_5_10","_",n_rep,"rep",".pdf",sep=""),sep=""),plot_grid_3(p.res_no_outliers$p.scatter,p.res_5_outliers$p.scatter,p.res_10_outliers$p.scatter), dpi=1200,width=18,height=8,unit="cm")
ggsave(paste(tmptitle,"GMRF_p1_est_no_5_10.pdf",sep=""),plot_grid_3(p.res_no_outliers$p.hist.p1,p.res_5_outliers$p.hist.p1,p.res_10_outliers$p.hist.p1), dpi=1200,width=18,height=8,unit="cm")
ggsave(paste(tmptitle,"GMRF_p2_est_no_5_10.pdf",sep=""),plot_grid_3(p.res_no_outliers$p.hist.p2,p.res_5_outliers$p.hist.p2,p.res_10_outliers$p.hist.p2), dpi=1200,width=18,height=8,unit="cm")
ggsave(paste(tmptitle,"GMRF_time_est_no_5_10.pdf",sep=""),plot_grid_3(p.res_no_outliers$p.time,p.res_5_outliers$p.time,p.res_10_outliers$p.time), dpi=1200,width=18,height=8,unit="cm")
ggsave(paste(tmptitle,"GMRF_time_hist_est_no_5_10.pdf",sep=""),plot_grid_3(p.res_no_outliers$p.time.hist+xlim(0,2)+ylim(0,61),p.res_5_outliers$p.time.hist+xlim(0,2)+ylim(0,61),p.res_10_outliers$p.time.hist+xlim(0,2)+ylim(0,61)), dpi=1200,width=18,height=8,unit="cm")

df0<-p.res_no_outliers$df
df5<-p.res_5_outliers$df
df10<-p.res_10_outliers$df

df0$outlier <- "0"
df5$outlier <- "5"
df10$outlier <- "10"
df.all <- rbind(df0,df5,df10)

df.all%>%group_by(outlier,method)%>%summarize(
  mean = mean(run.time, na.rm = TRUE),
  sd = sd(run.time, na.rm = TRUE)
)

df.all$method <- factor(df.all$method, levels=c('LL', 'Slog', 'SCRPS', 'Sroot','CRPS','rCRPS'))
df.all$outlier <- factor(df.all$outlier,levels=c("0","5","10"))

p.par.box_outliers<-plot_grid_2(ggplot(df.all,aes(x=outlier, y=par.2,fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = params_true[2])+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\log(\\kappa)$"))+xlab("Number of outliers")+
                                  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)),
                                ggplot(df.all,aes(x=outlier, y=par.1,fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = params_true[1])+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\log(\\tau)$"))+xlab("Number of outliers")+
                                  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)))
ggsave(paste(tmptitle,"resM2_log",".pdf",sep=""),p.par.box_outliers,dpi = 1200,width = 18,height = 8,units = 'cm')
p.par.box_outliers<-plot_grid_2(ggplot(df.all,aes(x=outlier, y=exp(par.2),fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = exp(params_true[2]))+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\kappa$"))+xlab("Number of outliers")+guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)),
                                ggplot(df.all,aes(x=outlier, y=exp(par.1),fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = exp(params_true[1]))+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\tau$"))+xlab("Number of outliers")+
                                  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)))
p.par.box_outliers
ggsave(paste(tmptitle,"resM2",".pdf",sep=""),p.par.box_outliers,dpi = 1200,width = 18,height = 8,units = 'cm')


# M3

res_0_outliers_nresp <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 0, outlier_val = 5,sigma_val=sigma_val,A=A,Atest=Atest,scoretypes=c("sroot","ll","slog","crps","scrps","rcrps") ) #no outliers 0.0002
res_10_outliers_nresp <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 10, outlier_val = 5,sigma_val=sigma_val,A=A,Atest=Atest,scoretypes=c("sroot","ll","slog","crps","scrps","rcrps") ) #no outliers 0.0002
res_10_outliers_nresp_large <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 10, outlier_val = 10,sigma_val=sigma_val,A=A,Atest=Atest,scoretypes=c("sroot","ll","slog","crps","scrps","rcrps") ) #no outliers 0.0002

save(res_0_outliers_nresp,res_10_outliers_nresp,file=paste(tmptitle,"resM3",".Rda",sep=""))

params_true <- c(log(sigma_val), params_true)
p.res_0_outliers_nresp <- plot_results(res_0_outliers_nresp,extended = TRUE)
p.res_10_outliers_nresp <- plot_results(res_10_outliers_nresp,extended = TRUE)
df0<-p.res_0_outliers_nresp$df
df10<-p.res_10_outliers_nresp$df

df0$outlier <- "0"
df10$outlier <- "10"
df.all <- rbind(df0,df10)

df.all%>%group_by(outlier,method)%>%summarize(
  mean = mean(run.time, na.rm = TRUE),
  sd = sd(run.time, na.rm = TRUE)
)

df.all$method <- factor(df.all$method, levels=c('LL', 'Slog', 'SCRPS', 'Sroot','CRPS','rCRPS'))
df.all$outlier <- factor(df.all$outlier,levels=c("0","10"))

p.par.box_outliers<-plot_grid_3(ggplot(df.all,aes(x=outlier, y=par.1,fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = params_true[1])+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\log(\\sigma_{\\epsilon})$"))+xlab("Number of outliers")+
                                  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)),
                                ggplot(df.all,aes(x=outlier, y=par.3,fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = params_true[2])+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\log(\\kappa)$"))+xlab("Number of outliers")+
                                  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)),
                                ggplot(df.all,aes(x=outlier, y=par.2,fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = params_true[3])+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\log(\\tau)$"))+xlab("Number of outliers")+
                                  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)))
ggsave(paste(tmptitle,"resM3_log",".pdf",sep=""),p.par.box_outliers,dpi = 1200,width = 18,height = 8,units = 'cm')

p.par.box_outliers<-plot_grid_3(ggplot(df.all,aes(x=outlier, y=exp(par.1),fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = exp(params_true[1]))+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\sigma_{\\epsilon}$"))+xlab("Number of outliers")+guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)),
                                ggplot(df.all,aes(x=outlier, y=exp(par.3),fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = exp(params_true[3]))+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\kappa$"))+xlab("Number of outliers")+guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)),
                                ggplot(subset(df.all,1/exp(par.2)<20),aes(x=outlier, y=exp(par.2),fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = exp(params_true[2]))+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\tau$"))+xlab("Number of outliers")+guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)))

p.par.box_outliers
ggsave(paste(tmptitle,"resM3_filtered",".pdf",sep=""),p.par.box_outliers,dpi = 1200,width = 18,height = 8,units = 'cm')

# M4 nonstationary


res_0_outliers_nresp_ns <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 0, outlier_val = 5,sigma_val=sigma_val,A=t(t(A)*mesh_sim$loc[,1]),Atest=t(t(Atest)*mesh_sim$loc[,1]),scoretypes=c("sroot","ll","slog","crps","scrps","rcrps")) #no outliers 0.0002

res_10_outliers_nresp_ns <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 10, outlier_val = 5,sigma_val=sigma_val,A=t(t(A)*mesh_sim$loc[,1]),Atest=t(t(Atest)*mesh_sim$loc[,1]),scoretypes=c("sroot","ll","slog","crps","scrps","rcrps") ) #no outliers 0.0002
#res_10_outliers_nresp <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 10, outlier_val = 5,sigma_val=0.5,A=t(t(A)*mesh_sim$loc[,1]),Atest=t(t(Atest)*mesh_sim$loc[,1]),scoretypes=c("sroot","ll","slog","crps","scrps") ) #no outliers 0.0002
save(res_0_outliers_nresp_ns,res_10_outliers_nresp_ns,file=paste(tmptitle,"resM4_ns",".Rda",sep=""))
p.res_0_outliers_nresp_ns <- plot_results(res_0_outliers_nresp_ns,extended = TRUE)
p.res_10_outliers_nresp_ns <- plot_results(res_10_outliers_nresp_ns,extended = TRUE)
df0<-p.res_0_outliers_nresp_ns$df
df10<-p.res_10_outliers_nresp_ns$df

df0$outlier <- "0"
df10$outlier <- "10"
df.all <- rbind(df0,df10)

df.all$method <- factor(df.all$method, levels=c('LL', 'Slog', 'SCRPS', 'Sroot','CRPS','rCRPS'))
df.all$outlier <- factor(df.all$outlier,levels=c("0","10"))

p.par.box_outliers<-plot_grid_3(ggplot(df.all,aes(x=outlier, y=par.1,fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = params_true[1])+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\log(\\sigma_{\\epsilon})$"))+xlab("Number of outliers")+
                                  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)),
                                ggplot(df.all,aes(x=outlier, y=par.3,fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = params_true[3])+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\log(\\kappa)$"))+xlab("Number of outliers")+
                                  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)),
                                ggplot(df.all,aes(x=outlier, y=par.2,fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = params_true[2])+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\log(\\tau)$"))+xlab("Number of outliers")+
                                  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)))
ggsave(paste(tmptitle,"resM4_log",".pdf",sep=""),p.par.box_outliers,dpi = 1200,width = 18,height = 8,units = 'cm')

p.par.box_outliers<-plot_grid_3(ggplot(df.all,aes(x=outlier, y=exp(par.1),fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = exp(params_true[1]))+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\sigma_{\\epsilon}$"))+xlab("Number of outliers")+guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)),
                                ggplot(df.all,aes(x=outlier, y=exp(par.3),fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = exp(params_true[3]))+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\kappa$"))+xlab("Number of outliers")+guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)),
                                ggplot(subset(df.all,1/exp(par.3)<20),aes(x=outlier, y=exp(par.2),fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = exp(params_true[2]))+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\tau$"))+xlab("Number of outliers")+guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)))

p.par.box_outliers
ggsave(paste(tmptitle,"resM4_filtered",".pdf",sep=""),p.par.box_outliers,dpi = 1200,width = 18,height = 8,units = 'cm')

###xy spatial

res_0_outliers_nresp_xy <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 0, outlier_val = 5,sigma_val=sigma_val,A=t(t(A)*sqrt(abs(mesh_sim$loc[,1]*mesh_sim$loc[,2]))),Atest=t(t(Atest)*mesh_sim$loc[,1]),scoretypes=c("sroot","ll","slog","crps","scrps","rcrps") ) #no outliers 0.0002

res_10_outliers_nresp_xy <- repeated_inference_norm_resp(spde,mesh_sim$n,n_rep,Q,n_outlier = 10, outlier_val = 5,sigma_val=sigma_val,A=t(t(A)*sqrt(abs(mesh_sim$loc[,1]*mesh_sim$loc[,2]))),Atest=t(t(Atest)*mesh_sim$loc[,1]),scoretypes=c("sroot","ll","slog","crps","scrps","rcrps") ) #no outliers 0.0002
save(res_0_outliers_nresp_xy,res_10_outliers_nresp_xy,file=paste(tmptitle,"resM4_xy",".Rda",sep=""))

p.res_0_outliers_nresp_xy <- plot_results(res_0_outliers_nresp_xy,extended = TRUE)
p.res_10_outliers_nresp_xy <- plot_results(res_10_outliers_nresp_xy,extended = TRUE)
df0<-p.res_0_outliers_nresp_xy$df
df10<-p.res_10_outliers_nresp_xy$df

df0$outlier <- "0"
df10$outlier <- "10"
df.all <- rbind(df0,df10)

df.all$method <- factor(df.all$method, levels=c('LL', 'Slog', 'SCRPS', 'Sroot','CRPS','rCRPS'))
df.all$outlier <- factor(df.all$outlier,levels=c("0","10"))

p.par.box_outliers<-plot_grid_3(ggplot(df.all,aes(x=outlier, y=par.1,fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = params_true[1])+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\log(\\sigma_{\\epsilon})$"))+xlab("Number of outliers")+
                                  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)),
                                ggplot(df.all,aes(x=outlier, y=par.3,fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = params_true[3])+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\log(\\kappa)$"))+xlab("Number of outliers")+
                                  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)),
                                ggplot(df.all,aes(x=outlier, y=par.2,fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = params_true[2])+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\log(\\tau)$"))+xlab("Number of outliers")+
                                  guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)))
ggsave(paste(tmptitle,"resM4xy_log",".pdf",sep=""),p.par.box_outliers,dpi = 1200,width = 18,height = 8,units = 'cm')

p.par.box_outliers<-plot_grid_3(ggplot(df.all,aes(x=outlier, y=exp(par.1),fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = exp(params_true[1]))+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\sigma_{\\epsilon}$"))+xlab("Number of outliers")+guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)),
                                ggplot(df.all,aes(x=outlier, y=exp(par.3),fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = exp(params_true[3]))+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\kappa$"))+xlab("Number of outliers")+guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)),
                                ggplot(subset(df.all,1/exp(par.3)<20),aes(x=outlier, y=exp(par.2),fill=method,color=method))+geom_boxplot(alpha=0.5)+geom_hline(yintercept = exp(params_true[2]))+
                                  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+ylab(TeX("$\\tau$"))+xlab("Number of outliers")+guides(colour = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)))

p.par.box_outliers
ggsave(paste(tmptitle,"resM4xy_filtered",".pdf",sep=""),p.par.box_outliers,dpi = 1200,width = 18,height = 8,units = 'cm')




#### traintest


score_sroot <- sapply(res_0_outliers_nresp$o_sroot[1:300], function(o) o$value)
score_ll <- res_0_outliers_nresp$score_ll

score_res_nresp_df_3 <- data.frame(val.sroot=score_sroot,val.ll=score_ll,outlier="no")
p.score.hist3<-ggplot(score_res_nresp_df_3,aes(x=(val.ll-val.sroot)/val.ll*100,color="#1B9E77",fill="#1B9E77"))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+theme(legend.position = "none")+xlab("Rel score diff")


score_res_nresp_df_4 <- data.frame(val.sroot=res_0_outliers_nresp$pred_sroot,val.ll=res_0_outliers_nresp$pred_ll,outlier="no")
p.score.hist4<-ggplot(score_res_nresp_df_4,aes(x=(val.ll-val.sroot)/val.ll*100,color="#1B9E77",fill="#1B9E77"))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+theme(legend.position = "none")+xlab("Rel score diff")

score_res_nresp_df_5 <- data.frame(val.sroot=res_0_outliers_nresp$rmse_sroot,val.ll=res_0_outliers_nresp$rmse_ll,outlier="no")
p.score.hist5<-ggplot(score_res_nresp_df_5,aes(x=(val.ll-val.sroot)/val.ll*100,color="#1B9E77",fill="#1B9E77"))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+theme(legend.position = "none")+xlab("Rel score diff")

score_res_nresp_df_6 <- data.frame(val.sroot=res_0_outliers_nresp$rmse_sroot_train,val.ll=res_0_outliers_nresp$rmse_ll_train,outlier="no")
p.score.hist6<-ggplot(score_res_nresp_df_6,aes(x=(val.ll-val.sroot)/val.ll*100,color="#1B9E77",fill="#1B9E77"))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+theme(legend.position = "none")+xlab("Rel score diff")


# p.score.hist.multi<-plot_grid_3_nolegend(p.score.hist3+xlim(c(-0.75,0.75)),p.score.hist4+xlim(c(-0.75,0.75)),p.score.hist5+xlim(c(-0.75,0.75)))
p.score.hist.multi<-plot_grid_4_nolegend(p.score.hist3+xlim(c(-0.75,0.75)),p.score.hist4+xlim(c(-0.75,0.75)),p.score.hist6+xlim(c(-0.75,0.75)),p.score.hist5+xlim(c(-0.75,0.75)))
p.score.hist.multi
ggsave(paste(tmptitle,"resM3_traintest_ll",".pdf",sep=""),p.score.hist.multi,dpi = 1200,width = 18,height = 8,units = 'cm')
ggsave(paste(tmptitle,"resM3_traintest_ll",".pdf",sep=""),p.score.hist.multi,dpi = 1200,width = 18,height = 12,units = 'cm')




score_sroot <- sapply(res_10_outliers_nresp$o_sroot[1:n_rep], function(o) o$value)
score_ll <- res_10_outliers_nresp$score_ll

score_res_nresp_df_3 <- rbind(score_res_nresp_df_3,data.frame(val.sroot=score_sroot,val.ll=score_ll,outlier="medium"))
p.score.hist3<-ggplot(score_res_nresp_df_3,aes(x=(val.ll-val.sroot)/val.sroot,color=outlier,fill=outlier))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+theme(legend.position = "none")


score_res_nresp_df_4 <- rbind(score_res_nresp_df_4,data.frame(val.sroot=res_10_outliers_nresp$pred_sroot,val.ll=res_10_outliers_nresp$pred_ll,outlier="medium"))
p.score.hist4<-ggplot(score_res_nresp_df_4,aes(x=(val.ll-val.sroot)/val.sroot,color=outlier,fill=outlier))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+theme(legend.position = "none")

score_res_nresp_df_5 <- rbind(score_res_nresp_df_5,data.frame(val.sroot=res_10_outliers_nresp$rmse_sroot,val.ll=res_10_outliers_nresp$rmse_ll,outlier="medium"))
p.score.hist5<-ggplot(score_res_nresp_df_5,aes(x=(val.ll-val.sroot)/val.sroot,color=outlier,fill=outlier))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+theme(legend.position = "none")

score_res_nresp_df_6 <- rbind(score_res_nresp_df_6,data.frame(val.sroot=res_10_outliers_nresp$rmse_sroot_train,val.ll=res_10_outliers_nresp$rmse_ll_train,outlier="medium"))
p.score.hist6<-ggplot(score_res_nresp_df_6,aes(x=(val.ll-val.sroot)/val.sroot,color=outlier,fill=outlier))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+theme(legend.position = "none")


score_sroot <- sapply(res_10_outliers_nresp_large$o_sroot[1:n_rep], function(o) o$value)
score_ll <- res_10_outliers_nresp_large$score_ll

score_res_nresp_df_3 <- rbind(score_res_nresp_df_3,data.frame(val.sroot=score_sroot,val.ll=score_ll,outlier="large"))
score_res_nresp_df_3$outlier <- as.factor(score_res_nresp_df_3$outlier)
score_res_nresp_df_3$outlier <-  factor(score_res_nresp_df_3$outlier,levels=c("no","medium","large"))
p.score.hist3<-ggplot(score_res_nresp_df_3,aes(x=(val.ll-val.sroot)/val.sroot,color=outlier,fill=outlier))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+xlab("Rel score diff")


score_res_nresp_df_4 <- rbind(score_res_nresp_df_4,data.frame(val.sroot=res_10_outliers_nresp_large$pred_sroot,val.ll=res_10_outliers_nresp_large$pred_ll,outlier="large"))
score_res_nresp_df_4$outlier <- as.factor(score_res_nresp_df_4$outlier)
score_res_nresp_df_4$outlier <-  factor(score_res_nresp_df_4$outlier,levels=c("no","medium","large"))
p.score.hist4<-ggplot(score_res_nresp_df_4,aes(x=(val.ll-val.sroot)/val.sroot,color=outlier,fill=outlier))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+xlab("Rel score diff")

score_res_nresp_df_5 <- rbind(score_res_nresp_df_5,data.frame(val.sroot=res_10_outliers_nresp_large$rmse_sroot,val.ll=res_10_outliers_nresp_large$rmse_ll,outlier="large"))
score_res_nresp_df_5$outlier <- as.factor(score_res_nresp_df_5$outlier)
score_res_nresp_df_5$outlier <-  factor(score_res_nresp_df_5$outlier,levels=c("no","medium","large"))
p.score.hist5<-ggplot(score_res_nresp_df_5,aes(x=(val.ll-val.sroot)/val.sroot,color=outlier,fill=outlier))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+xlab("Rel score diff")


score_res_nresp_df_6 <- rbind(score_res_nresp_df_6,data.frame(val.sroot=res_10_outliers_nresp_large$rmse_sroot_train,val.ll=res_10_outliers_nresp_large$rmse_ll_train,outlier="large"))
score_res_nresp_df_6$outlier <- as.factor(score_res_nresp_df_6$outlier)
score_res_nresp_df_6$outlier <-  factor(score_res_nresp_df_6$outlier,levels=c("no","medium","large"))
p.score.hist6<-ggplot(score_res_nresp_df_6,aes(x=(val.ll-val.sroot)/val.sroot,color=outlier,fill=outlier))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+xlab("Rel score diff")



#p.score.hist.multi<-plot_grid_3(p.score.hist3,p.score.hist4,p.score.hist5)
p.score.hist.multi<-plot_grid_4(p.score.hist3,p.score.hist4,p.score.hist6,p.score.hist5)
p.score.hist.multi
ggsave(paste(tmptitle,"resM3_traintest_outliers",".pdf",sep=""),p.score.hist.multi,dpi = 1200,width = 18,height = 8,units = 'cm')




######

p.score.hist3_normll<-ggplot(score_res_nresp_df_3,aes(x=100*(val.ll-val.sroot)/val.ll,color=outlier,fill=outlier))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+xlab("Rel score diff")
p.score.hist4_normll<-ggplot(score_res_nresp_df_4,aes(x=100*(val.ll-val.sroot)/val.ll,color=outlier,fill=outlier))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+xlab("Rel score diff")
p.score.hist5_normll<-ggplot(score_res_nresp_df_5,aes(x=100*(val.ll-val.sroot)/val.ll,color=outlier,fill=outlier))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+xlab("Rel score diff")
p.score.hist6_normll<-ggplot(score_res_nresp_df_6,aes(x=100*(val.ll-val.sroot)/val.ll,color=outlier,fill=outlier))+geom_histogram(alpha=0.5, position="identity")+
  scale_fill_brewer(palette = "Dark2")+  scale_color_brewer(palette = "Dark2")+xlab("Rel score diff")

# p.score.hist.multi_normll<-plot_grid_3(p.score.hist3_normll+xlim(c(-1,12)),p.score.hist4_normll+xlim(c(-1,12)),p.score.hist5_normll+xlim(c(-1,12)))
p.score.hist.multi_normll<-plot_grid_4(p.score.hist3_normll+xlim(c(-1,12)),p.score.hist4_normll+xlim(c(-1,12)),p.score.hist6_normll+xlim(c(-1,12)),p.score.hist5_normll+xlim(c(-1,12)))
p.score.hist.multi_normll
ggsave(paste(tmptitle,"resM3_traintest_outliers_normll",".pdf",sep=""),p.score.hist.multi_normll,dpi = 1200,width = 18,height = 12,units = 'cm')
