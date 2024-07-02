#time_res<-time_comparison(floor(2^(c(3:24)/2)))
time_res<-time_comparison(c(1:15)*400)
p1.time.comparison<-ggplot(data=time_res,aes(x=log(n),y=log(mean),color=expr))+geom_point()+
  geom_smooth(data=subset(time_res,log(n)>7),method = "lm", se = TRUE)+scale_color_brewer(labels=c("Sr","LL"),palette="Dark2")#+geom_line(aes(x=log(n),y=log(lq)))+geom_line(aes(x=log(n),y=log(uq)))


# time_res_2<-time_comparison(floor(2^(c(3:12)/2)),nlist2=floor(2^(c(3:12)/2)))
time_res_2<-time_comparison(floor(sqrt(c(1:20)*200)),nlist2=floor(sqrt(c(1:20)*200)))
p2.time.comparison<-ggplot(data=subset(time_res_2,log(n)>7),aes(x=log(n),y=log(mean),color=expr))+geom_point()+
  geom_smooth(method = "lm", se = TRUE)+scale_color_brewer(labels=c("Sr","LL"),palette="Dark2")#+geom_line(aes(x=log(n),y=log(lq)))+geom_line(aes(x=log(n),y=log(uq)))



time_res_3<-time_comparison(c(1:15)*400,Qsparse=FALSE)
p3.time.comparison<-ggplot(data=time_res_3,aes(x=log(n),y=log(mean),color=expr))+geom_point()+
  geom_smooth(data=subset(time_res_3,log(n)>7.5),method = "lm", se = TRUE)+scale_color_brewer(labels=c("Sr","LL"),palette="Dark2")

#p.time.comparison <- plot_grid_2(p1.time.comparison,p2.time.comparison)
p.time.comparison <- plot_grid_3(p1.time.comparison,p2.time.comparison,p3.time.comparison)
ggsave('time_comparison_3.pdf',p.time.comparison,dpi = 1200,width = 18,height = 10,units = 'cm')

save(time_res,time_res_3,file="time_comparison_regular_20240624_once_nlarge.Rda")

time_res$model <- "sparse"
#time_res_2$model <- "separable"
time_res_3$model <- "dense"
#time_res_all<-rbind(time_res,time_res_2,time_res_3)
time_res_all<-rbind(time_res,time_res_3)
levels(time_res_all$expr)<-c("Sr","LL")
time_res_all$group <- as.character(time_res_all$expr)
time_res_all$group <- sapply(c(1:nrow(time_res_all)), function(i) paste(time_res_all[i,]$model,time_res_all[i,]$group))
# p.time.comparison.all<-ggplot(data=subset(time_res_all,log(n)>7),aes(x=log(n),y=log(mean),color=group,shape=group))+geom_point()+
#   geom_smooth(method = "lm", se = TRUE, aes(fill=group))+scale_fill_brewer(palette="Paired")+scale_color_brewer(palette="Paired")+scale_shape_manual(values=c(19,17,19,17,19,17))
p.time.comparison.all<-ggplot(data=subset(time_res_all,log(n)>7.5),aes(x=log(n),y=log(mean),color=group,shape=group))+geom_point()+
  geom_smooth(method = "lm", se = TRUE, aes(fill=group))+scale_fill_manual(values=c("#DCE319FF","#95D840FF","#3CBB75FF","#1F968BFF","#2D708EFF","#39568CFF"))+scale_color_manual(values=c("#DCE319FF","#B8DE29FF","#29AF7FFF","#20A387FF","#2D708EFF","#33638DFF"))+scale_shape_manual(values=c(19,19,17,17,19,17))
p.time.comparison.all
ggsave('time_comparison_all_20240624.pdf',p.time.comparison.all,dpi = 1200,width = 18,height = 10,units = 'cm')




ggplot(data=time_res_all,aes(x=log(n),y=log(mean),color=group,shape=group))+geom_point()+
  geom_smooth(data=subset(time_res_all,log(n)>7),method = "lm", se = TRUE, aes(fill=group,linetype=group),alpha=0.3)+scale_fill_manual(values=c("#95D840FF","#2D708EFF","#3CBB75FF","#39568CFF"))+scale_color_manual(values=c("#95D840FF","#2D708EFF","#3CBB75FF","#39568CFF"))+scale_shape_manual(values=c(19,19,17,17,19,17))+scale_linetype_manual(values=c("solid","dotted","solid","dotted"))


subset(time_res,log(n)>7.5) %>%
  group_by(expr) %>%
  do({
    mod = lm(log(mean) ~ log(n), data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })


subset(time_res_2,log(n)>7) %>%
  group_by(expr) %>%
  do({
    mod = lm(log(mean) ~ log(n), data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })


subset(time_res_3,log(n)>8) %>%
  group_by(expr) %>%
  do({
    mod = lm(log(mean) ~ log(n), data = .)
    # print(confint(mod))
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })

