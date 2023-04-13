
require(RColorBrewer)
library(latex2exp)
library(cowplot)
library("randomcoloR")
library(latex2exp)
shapes <- c(3,4,6,1,7,8,5)
shapes_K <- c(21:24)
cols=c("#007413","#720044","#e166d1","#74a00d","#f54771","#715500","#017ad3","#c95038")
themet <-  theme_minimal()+   theme(title=element_text(size=15,face="plain"),
                   axis.title=element_text(size=18,face = "plain"), 
                   axis.text = element_text(size = 14),
                   legend.text = element_text(size=14),
                   legend.title=element_text(size=18), aspect.ratio = 1)

themeh <-  theme_minimal()+   theme(title=element_text(size=15,face="plain"),
                                    axis.title=element_text(size=18,face = "plain"), 
                                    axis.text = element_text(size = 14),
                                    legend.text = element_text(size=14),
                                    legend.title=element_text(size=18), aspect.ratio = 1,
                                    legend.direction="horizontal", legend.position = "top")
colsN =brewer.pal(8, "Dark2")
shapesW = c(3, 1,6,4,7,8,5)
# aux functions ----


make_legend <- function(fig, dir){
  if (dir=="H"){
    leg<-guides(shape = guide_legend(nrow = 1))+theme(legend.direction="horizontal", legend.position = "top")}
  else {
    leg<-guides(shape = guide_legend(ncol = 1))
  }
  fig <-ggplot()+
    geom_point(data=datat,aes(x=conn05*P, y=meanvol/2, shape=Pf, colour=Pf),size=3)+
    scale_shape_manual(values=shapes)+
    scale_color_manual(values=cols)+
    labs(y= TeX("\\Omega"), x = TeX(" \\kappa P"), shape="S", color="S:")+
    themet+
    leg
  
  fig
  leg <- get_legend(fig) 
  return(as_ggplot(leg))
}

# read and make stats for niche overlap figure
read_makestats_gamma <- function(){
  data3 <- readData("gamma/it_kappa",3,TRUE,1,10)
  data3e <- readData("gamma/it_kappa_gammaEx",3,TRUE,1,10)
  data3e$vol<-data3e$vol/2
  data3<- rbind(data3,data3e[,1:7])
  stats3 <- getStatsDiff(data3) 
  stats3$S<-"6"
  data5 <- readData("gamma/it_kappa",5,TRUE,1,10)
  data5e <- readData("gamma/it_kappa_gammaEx",5,TRUE,1,10)
  data5e$vol<-data5e$vol/2
  data5<- rbind(data5,data5e[,1:7])
  stats5 <- getStatsDiff(data5) 
  stats5$S<-"10"
  data10 <- readData("gamma/it_kappa",10,TRUE,1,10)
  data10e <- readData("gamma/it_kappa_gammaEx",10,TRUE,1,10)
  data10e$vol<-data10e$vol/2
  data10<- rbind(data10,data10e[,1:7])
  stats10 <- getStatsDiff(data10) 
  #stats10$meanvol <- 2*stats10$meanvol
  stats10$S<-"20"
  return(list(stats3, stats5, stats10))
}



# MCRM connectance vs volume -----

# kappa vs Omega for p = 6, 10 and 20 with fitted line
Omega_kappa_61020 <- function(datat, fitted_dS, three, fit){
  if (three == T){
    datat<-subset(datat,P==6|P==10|P==20, select=colnames(datat))
    colt <- cols[c(1,3,7)]
    shapt <- shapes[c(1,3,7)]
    fitted_dS<-fitted_dS[fitted_dS$S==6|fitted_dS$S==10|fitted_dS$S==20,]
  }else{
    colt <- cols
    shapt <- shapes 
  }
  if (fit==T){
    fit<-geom_line(data=fitted_dS, aes(x=x, y=y/2, colour=S,), show.legend = FALSE)
  } else {
    fit <- xlim(0,1)
  }
  fig <-ggplot()+
    geom_point(data=datat,aes(x=conn05, y=meanvol/2, shape=Pf, colour=Pf),size=3,show.legend = F)+
    scale_shape_manual(values=shapt)+
    fit+
    #geom_line(data=fitted_dS, aes(x=xS, y=y8))+
    ylim(0,.25)+
    #xlim(min(datat$SC)-.02,20)+
    scale_color_manual(values=colt)+
    labs(y= TeX("\\Omega"), x = TeX(" \\kappa"),  color="S")+
    themet
  return(fig)
}
# kappa vs Omega for all with fitted line
Omega_kappa_all <- function(datat, fitted_dS, ylim, xlim){
  fig <-ggplot()+
    geom_point(data=datat,aes(x=conn05, y=meanvol/2, shape=Pf, colour=Pf),size=3,show.legend = F)+
    scale_shape_manual(values=shapes)+
    geom_line(data=fitted_dS, aes(x=x, y=y/2, colour=S,), show.legend = FALSE)+
    #geom_line(data=fitted_dS, aes(x=xS, y=y8))+
    ylim(0,ylim)+
    xlim(xlim,1)+
    #xlim(min(datat$SC)-.02,20)+
    scale_color_manual(values=cols)+
    labs(y= TeX("\\Omega"), x = TeX(" \\kappa"),  color="S")+
    themet
  return(fig)
}

# volume per species - kappa vs Omega for p = 6, 10 and 20 with fitted line

Omega_per_sp_kappa <- function(datat){
  fig <-ggplot()+
    geom_point(data=datat[datat$P==6|datat$P==10|datat$P==20,],aes(x=conn05, y=(meanvol/2)^(1/(P/2)), shape=Pf, colour=Pf),size=3,show.legend = T)+
    geom_line(data=datat[datat$P==6|datat$P==10|datat$P==20,],aes(x=conn05, y=(meanvol/2)^(1/(P/2)), colour=Pf),alpha=0.15, linetype=2,show.legend = F)+
    scale_shape_manual(values=shapes[c(1,3,7)])+
    scale_color_manual(values=cols[c(1,3,7)])+
    labs(y= TeX("\\Omega^(1/n)"), x = TeX(" \\kappa"), shape="S", color="S")+
    themet
  return(fig)
}

# kappa vs Omega*P (data collapse)

Omega_kappaP <- function(datat, fitted){
  fig <-ggplot()+
    geom_point(data=datat,aes(x=conn05*P, y=meanvol/2, shape=Pf, colour=Pf),size=3,show.legend = F)+
    scale_shape_manual(values=shapes)+
    ylim(0,.25)+
    scale_color_manual(values=cols)+
    labs(y= TeX("\\Omega"), x = TeX(" \\kappa P"),  color="S")+
  themet+
    geom_line(data=fitted, aes(x=as.numeric(S)*x, y=y/2),color="black", show.legend = FALSE)
  return(fig)
}

Omega_kappa_all <- function(datat){
  fig <-ggplot()+
    geom_point(data=datat,aes(x=conn05, y=meanvol/2, shape=Pf, colour=Pf),size=3,show.legend = T)+
    geom_line(data=datat,aes(x=conn05, y=meanvol/2, colour=Pf),linetype=2,alpha=0.15,show.legend = T)+
    scale_shape_manual(values=shapes)+
    ylim(0,.25)+
    scale_color_manual(values=cols)+
    labs(y= TeX("\\Omega"), x = TeX(" \\kappa"),  color="S", shape="S")+
    themet#+
    #geom_line(data=fitted, aes(x=as.numeric(S)*x, y=y/2),color="black", show.legend = FALSE)
  return(fig)
}

logloga_all <- function(datat){
  fig <-ggplot()+
    geom_point(data=datat,aes(x=logSC, y=logV, shape=Pf, colour=Pf),size=3,show.legend = T)+
    geom_line(data=datat,aes(x=logSC, y=logV, colour=Pf),linetype=2,alpha=0.15,show.legend = T)+
    scale_shape_manual(values=shapes)+
    scale_color_manual(values=cols)+
    labs(y= TeX("log(1/\\Omega)"), x = TeX("log(\\kappa S)"),  color="S", shape="S")+
    themet#+
  #geom_line(data=fitted, aes(x=as.numeric(S)*x, y=y/2),color="black", show.legend = FALSE)
  return(fig)
}

Omega_kappa_all_f <- function(datat,fit){
  fig <-ggplot()+
    geom_point(data=datat,aes(x=conn05, y=meanvol/2, shape=Pf, colour=Pf),size=3,show.legend = T)+
    scale_shape_manual(values=shapes)+
    ylim(0,.25)+
    scale_color_manual(values=cols)+
    labs(y= TeX("\\Omega"), x = TeX(" \\kappa"),  color="S", shape="S")+
    themet+
    geom_line(data=fit, aes(x=x, y=y/2, colour=S), show.legend = FALSE)
  return(fig)
}

#make legend

make_legend_620 <- function(datat,dir){
  if (dir=="H"){
    leg<-guides(shape = guide_legend(nrow = 1))+theme(legend.direction="horizontal", legend.position = "top")}
  else {
    leg<-guides(shape = guide_legend(ncol = 1))
  }
  fig <-ggplot()+
    geom_point(data=datat,aes(x=conn05*P, y=meanvol/2, shape=Pf, colour=Pf),size=3)+
    scale_shape_manual(values=shapes)+
    scale_color_manual(values=cols)+
    labs(y= TeX("\\Omega"), x = TeX(" \\kappa P"), shape="S", color="S:")+
    themet+
    leg
  
  fig
    leg <- get_legend(fig) 
  return(as_ggplot(leg))
}


# make legend general ----

# Random competition -----

# random competition matrices ----


Omega_kappaN_random <- function(datat,coeft){
  dft<- data.frame("x"=seq(from=1, to=20,by=.1))
  dft$y <- 1/(coeft[1]*(dft$x^(coeft[2])))
  datat<-subset(datat,n==8|n==6|n==10|n==12|n==14|n==16|n==20, select=colnames(datat))
  datat$nf <- factor(datat$n, levels=c(6,8,10,12,14,16,20))
  ggplot(data= datat,aes(x=connr*n, y= meanvol/2)) +
    geom_point(aes(shape=nf,colour = nf),size=3,show.legend = F)+
    #geom_line(linetype="dashed", alpha = 0.1, show.legend = F)+
    labs(x=TeX("\\kappa S"), y=TeX("\\Omega"),colour="S",shape="S")+
    scale_shape_manual(values=shapes)+
    scale_color_manual(values=cols)+
    themet+
    geom_line(data=dft, aes(x=x, y= y/2), color="black")+
    #geom_point(data=data_bk, aes(x=n, y= meanvol), color="black", shape=2, size=3)+
    ylim(0,.25)
}

Omega_per_sp_kappa_random <- function(datat){
  datat<-subset(datat,n==6|n==10|n==20, select=colnames(datat))
  datat$nf <- factor(datat$n, levels=c(6,10,20))
  ggplot(data= datat,aes(x=connr, y= (meanvol/2)^(1/n))) +
    geom_point(aes(shape=nf,colour = nf),size=3,show.legend = F)+
    geom_line(aes(colour = nf),alpha=0.15, linetype=2,show.legend = F)+
    #geom_line(linetype="dashed", alpha = 0.1, show.legend = F)+
    labs(x=TeX("\\kappa"), y=TeX("\\Omega^(1/S)"),colour="S",shape="S")+
    scale_shape_manual(values=shapes[c(1,3,7)])+
    scale_color_manual(values=cols[c(1,3,7)])+
    geom_line(data=datat[datat$P==6|datat$P==10|datat$P==20,],aes(x=connr, y=(meanvol/2)^(1/(n)), colour=nf),alpha=0.15, linetype=2,show.legend = F)+
    themet 
}

Omega_kappa_61020_ran <- function(datat,coef){
  dft <- c()
  for (i in c(6,10,20)){
    dftt<- data.frame("x" = seq(from=0, to=1,by=.01))
    dftt$S <- i
    dftt$y <- 1/(coef[1]*((dftt$x*dftt$S)^(coef[2])))
    dft<-rbind(dft,dftt)
    }
  datat<-subset(datat,n==6|n==10|n==20, select=colnames(datat))
  datat$nf <- factor(datat$n, levels=c(6,10,20))
  fig <-ggplot()+
    geom_point(data=datat,aes(x=connr, y=meanvol/2, shape=nf, colour=nf),size=3,show.legend = F)+
    scale_shape_manual(values=shapes[c(1,3,7)])+
    geom_line(data=dft, aes(x=x, y=y/2, colour=as.character(S),), show.legend = FALSE)+
    #geom_line(data=fitted_dS, aes(x=xS, y=y8))+
    ylim(0,.25)+
    #xlim(min(datat$SC)-.02,20)+
    scale_color_manual(values=cols[c(1,3,7)])+
    labs(y= TeX("\\Omega"), x = TeX(" \\kappa"),  color="S")+
    themet
  return(fig)
}

# different number of species and resources

Omega_m <- function(datat){
  datat<- subset(datat, connr == 0.5| connr == 0.6|connr == 0.75|connr == 0.9 , select = colnames(datat) )
  datat <-subset(datat, m<=9 ,  select = colnames(datat) )
  fig<-ggplot(data=datat, aes(x=m, y=meanvol, colour=as.character(connr)))+
    geom_point(aes(shape=as.character(connr), fill=as.character(connr)),show.legend = F, size=3)+
    geom_line(alpha=.3,linetype=2, show.legend = F)+ 
    scale_shape_manual(values=shapes_K)+
    scale_color_manual(values=colsN)+
    scale_fill_manual(values=colsN)+
    labs(x=TeX("m"), y=TeX("\\Omega"),colour=TeX("\\Kappa"))+
    # scale_color_discrete(breaks = c(5:10,15))+
    themet
  fig
  return(fig)
}

Omega_psp_m <-  function(datat){
  datat<- subset(datat, connr == 0.5| connr == 0.6|connr == 0.75|connr == 0.9, select = colnames(datat) )
  datat <-subset(datat, m<=9 ,  select = colnames(datat) )
  fig<-ggplot(data=datat, aes(x=m, y=meanvol^(1/n), colour=as.character(connr)))+
  geom_point(aes(shape=as.character(connr), fill=as.character(connr)),show.legend = F, size=3)+geom_line(alpha=.3,linetype=2, show.legend = F)+ 
  scale_shape_manual(values=shapes_K)+
    scale_color_manual(values=colsN)+
    scale_fill_manual(values=colsN)+
  labs(x=TeX("m"), y=TeX("\\Omega^(1/n)"),colour=TeX("\\Kappa"))+
  themet
  return (fig)
}

make_leg_m <-  function(datat){
  datat<- subset(datat, connr == 0.5| connr == 0.6|connr == 0.75|connr == 0.9 & m<=9, select = colnames(datat) )
  fig<-ggplot(data=datat, aes(x=m, y=meanvol))+
    geom_point(aes(shape=as.character(connr),colour=as.character(connr), fill=as.character(connr)),show.legend = T, size=3)+ 
    scale_shape_manual(values=shapes_K)+
    scale_color_manual(values=colsN)+
    scale_fill_manual(values=colsN)+
    labs(x=TeX("m"), y=TeX("\\Omega^(1/n)"),colour=TeX("\\kappa"), shape=TeX("\\kappa"),fill=TeX("\\kappa"))+
    themet
  leg <- get_legend(fig) 
  return(as_ggplot(leg))
}


#nestedness
Omega_nest<-function(datat,S)  {
  fig<-ggplot(data=datat)+
    geom_point(aes( y=meanvol,x=nr, colour=as.character(connc), shape=as.character(connc), fill=as.character(connc)), size=3,show.legend = F)+
    geom_line(aes( y=meanvol,x=nr, colour=as.character(connc)), alpha=.3,linetype=2, show.legend = F)+
    scale_color_manual(values=colsN)+
    scale_fill_manual(values=colsN)+
    scale_shape_manual(values=shapes_K)+
    scale_x_continuous(breaks = seq(0,1,by=0.1) )+
    labs(x="NODFc", y=TeX("\\Omega"), title = paste("S=",S),colour = TeX("\\Kappa"))+
    themet
  return(fig)
}
make_leg_nest <- function(datat){
  fig<-ggplot(data=datat)+
    geom_point(aes( y=meanvol,x=nr, colour=as.character(connc), shape=as.character(connc), fill=as.character(connc)), size=3,show.legend = T)+
    scale_color_manual(values=colsN)+
    scale_fill_manual(values=colsN)+
    themet+
    scale_shape_manual(values=shapes_K)+
    guides(shape = guide_legend(nrow = 1))+
    theme(legend.direction="horizontal", legend.position = "top")+
    labs(x="NODFc", y=TeX("\\Omega"),colour = TeX("\\kappa :"), shape= TeX("\\kappa :"), fill= TeX("\\kappa :"))
  leg <- get_legend(fig) 
  return(as_ggplot(leg))
}



plotgamma  <-function(statsi,S){ 
  statsi$gamma <- statsi$diff
  statsi$gamma[statsi$diff == 0] <- 1
  statsi$gamma[statsi$diff == 3] <- round(1/(.6*5+1),2)
  statsi$gamma[statsi$diff == 6] <- round(1/(.6*10+1),2)
  statsi$gamma[statsi$diff == 1] <- round(1/(.6*2.5+1),2)
  statsi$gamma[statsi$diff == 4] <- round(1/(.6*7.5+1),2)
  #statsi$gamma <- factor(statsi$gamma, c(1,round(1/6,2), round(1/11,2)))
  statsi$gamma <- as.character(statsi$gamma)
  p0 <- statsi$conn05[2]
  p1 <- statsi$conn05[floor(length(unique(statsi$conn05))/4)]
  p2 <- statsi$conn05[floor(length(unique(statsi$conn05))/2)]
  p3 <- statsi$conn05[floor(3*length(unique(statsi$conn05))/4)]
  p4 <- statsi$conn05[length(unique(statsi$conn05))-2]
  statsi$connc <- as.character(statsi$conn05)
  plott <- ggplot(data=statsi[ statsi$conn05 == p2| statsi$conn05 == p4 | statsi$conn05 == p0 ,], aes(x=gamma, y=meanvol, group=connc, colour=connc))+
    geom_point(size=3,)+
    geom_line(linetype = "dashed", alpha = .6, show.legend = FALSE)  +
    ylim(0,.25)+
    #scale_color_manual(values=cols)+
    labs(y= TeX("\\Omega"), x = TeX("\\gamma"), title =paste("S=",as.character(S), sep=""), color=TeX("\\kappa"))+
    guides(shape = guide_legend(nrow = 1))+ themeh
  return(plott)
}
