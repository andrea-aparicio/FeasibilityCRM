setwd("/udd/spaap/Projects/Feasibility_package/")
source("code/FeasibilityVolume/functions-MACR.R")
source("code/FiguresVolume/Figure_functions.R")
setwd("data")

#Figure 3 ----

#-Consumer-resource model
statsZs<- read.csv("conn_Z0.001_Th0.8_stats.csv")
statsZs$SC <- statsZs$conn05*statsZs$P 
statsZs$Pf<-factor(statsZs$P, levels=c(6,8,10,12,14,16,20))
statsZs<- make_logs(statsZs)
coefZs <- data_fit(statsZs)
coefZs <- round(coefZs,1)
fittedZs<-make_fitted_data(coefZs,c(6,10,20))

#--three = T plots only 3 community sizes, fit = T plots the poowerlaw curve
Omega_kappa_61020(statsZs,fittedZs,three=T, fit=T) #fig 3.a
Omega_kappaP(statsZs, fittedZs) #fig 3.c 
Omega_per_sp_kappa(statsZs) #fig 3.b
#make_legend_620(statsZs,"V")

#-Random competition 
statsRcomp <- read.csv("stats_randComp_realk.csv")
Omega_kappaN_random(statsRcomp, c(1.5,.9)) #fig 3.d
Omega_per_sp_kappa_random(statsRcomp) #fig 3.e
Omega_kappa_61020_ran(statsRcomp, c(1.5,.9)) #fig 3.d

#Figure 4 ----
#-Different number of species and resources 
#--Linear growth
stats_ndiffmZs <- read.csv("stats_ndiffm_sz.csv")
Omega_m(stats_ndiffmZs) #fig 4.a
Omega_psp_m(stats_ndiffmZs) #fig 4.b
#--Logistic growth
stats_ndiffmZl <- read.csv("stats_ndiffm.csv")
Omega_m(stats_ndiffmZl) #fig 4.c
Omega_psp_m(stats_ndiffmZl) #fig 4.d
make_leg_m(stats_ndiffmZl)

#Figure 5 ----
#-Nestedness 

stats_nest_10 <- read.csv("stats_nest_10.csv") 
stats_nest_15 <- read.csv("stats_nest_15.csv")
stats_nest_20 <- read.csv("stats_nest_20.csv")

Omega_nest(stats_nest_10, 20) #fig 5.b.1
Omega_nest(stats_nest_15, 30) #fig 5.b.2
Omega_nest(stats_nest_20, 40) #fig 5.b.3


#Figure S1 ----
#-Different distributions
#--Tradeoff (Dirichlet distribution) ----

dataAllTrade <- c()
statsTrade <-c()
for (i in c(3,5,10)){
  datat <-  readData("Parameterizations/it_kappa_trad",i,TRUE,1,25)
  datat$S <- i*2
  datat$vol <- datat$vol/2
  statst <- getStats_b(datat)
  statst$S<-as.character(i*2)
  dataAllTrade <- rbind(dataAllTrade,datat)
  statsTrade <- rbind(statsTrade, statst)
}
statsTrade$SC <- as.numeric(statsTrade$S)*statsTrade$conn05

Omega_kappa_Diri <- ggplot()+
  geom_point(data=statsTrade,aes(x=conn05, y=meanvol, shape=S, colour=S),size=3,show.legend = T)+
  geom_line(data=statsTrade,aes(x=conn05, y=meanvol, colour=S), linetype = "dashed", alpha=.5, show.legend = T)+
  scale_shape_manual(values=shapes[c(3,7,1)])+
  ylim(0,.25)+
  scale_color_manual(values=cols[c(3,7,1)])+
  labs(y= TeX("\\Omega"), x = TeX(" \\kappa"),  color="S", shape="S", title = "Dirichlet distribution")+
  themet
Omega_kappaS_Diri <- ggplot()+
  geom_point(data=statsTrade,aes(x=SC, y=meanvol, shape=S, colour=S),size=3,show.legend = F)+
  geom_line(data=statsTrade,aes(x=SC, y=meanvol, colour=S), linetype = "dashed", alpha=.5, show.legend = F)+
  scale_shape_manual(values=shapes[c(3,7,1)])+
  ylim(0,.25)+
  #xlim(min(statsAllSC$SC)-.02,20)+
  scale_color_manual(values=cols[c(3,7,1)])+
  labs(y=TeX("\\Omega"), x = TeX(" \\kappa S"),  color="S")+
  themet

Omega_kappa_Diri #fig S1(1,2)
Omega_kappaS_Diri #fig S1(2,2)

#--Uniform distribution ----

dataAllUnif <- c()
statsUnif <-c()
for (i in c(3,5,10)){
  datat <-  readData("Parameterizations/it_kappa_unif",i,TRUE,1,25)
  datat$S <- i*2
  datat$vol <- datat$vol/2
  statst <- getStats_b(datat)
  statst$S<-as.character(i*2)
  dataAllUnif <- rbind(dataAllUnif,datat)
  statsUnif <- rbind(statsUnif, statst)
}
statsUnif$SC <-  as.numeric(statsUnif$S)*statsUnif$conn05

Omega_kappa_SC_Unif <- ggplot()+
  geom_point(data=statsUnif,aes(x=SC, y=meanvol, shape=S, colour=S),size=3,show.legend = F)+
  geom_line(data=statsUnif,aes(x=SC, y=meanvol, colour=S), linetype = "dashed", alpha=.5, show.legend = F)+
  scale_shape_manual(values=shapes[c(1,3,7)])+
  ylim(0,.25)+
  scale_color_manual(values=cols[c(1,3,7)])+
  labs(y= TeX("\\Omega"), x = TeX(" \\kappa S"),  color="S")+
  themet

Omega_kappa_Unif <- ggplot()+
  geom_point(data=statsUnif,aes(x=conn05, y=meanvol, shape=S, colour=S),size=3,show.legend = F)+
  geom_line(data=statsUnif,aes(x=conn05, y=meanvol, colour=S), linetype = "dashed", alpha=.5, show.legend = F)+
  scale_shape_manual(values=shapes[c(1,3,7)])+
  ylim(0,.25)+
  scale_color_manual(values=cols[c(1,3,7)])+
  labs(y= TeX("\\Omega"), x = TeX(" \\kappa"),  color="S", title="Uniform distribution")+
  themet

Omega_kappa_Unif #fig S1(1,1)
Omega_kappa_SC_Unif #fig S1(2,1)








#Figure S3 ----
sg<- read_makestats_gamma()
stats_g3 <- sg[1]
stats_g5 <- sg[2]
stats_g10 <- sg[3]
plotgamma(stats_g3[[1]],6) #figure S3.b1
plotgamma(stats_g5[[1]], 10) #figure S3.b2
plotgamma(stats_g10[[1]], 20) #figure S3.b3

#Figure S4 ----
Omega_kappa_61020(statsZs,fittedZs,three=F, fit=F)
#Figure S5 ----
statsZl<- read.csv("conn_Z1_Th0.8_stats.csv")
statsZl$conn05 <- statsZl$connr
statsZl$SC <- statsZl$conn05*statsZl$P 
statsZl$Pf<-factor(statsZl$P, levels=c(6,8,10,12,14,16,20))
statsZl<- make_logs(statsZl)
coefZl <- data_fit(statsZl)
coefZl <-round(coefZl,1)
coefZl[1] <- coefZl[1]
coefZl[2] <- coefZl[2]
fittedZl<-make_fitted_data(coefZl,c(6,10,20))
colnames(fittedZl)
fittedZl$y<-fittedZl$y/2
Omega_kappa_61020(statsZl,fittedZl,T,F) #Figure S5.a
Omega_per_sp_kappa(statsZl) #Figure S5.b
#Figure S7 ----
logloga_all(statsZs)
#Figure S8 ----
#Omega_kappa_all(statsZs) 
fittedZsall<-make_fitted_data(coefZs,2*c(3:8,10))
Omega_kappa_all_f(statsZs, fittedZsall)
#Figure S9 ----
statsZlK <- subset(statsZl, conn05>0.4, select = colnames(statsZl))
coefZlK <- data_fit(statsZlK)
coefZlK <-round(coefZlK,1)
fittedZlKall<-make_fitted_data(coefZlK,c(6,8,10,12,14,16,20))
fittedZlKall$y <- fittedZlKall$y/2
Omega_kappaP(statsZl, fittedZl)  #Figure S9.a
Omega_kappaP(statsZlK, fittedZlKall) #Figure S9.b













