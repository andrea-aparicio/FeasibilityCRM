pathgen <- "/udd/spaap/Projects/Feasibility_package/"
setwd(pathgen)

source("code/FeasibilityVolume/functions-MACR.R")

make_stats <- function(name, nv){
  stats <- c()
  for (n in nv){
    datat <- readData(name,n,F,1,20)
    datat <- na.omit(datat)
    statst<-getStats_b(datat)
    statst$P<-2*n
    stats <- rbind( stats,statst)
  }
  stats<-na.omit(stats)
  write.csv(stats,paste(pathgen,name,"stats.csv",sep=""))
  print(paste(pathgen,name,"stats.csv",sep=""))
}


make_stats_dir <- function(patht, namet, nv){
  datal <- c()
  for (i in nv){
    fname <- paste(namet,as.character(i),sep="")
    pathn <- paste(pathgen,patht,sep="")
    fnames <- list.files(path = pathn, pattern = fname)
    for (n in fnames){
      datat <- read.csv(paste(pathn,n,sep=""))
      datat$P <- i*2
      datal <- rbind(datal, datat)
      print(n)
    }
    datal <- datal[!is.na(datal$vol),]
  }
  datal$connr <- round(2*datal$conn,1)/2
  stats <- ddply(datal, .(P,connr), summarize, meanvol=mean(vol), servol=sd(vol)/sqrt(length(vol)))
  write.csv(stats,paste(pathn,namet,"stats.csv",sep=""))
  print(paste(pathn,namet,"stats.csv",sep=""))
}

make_stats_dir("data/connectance/", "conn_Z1_Th0.8_", c(3:8,10)) #read and summarize data for Fig S5

make_stats("data/connectance/conn_Z0.001_Th0.8_", c(3:8,10)) #read and summarize data for Fig 3

make_stats("data/connectance/conn_Z0_Th0.8_", c(3,5,10))


#make stats random competition
datat <- c()
for (i in 1:6){
  datat <- rbind(datat, read.csv(paste("data/random_compe/random_kappa",i,".csv", sep="")))
}
datat$connr <- round(2*datat$conn_calc,1)/2
stats=ddply(datat, .(n,connr), summarize, meanvol=mean(vol), servol=sd(vol)/sqrt(length(vol)))
write.csv(stats,"data/random_compe/random_kappa_stats.csv")


