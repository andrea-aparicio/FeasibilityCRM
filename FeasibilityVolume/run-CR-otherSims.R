setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("functions-MACR.R")

#interaction strength ~abs(N(0,stren))
stren <- .6
#efficiency ~  ~abs(N(0,ep_s))
ep_s <- .8
#zeta ~  ~abs(N(0,prop))
prop <-.001

#number of realizations
times <- 50
#number of repetitions
reps <- 20

#consumer-resource model with varying connectance and selecting distributions ----
#-select distribution ("uniform", "normal", "tradeoff")
distr <- "uniform"

for (k in c(3,5,10)){
  n<-k
  m<-k
  print("n")
  print(k)
  print("start")
  for (j in 1:reps){
    print(j)
    it_kappat <- iterate_log_dist(stren, times, ep_s,n,m,prop,distr)
    write.csv(it_kappat, paste("it_kappa_", distr, as.character(k), "_",as.character(j),".csv", sep=""))
  }
}
 
#consumer-resource model with varying connectance and nestedness ----

for (k in c(3,5,10)){
  n<-k
  m<-k
  print(k)
  print("count")
  for (j in 1:reps){
    print(j)
    it_kappan <- iterate_nest(stren, times, ep_s,n,m,prop)
    print("end")
    write.csv(it_kappan, paste("it_kappa_nest2", as.character(k), "_",as.character(j),".csv", sep=""))
  }
}


#consumer-resource model with varying connectance and magnitudes of Z ----

n=5
m=n

for (j in 1:reps){
  print(j)
  it_Zeta <-iterate_log_Z(stren, times, ep_s,n,m)
  write.csv(it_Zeta, paste("it_Zeta", as.character(n), "_",as.character(j),".csv", sep=""))
}

