setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
args = commandArgs(trailingOnly=T)
source("functions-MACR.R")

#read from console which kind of simulation is running (uncomment below line)
mode = args[1]
# or define if running directly ("MCRM_connectance", "random_competition, "n<m")
#mode = 
print(mode)
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

# consumer-resource model with varying connectance ----
if (mode == "MCRM_connectance") {
  #vector of number of consumers 
  kv <- c(3:8,10)
  for (k in kv){
    n<-k
    #n = m
    m<-k 
    print("n")
    print(n)
    print("start")
    for (j in 1:reps){
      name <- paste("/udd/spaap/Projects/Feasibility/data/MCRM/connectance/conn_Z",as.character(prop),"_Th",as.character(ep_s),"_",as.character(k), "_",as.character(j),".csv", sep="")
      it_kappat <- iterate_log(stren, times, ep_s,n,m,prop)
      write.csv(it_kappat, name)
      print(name)
    }
  }
} 


# random competition  ----

if (mode == "random_competition"){
  
  paramt = args[2]
  
  #if varying connectance
  if (paramt == "kappa"){
    connv <- seq(0.1,1,b=0.1)
  } else{
    connv <- 0.5
  }
  
  #if varying the value of main diagonal
  if (paramt == "diag"){
    diagv <- c(.1,seq(.2,1,by=.2),2:10)
  } else {
    diagv <- 1
  }
  
  #if varying the competition strength
  if (paramt=="stren"){
    strenv <- c(.1,seq(.2,1,by=.2),2:10)
  } else {
    strenv <- 1
  }
  
  #vector community sizes
  nv <- seq(2,20, by=1)
  
  for (j in 1:reps){
    name <- paste("/udd/spaap/Projects/Feasibility/data/random_compe/random_",
                  paramt,as.character(j),".csv", sep="")
    print(name)
    itt<-iterate_random_comp(strenv, times,nv, connv, diagv)
    write.csv(itt, name)  
  }
  
}

# number of species (n) <  number of resources (m) ----

if (mode == "n<m"){
  #set number of resources to 5
  n=5
  for (m in 5:10){
    print(m)
    for (j in 1:reps){
      print(j)
      itt<-iterate_log_ndm(stren, times, ep_s,n,m,prop)
    }
  }
}




