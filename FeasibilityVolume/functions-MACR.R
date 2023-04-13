# library(matlib)
require(Matrix)
library(ggplot2)
library(tidyr)
library(plyr)
library(maxnodf)
library(ggpubr)
library(mixtools)
library(mvtnorm)
library(sna)
library(igraph)



#Make matrices ----

#-remove elements of consumption matrix
#--with probability
RemoveLinksComp <- function(prob, vecC,n){
  A <- matrix(vecC, nrow=n, ncol=n)
  zeroesi <- c()
  zeroesj <- c()
  ms <- 1
  for (i in 1:(n)){
    for (j in ms:(n)){
      if (runif(1) <= (1-prob)){
        zeroesi <- rbind(zeroesi, i)
        zeroesj <- rbind(zeroesj, j)
      }
    }
    ms = ms+1
  }
  zeroes = cbind(rbind(zeroesi,zeroesj),rbind(zeroesj,zeroesi))
  if (!is.null(dim(zeroes))){
    for (i in 1:dim(zeroes)[1]){
      A[zeroes[i,1], zeroes[i,2]]<-0
    }
  }
  return(A)
}
#--to match exact connectance
#---for competition matrices
RemoveLinksCompKappa <- function(conn, vecC,n){
  avail <- c()
  ms <- 2
  for (i in 1:(n-1)){
    for (j in ms:n){
      avail <- rbind(avail,c(i,j))
    }
    ms = ms+1
  }
  avail
  nzeros <- round((1-conn)*dim(avail)[1])
  connected <- F
  if (dim(avail)[1]>1){
    while(connected == F){
      zeros_pos <- sample(1:dim(avail)[1], nzeros)
      n_zeros_pos <- setdiff(1:dim(avail)[1], zeros_pos)
      if (length(zeros_pos)!=0){
        zeros<- avail[zeros_pos,]
      } else {
        zeros=NULL
      }
      notzeros <- avail[n_zeros_pos,]
      notzeros <- rbind(notzeros,cbind(notzeros[,2],notzeros[,1]))
      graph <- graph_(notzeros, from_edgelist(directed=F))
      connected <- is_connected(graph)
    }
  } else {
    zeros=NULL
  }
  A <- matrix(vecC, nrow=n, ncol=n)
  if (!is.null(dim(zeros)[1])){
    for (i in 1:dim(zeros)[1]){
      A[zeros[i,1],zeros[i,2]]<-0
      A[zeros[i,2],zeros[i,1]]<-0
    }
}
  return(A)
}
#---for consumption matrices
RemoveLinksC <- function(conn, vecC, n, m){
  # if (conn<=n/(n*m)){
  #    Cn <- diag(diag(matrix(vecCt, ncol = n, nrow = m)))
  #  } else{
  #     zeroes <- sample(c(rep.int(1,floor(m*n*conn)), rep.i                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  nt(0,(m*n-floor(m*n*conn)))))
  #     vecCt[which(zeroes==0)] <- 0
  #     Cn <- matrix(vecCt, ncol = n, nrow = m)  
  #   }
  # }
  vecCt <- vecC
  zeroes <- sample(c(rep.int(1,floor(m*n*conn)), rep.int(0,(m*n-floor(m*n*conn)))))
  vecCt[which(zeroes==0)] <- 0
  Cn <- matrix(vecCt, ncol = n, nrow = m)  
  # }
  #   }
  # }
  return(Cn)
}

#-make efficiency vector (theta)
MakeT <- function(ep, nt,n){
  et <- abs(rnorm(nt, mean = 0, sd = ep))
  th <- diag(et,n)
}
#--make theta correlated with the competition matrix
MakeT_dep <- function(ep, nt,n, fact){
  et <- c()
  for (i in 1:nt){
    et <- rbind(et,abs(rnorm(1, mean = 0, sd = ep*fact[i])))
    }
  th <- diag(as.vector(et),n)
}
#-scramble edges to vary nestedness preserving connectance
adjust_nest_conn <- function(c, nestdes, difdes){
  m<-nrow(c)
  n<-ncol(c)
  nmat2 <- c
  dif <- abs(get_nest(nmat2)-nestdes)
  count <- 0
  temp <- 0.5
  while (dif>difdes & count<(n*m*m)){
    # print(count)
    nz=which(nmat2!=0,arr.ind = T)
    yz= which(nmat2==0, arr.ind = T)
    #print(nz)
    #print(yz)
    if (length(nz)>2 & length(yz)>2){
      nz = nz[nz[,1]!= nz[,2],] #don't move the main diagonal
      yz = yz[yz[,1]!= yz[,2],] #don't move the main diagonal
    }
    nmat <- nmat2
    if (length(nz)>2 & length(yz)>2){
      indmove <- nz[sample(1:length(nz[,1]),1, replace = F),]
      inddest <- yz[sample(1:length(yz[,1]),1, replace = F),]
      nmat[inddest[1], inddest[2]] <- nmat[indmove[1],indmove[2]]
      nmat[indmove[1],indmove[2]] <- 0
    }
    if (rankMatrix(nmat)==n){
      count<-count+1
      difn <- abs(get_nest(nmat)-nestdes)
      if (difn<dif | runif(1,0,1)<=temp){
        count<-0
        temp <- temp/2
        # print("ch")
        dif<-difn
        nmat2<-nmat
      }
    }
  }
  return(list("mat"=nmat2, "nest"=get_nest(nmat2)))
}
#--decreasing only
decrease_nest_conn <- function(c){
  m<-nrow(c)
  n<-ncol(c)
  nmat2 <- c
  nestc <- get_nest(c)
  nestnmat <- nestc
  flag <- 0
  while (nestnmat>= nestc & flag < n*m*n){
    #non zero elements
    nz=which(nmat2!=0,arr.ind = T)
    #zero elements
    yz= which(nmat2==0, arr.ind = T)
      # if (length(nz)>2 & length(yz)>2){
      #   nz = nz[nz[,1]!= nz[,2],] #don't move the main diagonal
      #   yz = yz[yz[,1]!= yz[,2],] #don't move the main diagonal
      # }
    #nmat <- nmat2
    indmove <- nz[sample(1:length(nz[,1]),1, replace = F),]
    inddest <- yz[sample(1:length(yz[,1]),1, replace = F),]
    nmat2[inddest[1], inddest[2]] <- nmat2[indmove[1],indmove[2]]
    nmat2[indmove[1],indmove[2]] <- 0
    if (rankMatrix(nmat2)[[1]]==n){
      # print("fr")
      #count<-count+1
      nestnmat<- get_nest(nmat2)
      # print(nestnmat)
      # print(nestc)
      flag <- flag+1
      #nmat2<-nmat
    }
  }
  return(list("mat"=nmat2, "nest"=get_nest(nmat2), "flag"=flag))
}



#Volume calculation ----
#-find generating vectors
#--function retrieved from https://github.com/clsong/Ecology_Song-Saavedra_2018 and adapted
spanned_vectors <- function(A){
  G <- matrix(0, ncol=ncol(A), nrow=nrow(A))
  for(k in 1:nrow(A)) G[,k] <- A[,k]/sqrt(sum(A[,k]^2))
  G
}
#-calculate feasibility volume
#--function retrieved from https://github.com/clsong/Ecology_Song-Saavedra_2018
Omega <- function(alpha) {
  S <- nrow(alpha)
  omega <- function(S, Sigma) {
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- try(pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma))
    if (class(d)=="numeric"){
      out <- d[1]^(1 / S)
    } else {
      out <- NaN
    }
    return(out)
  }
  f <- function(m) class(try(solve(t(m) %*% m), silent = T)) == "matrix"
  
  if (f(alpha)[1] == FALSE) {
    return(0)
  } else {
    Sigma <- solve(t(alpha) %*% alpha)
    return(omega(S, Sigma))
  }

}
#-build matrix A, calculate volume, connectance and nestedness
get_vol_log <- function(C2, Z,th,n,m,nest){
  A <- rbind(cbind(-Z,-C2),cbind(th %*% t(C2), matrix(0, nrow=n, ncol=n)))
  Vt<- Omega(spanned_vectors(A))
  real_conn<-  sum(C2!=0)/(n*m)
  if (nest == T){
    nest<-get_nest(C2)
  }else{nest<-0}
  return(list("vol"=Vt, "con"=real_conn, "nest"=nest))
}
#-calculate volume and connectance (used for competition matrices)
get_vol_comp <- function(A,n){
  Vt<- Omega(spanned_vectors(A))
  real_conn <-  sum(A!=0)/(n*n)
  return(list("vol"=Vt, "con"=real_conn))
}
#-calculate nestedness 
#--manually
get_nest_o <- function(B){
  B[B>0]<-1
  n <- nrow(B)
  m <- ncol(B)
  ki <- rep(0,n)
  ka <- rep(0,m)
  P <- matrix(0,nrow=n,ncol=n)
  Q <- matrix(0,nrow=m,ncol=m)
  Pt <- matrix(0,nrow=n,ncol=n)
  Qt <- matrix(0,nrow=m,ncol=m)
  for (i in 1:n){
    ki[i] <- sum(B[i,]>0)
    for (j in 1:n){
      temp<-0
      for (a in 1:m){
        temp=temp+(B[i,a]*B[j,a])
      }
      P[i,j] = temp
    }
  }
  for (a in 1:m){
    ka[a] <- sum(B[,a]>0)
    for (b in 1:m){
      temp<-0
      for (j in 1:n){
        temp = temp + (B[j,a]*B[j,b])
      }
      Q[a,b] = temp
    }
  }
  for (i in 1:n){
    for (j in 1:n){
      if (i!=j){
        Pt[i,j] <- P[i,j]/min(ki[i], ki[j])
      }
    }
  }
  for (a in 1:m){
    for (b in 1:m){
      if (a!=b){
        Qt[a,b] <- Q[a,b]/min(ka[a], ka[b])
      }
    }
  }
  NODF <- (sum(Qt[upper.tri(Qt)] + Pt[upper.tri(Pt)]))/(((n*(n-1)/2))+((m*(m-1))/2))
  return(NODF)
}
#--using maxnodf package
get_nest <- function(B){
  B[B!=0]<-1
  NODF <- nodf_cpp(B)
  return(NODF)
}

# Data generation: generate communities and calculate volume and other properties  ----

#-consumer-resource model with varying connectance
#--to vary niche overlap, follow instructions for diagm below
iterate_log <- function(stren, times, ep_s,n,m,prop){
  #number of nonzero elements in C 
  seqs = append(0,seq(5,n*m,ceiling(n*m/10)))
  #number of iterations
  runs=seq(1,times,1)
  #result vector
  res <- c()
  for(run in runs) {
    #consumption rates
    vecC <-  abs(rnorm(n*m, mean = 0, sd = stren))
    #zeta
    Z <- prop*diag(abs(rnorm(m, mean = 0, sd = stren))) 
    #efficiancy rate
    th <- MakeT(ep_s,n,n)
    for (connv in seqs) {
      #calculate connectance
      conn = connv/(n*m)
      #diagm is the mean for the species favorite resource consumption rate. 
      # keep at 0 for maximum niche overlap change to diag>0 to decrease the niche overlap
      for(diagm in 0){
        Vt = NaN
        while(is.nan(Vt)==TRUE||Vt==0){
          C <- RemoveLinksC(conn, vecC, n, m)
          #comment the following line for max niche overlap, and uncomment (***)
          C2<-C+diag(abs(rnorm(n, mean = diagm, sd = stren)))
          C2<-C
          #add nonzero diagonal elements so that every species has a favorite resource
          #diag(C2)<-abs(rnorm(n, mean = diagm, sd = stren)) (***)
          #auxiliary to calculate difference in favorite and non favorite consumption rate
          Ct <- C2+diag(rep(100,n))
          meanO <- mean(Ct[c(Ct!=0 & Ct<100) ])
          meanD <- mean(diag(C2))
          #assemble A matrix, calculate volume, connectance and nestedness
          temp<- get_vol_log(C2,Z,th,n,m, F)
          Vt=temp$vol
          rest <-  c(stren, Vt,meanO, meanD,temp$con, diagm,temp$nest,10)
        }
        res <- rbind(res,rest)
      }
    }
  }
  resdf = data.frame(res)
  colnames(resdf)<-c("stren","vol", "meanO", "meanD", "conn","diag","nest","desnest")
  return(resdf)
}

#-random competition matrices selecting distribution
iterate_random_comp  <- function(strenv, times,nv, connv, diagv){
  runs=seq(1,times,1)
  res <- c()
  for(run in runs) {
    print(run)
    for (n in nv){
      connv<- connv[connv>(n/(n*n))]
      for (stren in strenv){
        for (conn in connv){
          #print(conn)
          for (diagt in diagv){
            vecC <-  -abs(rnorm(n*n, mean = 0, sd = stren))
            Vt = NaN
            while(is.nan(Vt)==TRUE||Vt==0){
              A <- RemoveLinksCompKappa(conn, vecC, n)
              A=A-diag(diag(A))-diag(rep(diagt,n))
              temp<- get_vol_comp(A,n)
              Vt=temp$vol
              rest <-  c(n, stren, conn, diagt, Vt,temp$con)
            }
            res <- rbind(res,rest)
          }
        }
      }
    }
  }
  resdf = data.frame(res)
  colnames(resdf)<-c("n","stren", "conn", "diagt", "vol","conn_calc")
  return(resdf)
}

#-number of species (n) <  number of resources (m)
iterate_log_ndm<- function(stren, times, ep_s,n,m,prop){
  seqs = seq(m,n*m,ceiling(n*m/10))
  runs=seq(1,times,1)
  res <- c()
  for(run in runs) {
    for (connv in seqs){
      vecC <-  abs(rnorm(n*m, mean = 0, sd = stren))
      Z <- prop*diag(abs(rnorm(m, mean = 0, sd = stren)))       
      th <- MakeT(ep_s,n,n)
      conn = connv/(n*m)
      Vt = NaN
      while(is.nan(Vt)==TRUE||Vt==0){
        zero_el <- n
        while (zero_el!=0){
          C <- RemoveLinksC(conn, vecC, n, m)
          zero_el <- sum(colSums(C)==0,rowSums(C)==0)
        }
        C2<-C
        diag(C2)<-abs(rnorm(n, mean = diagm, sd = stren))
        
        temp<- get_vol_log(C2,Z,th,n,m, F)
        temp
        Vt=temp$vol
        rest <-  c(stren, Vt,temp$con)
      }
      res <- rbind(res,rest)
    }
  }
  resdf = data.frame(res)
  colnames(resdf)<-c("stren","vol", "conn")
  return(resdf)
}

#-consumer-resource model with varying connectance, selecting distribution
iterate_log_dist <- function(stren, times, ep_s,n,m,prop, distr){
  seqs = append(0,seq(5,n*m,ceiling(n*m/10)))
  runs=seq(1,times,1)
  res <- c()
  # res=matrix(data=0, nrow=length(runs)*length(seqs), ncol=2)
  # cou = 1
  for(run in runs) {
    #print("run")
    #print(run)
    if (distr == "normal"){
      vecC <-  abs(rnorm(n*m, mean = 0, sd = stren))
    } else if (distr == "uniform"){
      vecC <-  runif(n*m, min = 0.01, max = 2*stren)
    } else if (distr == "tradeoff"){
      vecC <- rep(1,n*m)
    }
    Z <- diag(abs(rnorm(n, mean = 0, sd = prop)))        
    th <- MakeT(ep_s,n,n)
    for (connv in seqs) {
      #print("connV")
      #print(connv)
      conn = connv/(n*m)
      for(diagm in 0){
        #for(diagm in c(stren*2.5, stren*7.5)){
        #for(diagm in c(stren*10)){ #!!!!!!!!!!!!!!!!!!!!!!!!!
        Vt = NaN
        while(is.nan(Vt)==TRUE||Vt==0){
          C <- RemoveLinksC(conn, vecC, n, m)
          
          if (distr == "tradeoff"){
            rows_nz <- rowSums(C)
            for (rowt in 1:n){
              vals_n <- rows_nz[rowt]
              #generate as many new values as non-zero elements in the row of C
              newr <- 2*stren*rdiric(n=1,shape= runif(vals_n,min=1,max=10))
              print(newr)
              coltt <- 1
              for (colt in 1:n){
                if (C[rowt,colt]!=0){
                  C[rowt,colt]<-newr[coltt]
                  coltt<-coltt+1
                }
              }             
            }
          }
          
          #print("d")
          #print(diagm)
          #C2<-C+diag(abs(rnorm(n, mean = diagm, sd = stren)))
          C2<-C
          diag(C2)<-abs(rnorm(n, mean = diagm, sd = stren))
          #print(C2)
          Ct <- C2+diag(rep(100,n))
          meanO <- mean(Ct[c(Ct!=0 & Ct<100) ])
          meanD <- mean(diag(C2))
          # print(meanD)
          temp<- get_vol_log(C2,Z,th,n,m)
          # A <- rbind(cbind(-Z,-C2),cbind(th %*% t(C2), matrix(0, nrow=n, ncol=n)))
          # Vt<- Omega(spanned_vectors(A))/2
          Vt=temp$vol
          # real_conn<-  sum(C2!=0)/(n*m)
          rest <-  c(stren, Vt,meanO, meanD,temp$con, diagm,temp$nest,10)
        }
        
        res <- rbind(res,rest)
        if (conn > 0 & conn < .9) {
          for (desnest in seq(.2,1,.2)){
            #print(desnest)
            #print(C2)
            Cn <- adjust_nest_conn(C2, desnest, .05)
            temp<- get_vol_log(Cn$mat,Z,th,n,m)
            rest <-  c(stren, temp$vol,meanO, meanD,temp$con, diagm,temp$nest,desnest)
            # cou=cou+1
            res <- rbind(res,rest)
          }
        }
      }
    }
  }
  resdf = data.frame(res)
  colnames(resdf)<-c("stren","vol", "meanO", "meanD", "conn","diag","nest","desnest")
  return(resdf)
}

#-consumer-resource model with varying connectance and nestedness
iterate_nest <- function(stren, times, ep_s,n,m,prop){
  seqs = c(0.4, .5, 0.7, 0.9)
  maxNestl <-  vector(mode='list', length=length(seqs))
  names(maxNestl)<-as.character(seqs)
  for (conn in seqs){
    rcbt <- n-1
    while (rcbt<n){
    Cbt <- RemoveLinksC(conn, rep(1,n*m),n,m)
    rcbt <- rankMatrix(Cbt)
    }
    Cbtn <- maxnodf(Cbt,2)
    maxNestl[[as.character(conn)]]<-Cbtn
  }
  print("list built")
  runs=seq(1,times,1)
  res <- c()
  for(run in runs) {
    print(run)
    vecC <-  abs(rnorm(n*m, mean = 0, sd = stren))
    matC <- matrix(vecC, nrow=n, ncol=n)
    Z <- diag(abs(rnorm(n, mean = 0, sd = prop)))        
    th <- MakeT(ep_s,n,n)
    for (conn in seqs) { #change connectance
      print(c("conn",conn))
      nestMax <- maxNestl[[as.character(conn)]]$max_nodf
      Vt = NaN
      while(is.nan(Vt)==TRUE||Vt==0){
        # rankt <- n-1
        # while (rankt<n){
        #   Ctt <- RemoveLinksC(conn, vecC, n, m)
        #   rankt <- rankMatrix(Ctt)[[1]]
        # }
        # nestMax <- maxnodf(Ctt,1)
        Ct <- maxNestl[[as.character(conn)]]$max_nodf_mtx * matC
        flag <- 0
        while (rankMatrix(Ct)[[1]]<n & flag < n*m*m){
          Ctd<-decrease_nest_conn(Ct)
          Ct<-Ctd$mat
          flag <- Ctd$flag
        } 
        temp<- get_vol_log(Ct,Z,th,n,m) ####
        Vt=temp$vol
        nestCt <- temp$nest
        rest <-  c(stren,  Vt,     temp$con, nestCt, nestMax, get_nest_o(Ct))
        #       c("stren","vol",  "conn",     "nest",     "maxnodf",        "nest_o")
      }
      res <- rbind(res,rest)
      while (nestCt > 0.2 & flag < n*m*m){
        Ctd<-decrease_nest_conn(Ct)
        Ct<-Ctd$mat
        temp<- get_vol_log(Ct,Z,th,n,m) ####
        Vt=temp$vol
        nestCt <- temp$nest
        rest <-  c(stren,  Vt,     temp$con, nestCt, nestMax, get_nest_o(Ct))
        res <- rbind(res,rest)
        flag <- Ctd$flag
      }
    }
  }
  resdf = data.frame(res)
  colnames(resdf)<-c("stren","vol",  "conn","nest","maxnodf","nest_o")
  return(resdf)
}

#-consumer-resource model with varying connectance and different magnitude of Z
iterate_log_Z <- function(stren, times, ep_s,n,m){
  seqs = append(0,seq(5,n*m,ceiling(n*m/10)))
  runs=seq(1,times,1)
  res <- c()
  # res=matrix(data=0, nrow=length(runs)*length(seqs), ncol=2)
  # cou = 1
  for(run in runs) {
    vecC <-  abs(rnorm(n*m, mean = 0, sd = stren))
    th <- MakeT(ep_s,n,n)
    for (connv in seqs) {
      # print(connv)
      conn = connv/(n*m)
      Vt = NaN
      while(is.nan(Vt)==TRUE||Vt==0){
        C <- RemoveLinksC(conn, vecC, n, m)
        C<-C+diag(abs(rnorm(n, mean = 0, sd = stren)))
        prop <- .001
        Z <- diag(abs(rnorm(n, mean = 0, sd = prop)))        
        # print(diag(Z))
        A <- rbind(cbind(-Z,-C),cbind(th %*% t(C), matrix(0, nrow=n, ncol=n)))
        Vt<- Omega(spanned_vectors(A))/2
        real_conn<-  sum(C!=0)/(n*m)
        rest <-  c(stren, Vt,real_conn, prop)
      }
      res<-rbind(res,rest)
      for (prop in c(.1, 1)){
        Vt = NaN
        ct=0
        while(is.nan(Vt)==TRUE||Vt==0){
          # print(prop)
          Z <- diag(abs(rnorm(n, mean = 0, sd = prop)))        
          # print(diag(Z))
          A <- rbind(cbind(-Z,-C),cbind(th %*% t(C), matrix(0, nrow=n, ncol=n)))
          Vt<- Omega(spanned_vectors(A))/2
          real_conn<-  sum(C!=0)/(n*m)
          ct=ct+1
          if (ct==500){Vt<-100}
          rest <-  c(stren, Vt,real_conn, prop)
        }
        res<-rbind(res,rest)
      }
    }
    # res[cou,]<- c(real_conn,Vt)
    # cou=cou+1
  }
  resdf = data.frame(res)
  colnames(resdf)<-c("stren","vol", "conn","prop")
  return(resdf)
}


# Read data from disk ----
readData <- function(name,n,mean,ini,fin){
  datat <- c()
  if (length(ini)==1){
    ran=ini:fin
  }
  else { ran=ini}
  for (i in ran){
    datat <- rbind(datat, read.csv(paste(name,as.character(n),"_",as.character(i),".csv", sep="")))
  }
  if (mean==TRUE){ datat$meanO[is.na(datat$meanO)==TRUE]=0 }
  return(datat)
}

readDataStren <- function(name,mean,ini,fin){
  datat <- c()
  if (length(ini)==1){
    ran=ini:fin
  }
  else { ran=ini}
  for (i in ran){
    datat <- rbind(datat, read.csv(paste(name,"_",as.character(i),".csv", sep="")))
  }
  if (mean==TRUE){ datat$meanO[is.na(datat$meanO)==TRUE]=0 }
  return(datat)
}


# Obtain statistics for figures ----
getStatsDiff <- function(datai){
  statsAll <- c()
  #define which niche overlap to use, make diagv=0 for maximal
  diagv <- c(0,1.5,3,4.5,6)
  
  for (i in diagv){
    datat = datai[datai$diag==i,]
    datat$conn05 <- round(2*datat$conn,1)/2
    stats=ddply(datat, .(conn05), summarize,  
                meanvol=mean(vol), servol=sd(vol)/sqrt(length(vol)))
    stats$diff <- as.character(round(mean(datat$meanD-datat$meanO)))
    statsAll <- rbind(statsAll,stats)
  }
  return(statsAll)
}

getStats_b <- function(datat){
  statsAll <- c()
  datat <- na.omit(datat)
  datat$conn05 <- round(2*datat$conn,1)/2
  stats=ddply(datat, .(conn05), summarize,
              meanvol=mean(vol), servol=sd(vol)/sqrt(length(vol)))
    statsAll <- rbind(statsAll,stats)
  return(statsAll)
}

getStats_FixedP <- function(datat){
  statsAll <- c()
  datat <- na.omit(datat)
  datat$conn05 <- round(2*datat$conn,1)/2
  stats=ddply(datat, .(conn05), summarize,
              meanvol=mean(vol), servol=sd(vol)/sqrt(length(vol)))
  statsAll <- rbind(statsAll,stats)
  return(statsAll)
}


treatKappa <- function(data,S){
  data$SC <- S*data$conn
  data$S2C <- (S/2)^2*data$conn
  data$sqSC <- sqrt(S*data$conn)
  data$sqCoC <- sqrt(S/data$conn)
  data$norm <- (data$conn - min(data$conn))/max(data$conn - min(data$conn))
  data$S <- S
  return(data)
}


# plot ----
plotDiff <-function(statsi,j){ 
  plott <- ggplot(data=statsi[statsi$diff == j,], aes(x=conn05, y=meanvol))+
    geom_point(size=2)+
    geom_line()+
    ylim(0,.5)+
    scale_color_manual(values=cols)+
    labs(y= "volume", x = TeX("connectance"), title ="")+
    theme(title=element_text(size=15,face="plain"),
          axis.title=element_text(size=15,face = "plain"), 
          axis.text = element_text(size = 11.5),
          legend.text = element_text(size=11.5))
  return(plott)
}

plotDiffAll <-function(statsi,S){ 
  statsi$gamma <- statsi$diff
  statsi$gamma[statsi$gamma == 0] <- 1
  statsi$gamma[statsi$gamma == 3] <- round(1/(.6*5+1),2)
  statsi$gamma[statsi$gamma == 6] <- round(1/(.6*10+1),2)
  statsi$gamma[statsi$gamma == 1] <- round(1/(.6*2.5+1),2)
  statsi$gamma[statsi$gamma == 4] <- round(1/(.6*7.5+1),2)
  statsi$gamma <- factor(statsi$gamma, c(1,round(1/6,2), round(1/11,2)))
  plott <- ggplot(data=statsi, aes(x=conn05, y=meanvol, colour=gamma))+
    geom_point(size=3,)+
    geom_line(linetype = "dashed", alpha = .6, show.legend = FALSE)+
    ylim(0,.5)+
    #scale_color_manual(values=cols)+
    labs(y= "volume", x = TeX("\\kappa"), title =paste("S=",as.character(S), sep=""), color=TeX("\\gamma"))+
    theme(title=element_text(size=15,face="plain"),
          axis.title=element_text(size=18,face = "plain"), 
          axis.text = element_text(size = 11.5),
          legend.text = element_text(size=11.5),
          legend.title=element_text(size=19), aspect.ratio = 1)
  return(plott)
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
    labs(y= "volume", x = TeX("\\gamma"), title =paste("S=",as.character(S), sep=""), color=TeX("\\kappa"))+
    theme(title=element_text(size=15,face="plain"),
          axis.title=element_text(size=18,face = "plain"), 
          axis.text = element_text(size = 11.5),
          legend.text = element_text(size=11.5),
          legend.title=element_text(size=19), aspect.ratio = 1)
  return(plott)
}

plotDiffS <-function(statsi, whichcols){ 
  plott <- ggplot(data=statsi, aes(x=conn05, y=meanvol, colour=S))+
    geom_point(aes(x=conn05, y=meanvol, shape=S),size=3)+
    scale_shape_manual(values=shapesW[whichcols])+
    geom_line(alpha=.3, show.legend = FALSE, linetype="longdash")+
    ylim(0,.25)+
    scale_color_manual(values=cols[whichcols])+
    labs(y= "volume", x = TeX("\\kappa"))+
    theme(title=element_text(size=15,face="plain"),
          axis.title=element_text(size=18,face = "plain"), 
          axis.text = element_text(size = 11.5),
          legend.text = element_text(size=11.5),
          legend.title=element_text(size=19),aspect.ratio = 1)
  return(plott)
}

plotDiffSZ <-function(statsi, whichcols){ 
  statsi$Z<-as.character(statsi$Z)
  plott <- ggplot(data=statsi, aes(x=conn05, y=meanvol,colour=Z))+
    geom_point(aes(x=conn05, y=meanvol, colour=Z,shape=S),size=3)+
    #scale_shape_manual(values=shapesW[whichcols])+
    geom_line(aes(x=conn05, y=meanvol, colour=Z,shape=S),alpha=.3, show.legend = FALSE, linetype="longdash")+
    ylim(0,.25)+
    #scale_color_manual(values=cols[whichcols])+
    labs(y= "volume", x = TeX("\\kappa"))+
    theme(title=element_text(size=15,face="plain"),
          axis.title=element_text(size=18,face = "plain"), 
          axis.text = element_text(size = 11.5),
          legend.text = element_text(size=11.5),
          legend.title=element_text(size=19),aspect.ratio = 1)
  return(plott)
}


plotDiffZ <-function(statsi, St,xlab){ 
  tit <- paste("=",St,sep="")
  statsi$S<-as.character(statsi$S)
  plott <- ggplot(data=statsi[statsi$Z == St,], aes(x=conn05, y=meanvol, colour=S))+
    geom_point(aes(x=conn05, y=meanvol, shape=S),size=3)+
    scale_shape_manual(values=shapesW[c(1,3,7)])+
    geom_line(alpha=.3, show.legend = FALSE, linetype="longdash")+
    ylim(0,.25)+
    scale_color_manual(values=cols[c(1,3,7)])+
    labs(y= "volume", x = xlab,shape = "P", colour = "P")+
    ggtitle(tit)+
    theme(title=element_text(size=15,face="plain"),
          axis.title=element_text(size=18,face = "plain"), 
          axis.text = element_text(size = 11.5),
          legend.text = element_text(size=11.5),
          legend.title=element_text(size=19),aspect.ratio = 1)
  return(plott)
}

plotDiffZc <-function(statsi, whichcols, xl){ 
  statsi$Z<-as.character(statsi$Z)
  plott <- ggplot(data=statsi, aes(x=conn05, y=meanvol, colour=S))+
    geom_point(aes(x=conn05, y=meanvol, shape=Z),size=3)+
    scale_shape_manual(values=c(15,16,17,18, 9, 12))+
    #geom_line(alpha=.3, show.legend = FALSE, linetype="longdash")+
    ylim(0,.25)+
    scale_color_manual(values=cols[whichcols])+
    labs(y= "volume", x = xl, shape = TeX("\\zeta "), colour = "P")+
    #ggtitle(paste("S=",as.character(St),sep=""))+
    theme(title=element_text(size=15,face="plain"),
          axis.title=element_text(size=18,face = "plain"), 
          axis.text = element_text(size = 11.5),
          legend.text = element_text(size=11.5),
          legend.title=element_text(size=19),aspect.ratio = 1)
  return(plott)
}


plot_treat <- function(datai,diag,treat,whichS){
  datat<-datai[datai$diag==diag,]
  datat$conn05 <- datat[[treat]]
  if (treat == "SC"){
    datat$conn05 <- round(2*datat$conn05)/2
  } else {
    if (treat == "sqSC"){
      datat$conn05 <- round(5*datat$conn05)/5
    }
  }
  statsAll=c()
  for (i in whichS){
    stats=ddply(datat[datat$S==i,], .(conn05), summarize,  
                meanvol=mean(vol), servol=sd(vol)/sqrt(length(vol)))
    stats$S <- as.character(2*i)
    statsAll<-rbind(statsAll, stats)
  }
  plott <- ggplot(data=statsAll, aes(x=conn05, y=meanvol, colour=S))+
    geom_line()+
    geom_point(aes(x=conn05, y=meanvol, shape=S),size=3)+
    #scale_color_manual(values=cols)+
    scale_shape_manual(values=c(3, 6,4,1,7,8,5))+
    #labs(y= "volume", x = TeX("\\kappa S"), title =TeX(paste(" $ \\Delta = $", as.character(diag),sep="")), color="S")+
    labs(y= "volume", x = TeX("\\kappa S"),  color="S")+
    theme(title=element_text(size=15,face="plain"),
          axis.title=element_text(size=18,face = "plain"), 
          axis.text = element_text(size = 11.5),
          legend.text = element_text(size=11.5),
          legend.title=element_text(size=15))
  return(plott)
}

plotNestedNODF <- function(nestAllt,diagt,S){
  datat <- nestAllt[nestAllt$diag==diagt,]
  #datat$nr <- round(2*datat$nest,1)/2
  
  # if (S==3){
  datat$nr <- round(datat$nest,1)
  #}
 
#  if (S==10){
  datat$conn <- round(2*datat$conn,1)/2
 # }
  statsAll=c()
  for (i in unique(datat$nr)){
    datatt <- datat[datat$nr==i,]
    stats=ddply(datatt, .(conn), summarize,  
                meanvol=mean(vol), servol=sd(vol)/sqrt(length(vol)))
    stats$nr <- as.character(i)
    statsAll<-rbind(statsAll, stats)
  }
  
  
  
  
  datatt <- datat[datat$desnest==10,]
  stats=ddply(datatt, .(conn), summarize,  
              meanvol=mean(vol), servol=sd(vol)/sqrt(length(vol)))
  
  #stats$expr <- round(stats$expnest,1)
  
  vt<- c()
  for (i in stats$conn){
    vt<- rbind(vt,statsAll[statsAll$conn == i & statsAll$nr == stats$expr[stats$conn==i],])
  }
  
  if (S==5|S==10){
    statsAlls <- subset(statsAll, nr == "0.2" | nr=="0.4" |  nr=="0.6" | nr=="0.8" | nr=="1")
    vts <- subset(vt,nr == "0.2" | nr=="0.4" |  nr=="0.6" | nr=="0.8" | nr=="1")
  }
  
  if (S==3){
    statsAlls <- statsAll
  }
  # 
  p0 <- statsAll$conn[2]
  p1 <- statsAll$conn[floor(length(unique(statsAll$conn))/4)]
  p2 <- statsAll$conn[floor(length(unique(statsAll$conn))/2)]
  p3 <- statsAll$conn[floor(3*length(unique(statsAll$conn))/4)]
  p4 <- statsAll$conn[length(unique(statsAll$conn))-2]
  # 
  # p1 <- statsAlls$conn[2]
  # p2 <- statsAlls$conn[floor(length(unique(statsAlls$conn))/2)]
  # p3 <- statsAlls$conn[length(unique(statsAlls$conn))-2]
  statsAll$connc <- as.character(statsAll$conn)
  
  
  if (diagt==0){
    
    plott <- ggplot()+
      #geom_line(data=stats, aes(x=conn, y=meanvol),color= "#e166d1",size=5, alpha=0.3,show.legend = FALSE)+
      geom_point(data=statsAll[statsAll$conn == p1 | statsAll$conn == p4 | statsAll$conn == p3 ,], aes(x=nr, y=meanvol, colour=connc, group=connc),size=3)+
      geom_line(data=statsAll[statsAll$conn == p1 | statsAll$conn == p4 | statsAll$conn == p3 ,], aes(x=nr, y=meanvol, colour=connc, group=connc),alpha=.6, show.legend = FALSE, linetype="longdash")+
      ylim(0,.25)+
      scale_color_manual(values=colsN)+
      labs(y= "volume", x = TeX("NODF"), title =paste("S=", as.character(S*2),sep=""), color=TeX("\\kappa"))+
      #labs(y= "volume", x = TeX("\\kappa"), title =TeX(paste(" $ \\Delta = $", as.character(diagt),sep="")), color="NODF")+
      theme(title=element_text(size=15,face="plain"),
            axis.title=element_text(size=18,face = "plain"), 
            axis.text = element_text(size = 11.5),
            legend.text = element_text(size=11.5),
            legend.title=element_text(size=15),
            aspect.ratio=1/1)
  }
  else{
    plott <- ggplot()+
      #geom_line(data=stats, aes(x=conn, y=meanvol),color="gray",size=3, alpha=0.7,show.legend = FALSE)+
      geom_point(data=statsAlls[statsAlls$conn == p1 | statsAlls$conn == p2 | statsAlls$conn == p3 ,], aes(x=nr, y=meanvol, colour=connc, group=connc),size=3)+
      geom_line(data=statsAlls[statsAlls$conn == p1 | statsAlls$conn == p2 | statsAlls$conn == p3 ,], aes(x=nr, y=meanvol, colour=connc),alpha=.6, show.legend = FALSE, linetype="longdash")+
      ylim(0,.5)+
      scale_color_manual(values=colsN)+
      labs(y= "volume", x = TeX("\\kappa"), title ="S=10", color="NODF")+
      #labs(y= "volume", x = TeX("\\kappa"), title =TeX(paste(" $ \\Delta = $", as.character(diagt),sep="")), color="NODF")+
      theme(title=element_text(size=15,face="plain"),
            axis.title=element_text(size=18,face = "plain"), 
            axis.text = element_text(size = 11.5),
            legend.text = element_text(size=11.5),
            legend.title=element_text(size=15),
            aspect.ratio=1/1)
  }
  return(plott)
}
get_statsNest <- function(datat,nest_col){
  datat$nr <- round(datat[[nest_col]],1)
  datat$conn <- round(2*datat$conn,1)/2
  datat$connf <- as.factor(datat$conn)
  
  statsAll=c()
  for (i in unique(datat$nr)){
    datatt <- datat[datat$nr==i,]
    stats=ddply(datatt, .(conn), summarize,  
                meanvol=mean(vol), servol=sd(vol)/sqrt(length(vol)))
    stats$nr <- i
    statsAll<-rbind(statsAll, stats)
  }
  statsAll$connc <- as.character(statsAll$conn)
  return(statsAll)
}

calcNODFc <- function(datat,n){
  #get unique connectances
  connvec <- unique(datat$conn,n)
  #value strorage
  nestMax <- c()
  #aux values for matrix
  vecC <-  abs(rnorm(n*n, mean = 0, sd = 1))
  #get max nodf for every connectance
  for (i in connvec){ 
    #build any non-singular matrix C
    Ct <- matrix(rep(1,4),nrow = 2,ncol=2)
    while (rankMatrix(Ct)[1]<n){
    Ct <- RemoveLinksC(i, vecC, n,n)
    }
    #get maxnodf and store value
    maxth <- maxnodf(Ct,1)$max_nodf
    nestMax<- append(nestMax,maxth)
  }
  #create column
  datat$nodfmax <- 1
  #insert corresponding maxnodf in every line
  for (i in 1:length(datat[[1]])){
    conth <- datat$conn[[i]]
    maxth <- nestMax[]
    datat$nodfmax[[i]] <- nestMax[connvec == conth]
  }
  #calculate nodfc
  datat$nodfc <- (datat$nest)/(datat$conn* log(n) *datat$nodfmax)
  return(datat)
}

plotNested <- function(nestAllt,diagt,S){
  datat <- nestAllt[nestAllt$diag==diagt,]
  datat <- na.omit(datat)
  
  # if (S==3){
  datat$nr <- round(datat$nest,1)
  #}
  
  #if (S==10){
  datat$conn <- round(2*datat$conn,1)/2
  #}
  statsAll=c()
  for (i in unique(datat$nr)){
    datatt <- datat[datat$nr==i,]
    stats=ddply(datatt, .(conn), summarize,  
                meanvol=mean(vol), servol=sd(vol)/sqrt(length(vol)))
    stats$nr <- as.character(i)
    statsAll<-rbind(statsAll, stats)
  }
  
  
  
  
  datatt <- datat[datat$desnest==10,]
  stats=ddply(datatt, .(conn), summarize,  
              meanvol=mean(vol), servol=sd(vol)/sqrt(length(vol)))
  
  #stats$expr <- round(stats$expnest,1)
  
  vt<- c()
  for (i in stats$conn){
    vt<- rbind(vt,statsAll[statsAll$conn == i & statsAll$nr == stats$expr[stats$conn==i],])
  }
  
  #if (S==5|S==10){
  statsAlls <- subset(statsAll, nr == "0.2" | nr == "0.3" |nr == "0.5" |nr == "0.7" |nr == "0.9" |    nr=="0.4" |  nr=="0.6" | nr=="0.8" | nr=="1")
  vts <- subset(vt,nr == "0.2" | nr=="0.4" |  nr=="0.6" | nr=="0.8" | nr=="1")
  #}
  
  # if (S==3){
  # statsAlls <- statsAll
  # }
  # 
  if (diagt==0){
    plott <- ggplot()+
      #geom_line(data=stats, aes(x=conn, y=meanvol),color= "#e166d1",size=5, alpha=0.3,show.legend = FALSE)+
      geom_point(data=statsAlls, aes(x=conn, y=meanvol, colour=nr),size=3)+
      geom_line(data=statsAlls, aes(x=conn, y=meanvol, colour=nr),alpha=.6, show.legend = FALSE, linetype="longdash")+
      ylim(0,.25)+
      # scale_color_manual(values=colsN)+
      labs(y= "volume", x = TeX("\\kappa"), title =paste("S=", as.character(S*2),sep=""), color="NODF")+
      #labs(y= "volume", x = TeX("\\kappa"), title =TeX(paste(" $ \\Delta = $", as.character(diagt),sep="")), color="NODF")+
      theme(title=element_text(size=15,face="plain"),
            axis.title=element_text(size=18,face = "plain"), 
            axis.text = element_text(size = 11.5),
            legend.text = element_text(size=11.5),
            legend.title=element_text(size=15),
            aspect.ratio=1/1)
  }else{
    plott <- ggplot()+
      #geom_line(data=stats, aes(x=conn, y=meanvol),color="gray",size=3, alpha=0.7,show.legend = FALSE)+
      geom_point(data=statsAlls, aes(x=conn, y=meanvol, colour=nr),size=3)+
      geom_line(data=statsAlls, aes(x=conn, y=meanvol, colour=nr),alpha=.5, show.legend = FALSE)+
      ylim(0,.25)+
      scale_color_manual(values=colsN)+
      labs(y= "volume", x = TeX("\\kappa"), title ="S=10", color="NODF")+
      #labs(y= "volume", x = TeX("\\kappa"), title =TeX(paste(" $ \\Delta = $", as.character(diagt),sep="")), color="NODF")+
      theme(title=element_text(size=15,face="plain"),
            axis.title=element_text(size=18,face = "plain"), 
            axis.text = element_text(size = 11.5),
            legend.text = element_text(size=11.5),
            legend.title=element_text(size=15),
            aspect.ratio=1/1)
  }
  return(plott)
}

# data fit----

make_logs <- function(datat){
  datat$logV <-  log(1/datat$meanvol,10)
  datat$logSC <- log(datat$SC,10)
  return(datat)
}

data_fit<-function(datat){
  ls0 <- lm(logSC~logV, data=datat)
  #summary(ls0)
  #confint(ls0)
  a=10^ls0$coefficients[1]
  p=ls0$coefficients[2]
  return(c(a,p))
}

make_fitted_data <- function(coef, nv){
  a<-coef[1]
  p<-coef[2]
  x=seq(.1,20,.01)
  fitted_d = data.frame(x=x)
  
  fitted_d$y <-  1/(a*(x^p))
  
  fitted_dS=c()
  xS=seq(.05,1,.001)
  for (i in nv){
    fitted_dSt = data.frame(x=xS)
    fitted_dSt$y <- 1/(a*((i*xS)^p))
    fitted_dSt$S <- as.character(i)
    fitted_dS <- rbind(fitted_dS,fitted_dSt)}
  return(fitted_dS)
}



# 
# 
# colsN =brewer.pal(8, "Dark2")
# shapesW = c(3, 1,6,4,7,8,5)
# #shapesW <- c(3,4,6,1,7,8,5)
