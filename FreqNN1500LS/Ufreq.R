
rm(list=ls()) #clear screen

##set working dir
setwd("C:/Users/dell/Desktop/FreqNN1500LS")

library(MASS)
library(regDIF)


##true values
set.seed(10)
N=1500
J=10
P=2
CNUM=100

nu0=c(1.0,0.5,0,-0.5,-1.0,1.0,0.5,0,-0.5,-1.0) 
kappa0=matrix(c(
0,0,
0,0,
0,0,
0,0,
0,0,
0,0,
0,0,
0.15,0,
0,0.15,
0,0),nrow=J,ncol=P,byr=TRUE)
lam0=c(0.6,0.8,1.0,1.2,1.4, 0.6,0.8,1.0,1.2,1.4)  
omega0=matrix(c(
0,0,
0,0,
0,0,
0,0,
0,0,
0,0,
0,0,
-0.1,0,
0,-0.1,
0,0),nrow=J,ncol=P,byr=TRUE)      
gamma0=c(0.5,0)
beta0=c(0,0.2)


##generated data
for(CIR in 1:CNUM){
  x1 <- rbinom(N, size=1, prob=0.5)          
  x2 <- rnorm(N, mean=0, sd=0.5) 
  x <- cbind(x1,x2)
  eta <- rep(0,N)
  y <- matrix(0,nrow=N,ncol=J)

  for(i in 1:N){
    mut <- x[i,]%*%gamma0
    sig2t <- exp(x[i,]%*%beta0)
    eta[i] <- rnorm(1, mut, sqrt(sig2t))
    for(j in 1:J){
      pp <- exp((nu0[j]+kappa0[j,]%*%x[i,])+(lam0[j]+omega0[j,]%*%x[i,])*eta[i])/(1+exp((nu0[j]+kappa0[j,]%*%x[i,])+(lam0[j]+omega0[j,]%*%x[i,])*eta[i]))
      y[i,j] <- rbinom(1, size=1, prob=pp)
    }
  }
  write(y, file="Y.txt", ncol=dim(y)[1], append=T)
  write(x, file="X.txt", ncol=dim(x)[1], append=T)
  print(CIR)
}


for(CIR in 1:CNUM){
  bt <- proc.time()
  y <- matrix(0, nrow=N, ncol=J)
  x <- matrix(0, nrow=N, ncol=P)
  y <- matrix(scan("Y.txt", skip=(CIR-1)*J, nlines=J), nrow=N, ncol=J)
  x <- matrix(scan("X.txt", skip=(CIR-1)*P, nlines=P), nrow=N, ncol=P)

  aa <- regDIF(item.data=y, pred.data=x, pen.type="lasso")
  
  op <- which( aa$bic == min(na.omit(aa$bic)) )

  res <- c(coef(aa)$base[1:J,op], coef(aa)$dif[1:(2*J),op], coef(aa)$base[(J+1):(2*J),op], coef(aa)$dif[(2*J+1):(4*J),op], coef(aa)$impact[,op])
  
  write(res, file="res.txt", ncol=length(res), append=TRUE, sep="\t")

  print(CIR)
  et<-proc.time()
  print((et-bt)[3])
  write((et-bt)[3], file="time.txt", ncol=1, append=TRUE, sep="\t")
}

date()

save.image(paste("Freq",".RData",sep=""))



