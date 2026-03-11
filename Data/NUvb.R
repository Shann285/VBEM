#clear screen
rm(list=ls())

##set working dir
setwd("C:/Users/dell/Desktop/Data")

##used packages
library(MASS)
library(statmod)
library(rstan)

bt <- proc.time()
##user-defined function of tau
tauksi <- function(ksi){
   aa=ksi
   aa[ksi==0] <- 1/8
   aa[ksi!=0] <- (1/(2*aa[ksi!=0]))*(1/(1+exp(-aa[ksi!=0]))-1/2)
   return(aa)
}

##Bayesian GLM with gamma response
stanm <- stan_model(model_code = "
data {
  int<lower=0> N;      // sample size
  int<lower=0> P;      // dimension of beta
  vector[N] w;         // response
  matrix[N, P] x;      // design matrix
  real Nmupr;
  real<lower=0> Nsigpr;
}

parameters {
  vector[P] beta;     //regression coefficients
}

transformed parameters {
  vector[N] mu;
  for (i in 1:N) {
    mu[i] = exp(x[i,]*beta);
  }
}

model {
  beta ~ normal(Nmupr, Nsigpr);
  // gamma distribution(shape=0.5)
  for (i in 1:N) {
    w[i] ~ gamma(0.5, 0.5 / mu[i]);  
  }
}
", verbose=T )


##true values
data1 <- read.csv(file="b21f.csv", header = TRUE)
head(data1)
dim(data1)

data11 <- data1[(data1$age>=16)&(data1$age<=25),]
head(data11)
dim(data11)

y <- as.matrix(data11[,c(3:21)])
x <- as.matrix(cbind(data11[,2],scale(data11[,1])))

N=dim(y)[1]
J=dim(y)[2]
P=dim(x)[2]
NMCMC=2000

##prior values
set.seed(10)
Nmupr=0
Nsigpr=2
Nsig2pr=2^2
MNmupr=rep(0,P)
MNsigpr=diag(2^2,P)
apr=10
bpr=1


##initial values
munu <- rep(0.1,J)
sig2nu <- rep(0.1,J)
mukappa <- matrix(0.1,nrow=J,ncol=P)
sig2kappa <- matrix(0.1,nrow=J,ncol=P)
mukdel <- matrix(0.1,nrow=J,ncol=P)
mukalp <- matrix(0.1,nrow=J,ncol=P)
mulam <- rep(0.1,J)
sig2lam <- rep(0.1,J)
muomega <- matrix(0.1,nrow=J,ncol=P)
sig2omega <- matrix(0.1,nrow=J,ncol=P)
muodel <- matrix(0.1,nrow=J,ncol=P)
muoalp <- matrix(0.1,nrow=J,ncol=P)
ksi <- y
mueta <- rnorm(N,0,1)   
sig2eta <- rep(0.1,N)
mugamma <- rep(0.1,P)             
sig2gamma <- diag(0.1,P)                                 
mubeta <- rep(0.1,P)           
sig2beta <- diag(0.1,P)                                
M <- 1000                               

##VB iteration
for(g in 1:NMCMC){
  ppold <- 1/(1+exp(-1*( matrix(rep(munu,each=N),nrow=N)+ x%*%t(mukappa)+mueta%*%t(mulam)+(x%*%t(muomega))*matrix(rep(mueta,J),nrow=N) ) )) 

  sig2eta <- 1/( 2*rowSums(tauksi(ksi)*( matrix(rep((sig2lam+mulam^2),N),nrow=N,byr=TRUE)+2*matrix(rep(mulam,N),nrow=N,byr=TRUE)*(x%*%t(muomega))+(x^2)%*%t(sig2omega)+(x%*%t(muomega))^2 ))
                  + 1/exp(x%*%mubeta-0.5*diag(x%*%sig2beta%*%t(x))) )
  mueta <- sig2eta*( rowSums((y-0.5)*(matrix(rep(mulam,N),nrow=N,byr=TRUE)+x%*%t(muomega)))-2*rowSums( tauksi(ksi)*(matrix(rep(munu,N),nrow=N,byr=TRUE)+x%*%t(mukappa))*(matrix(rep(mulam,N),nrow=N,byr=TRUE)
                 +x%*%t(muomega))) + (x%*%mugamma)/exp(x%*%mubeta-0.5*diag(x%*%sig2beta%*%t(x))) ) 

  sig2gamma <- chol2inv(chol(  t(x)%*%diag(as.vector(1/exp(x%*%mubeta-0.5*diag(x%*%sig2beta%*%t(x)))))%*%x +  chol2inv(chol(MNsigpr)) ))
  mugamma <-  sig2gamma%*%(  t(x)%*%diag(as.vector(1/exp(x%*%mubeta-0.5*diag(x%*%sig2beta%*%t(x)))))%*%mueta + chol2inv(chol(MNsigpr))%*%MNmupr )

  yhat <- (sig2eta+mueta^2)-2*mueta*(x%*%mugamma)+diag(x%*%(sig2gamma+mugamma%*%t(mugamma))%*%t(x))
  dbeta2 <- list(N=N,P=P,w=as.vector(yhat),x=x,Nmupr=Nmupr,Nsigpr=Nsigpr)
  fit_opt <- optimizing(
    stanm, 
    data=dbeta2,                   
    hessian = TRUE                      
  )
  mubeta <- fit_opt$par[1:P]
  sig2beta <- chol2inv(chol( t(x)%*%diag(as.vector(0.5*yhat*exp(-x%*%mubeta)))%*%x + chol2inv(chol(MNsigpr)) ))

  for(j in 1:J){
    sig2nu[j] <- 1/(2*sum(tauksi(ksi[,j])) +1/Nsig2pr)
    munu[j] <- sig2nu[j]*( sum(y[,j]-0.5)-2*sum( tauksi(ksi[,j])*(x%*%mukappa[j,]+mulam[j]*mueta+(x%*%muomega[j,])*mueta) )+ Nmupr/Nsig2pr )

    sig2lam[j] <- 1/(2*sum(tauksi(ksi[,j])*(sig2eta+mueta^2)) + 1/Nsig2pr)
    mulam[j] <- sig2lam[j]*( sum((y[,j]-0.5)*mueta)- 2*sum( tauksi(ksi[,j])*((munu[j]+x%*%mukappa[j,])*mueta+(x%*%muomega[j,])*(sig2eta+mueta^2))) + 6/Nsig2pr )

    for(pp in 1:P){
      sig2kappa[j,pp] <- 1/( 2*sum(tauksi(ksi[,j])*(x[,pp]^2))+ mukdel[j,pp] )
      mukappa[j,pp] <- sig2kappa[j,pp]*( sum((y[,j]-0.5)*x[,pp])-2*sum(tauksi(ksi[,j])*x[,pp]*(munu[j]+x[,-pp]*mukappa[j,-pp]+mulam[j]*mueta+x%*%muomega[j,]*mueta)) ) 
      mukdel[j,pp] <- sqrt(mukalp[j,pp]/(sig2kappa[j,pp]+mukappa[j,pp]^2))
      atau2 <- 1/rinvgauss(M,mean=mukdel[j,pp],shape=mukalp[j,pp])
      mukalp[j,pp] <- (apr+1)/(bpr+0.5*mean(atau2))

      sig2omega[j,pp] <- 1/( 2*sum(tauksi(ksi[,j])*( (x[,pp]^2)*(sig2eta+mueta^2) ))+ muodel[j,pp] )
      muomega[j,pp] <- sig2omega[j,pp]*( sum((y[,j]-0.5)*(x[,pp]*mueta))-2*sum(tauksi(ksi[,j])*x[,pp]*(munu[j]*mueta+(x%*%mukappa[j,])*mueta+mulam[j]*(sig2eta+mueta^2)+(x[,-pp]*muomega[j,-pp])*(sig2eta+mueta^2) )) ) 
      muodel[j,pp] <- sqrt(muoalp[j,pp]/(sig2omega[j,pp]+muomega[j,pp]^2))
      btau2 <- 1/rinvgauss(M,mean=muodel[j,pp],shape=muoalp[j,pp])
      muoalp[j,pp] <- (apr+1)/(bpr+0.5*mean(btau2))
    }
  }

  for(i in 1:N){
    for(j in 1:J){
      ksi[i,j] <- sqrt( (sig2nu[j]+munu[j]^2)+2*munu[j]*mulam[j]*mueta[i]+(sig2lam[j]+mulam[j]^2)*(sig2eta[i]+mueta[i]^2) +(t(x[i,])%*%(diag(sig2omega[j,])+muomega[j,]%*%t(muomega[j,]))%*%x[i,])*(sig2eta[i]+mueta[i]^2)+
                    2*munu[j]*(mukappa[j,]%*%x[i,])+2*mulam[j]*mueta[i]*(mukappa[j,]%*%x[i,]) + t(x[i,])%*%(diag(sig2kappa[j,])+mukappa[j,]%*%t(mukappa[j,]))%*%x[i,] +
                    2*(muomega[j,]%*%x[i,])*(munu[j]*mueta[i]+mulam[j]*(sig2eta[i]+mueta[i]^2)+mukappa[j,]%*%x[i,]*mueta[i]) )
      if(y[i,j]==0){ksi[i,j] <- -ksi[i,j]}
    }
  }
  
  L2new <- c(munu,as.vector(t(mukappa)),mulam, as.vector(t(muomega)),mugamma,mubeta, sig2nu,as.vector(sig2kappa),sig2lam,as.vector(sig2omega),as.vector(mukdel),as.vector(mukalp),as.vector(muodel),as.vector                                       (muoalp),sig2gamma,sig2beta)

  ppnew <- 1/(1+exp(-1*( matrix(rep(munu,each=N),nrow=N)+ x%*%t(mukappa)+mueta%*%t(mulam)+(x%*%t(muomega))*matrix(rep(mueta,J),nrow=N) ) ))

  if( max(abs(ppnew-ppold))<0.0001 ){break}
}


et<-proc.time()
print((et-bt)[3])

round(L2new[1:(2*J+2*J*P+2*P)],3)

save.image(paste("DatavbNU",".RData",sep=""))




