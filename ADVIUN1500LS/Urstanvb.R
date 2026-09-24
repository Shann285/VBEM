
rm(list=ls()) #clear screen

##set working dir
setwd("C:/Users/dell/Desktop/ADVIUN1500LS")

##used packages
library(MASS)
library(cmdstanr)

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
lam0=c(0.6,0.8,1.0,1.2,1.4,0.6,0.8,1.0,1.2,1.4) 

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
      pp <- exp((nu0[j]+kappa0[j,]%*%x[i,])+lam0[j]*eta[i])/(1+exp((nu0[j]+kappa0[j,]%*%x[i,])+lam0[j]*eta[i]))
      y[i,j] <- rbinom(1, size=1, prob=pp)
    }
  }
  write(y, file="Y.txt", ncol=dim(y)[1], append=T)
  write(x, file="X.txt", ncol=dim(x)[1], append=T)
  print(CIR)
}

mod <- cmdstan_model("MU.stan")

for(CIR in 1:CNUM){
  bt <- proc.time()
  y <- matrix(0, nrow=N, ncol=J)
  x <- matrix(0, nrow=N, ncol=P)
  y <- matrix(scan("Y.txt", skip=(CIR-1)*J, nlines=J), nrow=N, ncol=J)
  x <- matrix(scan("X.txt", skip=(CIR-1)*P, nlines=P), nrow=N, ncol=P)

  # model input and starting values
  fa.data <- list(N = N, J = J, Y = y, P = P, X = x, Nmupr = 0, Nsigpr = 2, alpha_a = 10, alpha_b = 1)

  init_model = function(){
    init.values <- list( 
      inu = rep(0.1, J) + runif(J,0,.1),
      lam = rep(0.1, J) + runif(J,0,.1),
      kap = matrix((rep(0.1, (J*P)) + runif((J*P),0,.1) ), nrow=J, ncol=P),
      mu_imp = rep(0.1, P) + runif(P,0,.1),
      phi_imp = rep(0.1, P)  + runif(P,0,.1) )
    return(init.values)
  }

  # Fit the model using CmdStan
  stan_ssp <- mod$variational(
    data = fa.data,
    init =  init_model,
    grad_samples = 10,
    elbo_samples = 500,
    tol_rel_obj = 0.001
  )

  # Generate summary and convergence diagnostics
  aa <- stan_ssp$summary(
     variables = c("inu","kap","lam", "mu_imp", "phi_imp"),
     posterior::default_summary_measures()[1:4],
     quantiles = ~ posterior::quantile2(., probs = c(0.025, 0.975)),
     posterior::default_convergence_measures()
  )

  write(aa$rhat, file="Rhh.txt", ncol=nrow(aa), append=TRUE, sep="\t")
  write(aa$mean, file="amean.txt", ncol=nrow(aa), append=TRUE, sep="\t")
  write(aa$sd, file="asd.txt", ncol=nrow(aa), append=TRUE, sep="\t")
  write(aa$`q2.5`, file="ainl.txt", ncol=nrow(aa), append=TRUE, sep="\t")
  write(aa$`q97.5`, file="ainr.txt", ncol=nrow(aa), append=TRUE, sep="\t")

  print(CIR)
  et <-proc.time()
  print((et-bt)[3])
  write((et-bt)[3], file="time.txt", ncol=1, append=TRUE, sep="\t")
}

date()

save.image(paste("Rstanvb",".RData",sep=""))
