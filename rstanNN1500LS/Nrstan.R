#clear screen
rm(list=ls()) 

##set working dir
setwd("C:/Users/dell/Desktop/rstan/NN1500LS")

##used packages
library(MASS)
library(rstan)

# compile the stan model
stan_m <- stan_model(model_code = "
data {
int<lower=1> N; // sample size
int<lower=1> J; // number of items
int Y[N,J]; // data matrix of order [N,J]
int<lower=1> P; // number of covariates 
matrix[N,P] X;
real Nmupr;
real<lower=0> Nsigpr;
real alpha_a;
real alpha_b;
}

parameters {
vector[J] inu; // baseline intercept
vector<lower=0>[J] lam; // baseline slope
matrix[J,P] kap;  // intercept DIF coefficients
matrix[J,P] lom;   // slope DIF coefficients
vector[N] fac_dist_helper; // helper for non-centered sampling
vector[P] mu_imp; // factor mean impact
vector[P] phi_imp; // factor sd impact
matrix<lower=0>[J,P] alpha_i; // Laplace variance on intercept DIF
matrix<lower=0>[J,P] alpha_l; // Laplace variance on slope DIF
}

transformed parameters {
vector[N] fac_scor; // person factor scores
matrix[N,J] mu;

// build non-centered factor score 
for(i in 1:N){
  fac_scor[i] <- 0 + X[i]*mu_imp + exp(0.5*(X[i]*phi_imp)) * fac_dist_helper[i];
  for(j in 1:J){
    mu[i,j] <- inu[j] + X[i]*kap[j]' + (lam[j] + X[i]*lom[j]') * fac_scor[i];}
}
}

model {
// the priors
to_vector(alpha_i) ~ gamma(alpha_a, alpha_b);
to_vector(alpha_l) ~ gamma(alpha_a, alpha_b);
inu ~ normal(Nmupr, Nsigpr);
lam ~ normal(Nmupr, Nsigpr);
fac_dist_helper ~ normal(0, 1);
to_vector(kap) ~ double_exponential(0, 1 ./ sqrt(to_vector(alpha_i)));
to_vector(lom) ~ double_exponential(0, 1 ./ sqrt(to_vector(alpha_l)));

mu_imp ~ normal(Nmupr, Nsigpr);
phi_imp ~ normal(Nmupr, Nsigpr);

// The likelihood
for(j in 1:J){
  Y[ ,j] ~ bernoulli_logit(mu[ ,j]);
}
}
", verbose = T )


##true values
set.seed(10)
N=1500
J=10
P=2
CNUM=100

nu0=rep(0,J)
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
lam0=rep(1,J)   
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

##replications
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
      lom = matrix((rep(0.1, (J*P)) + runif((J*P),0,.1) ), nrow=J, ncol=P),
      mu_imp = rep(0.1, P) + runif(P,0,.1),
      phi_imp = rep(0.1, P)  + runif(P,0,.1) )
    return(init.values);
  }

  stan_ssp <- sampling(stan_m, data = fa.data, pars =c("inu","kap","lam","lom", "alpha_i","alpha_l", "mu_imp", "phi_imp"), chains = 3, iter = 5000, init = init_model, cores = 3)

  aa <- summary(stan_ssp, probs = c(0.025, 0.975), pars = c("inu","kap","lam","lom", "mu_imp", "phi_imp") )

  write(get_num_divergent(stan_ssp), file="diver.txt", ncol=1, append=TRUE, sep="\t")
  write(aa$summary[,c(7)], file="Rhh.txt", ncol=length(aa$summary[,c(7)]), append=TRUE, sep="\t")
  write(aa$summary[,c(1)], file="amean.txt", ncol=length(aa$summary[,c(1)]), append=TRUE, sep="\t")
  write(aa$summary[,c(3)], file="asd.txt", ncol=length(aa$summary[,c(3)]), append=TRUE, sep="\t")
  write(aa$summary[,c(4)], file="ainl.txt", ncol=length(aa$summary[,c(4)]), append=TRUE, sep="\t")
  write(aa$summary[,c(5)], file="ainr.txt", ncol=length(aa$summary[,c(5)]), append=TRUE, sep="\t")

  print(CIR)
  et<-proc.time()
  print((et-bt)[3])
  write((et-bt)[3], file="time.txt", ncol=1, append=TRUE, sep="\t")
}

date()
save.image(paste("Rstan",".RData",sep=""))





