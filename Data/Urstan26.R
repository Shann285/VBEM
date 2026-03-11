#clear screen
rm(list=ls()) 

##set working dir
setwd("C:/Users/dell/Desktop/Data")

##used packages
library(MASS)
library(rstan)

bt<-proc.time()
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
vector[N] fac_dist_helper; // helper for non-centered sampling
vector[P] mu_imp; // factor mean impact
vector[P] phi_imp; // factor sd impact
matrix<lower=0>[J,P] alpha_i; // Laplace variance on intercept DIF
}

transformed parameters {
vector[N] fac_scor; // person factor scores
matrix[N,J] mu;

// build non-centered factor score 
for(i in 1:N){
  fac_scor[i] <- 0 + X[i]*mu_imp + exp(0.5*(X[i]*phi_imp)) * fac_dist_helper[i];
  for(j in 1:J){
    mu[i,j] <- inu[j] + X[i]*kap[j]' + lam[j] * fac_scor[i];}
}
}

model {
// the priors
to_vector(alpha_i) ~ gamma(alpha_a, alpha_b);
inu ~ normal(Nmupr, Nsigpr);
lam ~ normal(Nmupr, Nsigpr);
fac_dist_helper ~ normal(0, 1);
to_vector(kap) ~ double_exponential(0, 1 ./ sqrt(to_vector(alpha_i)));

mu_imp ~ normal(Nmupr, Nsigpr);
phi_imp ~ normal(Nmupr, Nsigpr);

// The likelihood
for(j in 1:J){
  Y[ ,j] ~ bernoulli_logit(mu[ ,j]);
}
}
", verbose = T )


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


##model input and starting values
set.seed(10)
fa.data <- list(N = N, J = J, Y = y, P = P, X = x, Nmupr = 0, Nsigpr = 2, alpha_a = 10, alpha_b = 1)

init_model = function(){
  init.values <- list( 
    inu = rep(0.1, J) + runif(J,0,.1),
    lam = rep(0.1, J) + runif(J,0,.1),
    kap = matrix((rep(0.1, (J*P)) + runif((J*P),0,.1) ), nrow=J, ncol=P),
    mu_imp = rep(0.1, P) + runif(P,0,.1),
    phi_imp = rep(0.1, P)  + runif(P,0,.1) )
  return(init.values);
}

stan_ssp <- sampling(stan_m, data = fa.data, pars =c("inu","kap","lam", "alpha_i","mu_imp", "phi_imp"), chains = 3, iter = 5000, init = init_model, cores = 3)

aa <- summary(stan_ssp, probs = c(0.025, 0.975), pars = c("inu","kap","lam","mu_imp", "phi_imp") )

et<-proc.time()
print((et-bt)[3])


get_num_divergent(stan_ssp)

sum(aa$summary[,c(7)]>=1.05)  

round(aa$summary[,c(1)],3)

save.image(paste("DataURstan",".RData",sep=""))




