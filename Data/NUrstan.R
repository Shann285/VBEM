
rm(list=ls()) #clear screen

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
int<lower=0, upper=1> Y[N*J]; // data vector of order N*J
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

model {
vector[N] fac_scor = X*mu_imp + exp(0.5*(X*phi_imp)) .* fac_dist_helper;
matrix[J,N] mu = rep_matrix(inu, N) + (kap*X') + (lam*fac_scor') + lom * (X .* rep_matrix(fac_scor, P))';

// the priors
to_vector(alpha_i) ~ gamma(alpha_a, alpha_b);
to_vector(alpha_l) ~ gamma(alpha_a, alpha_b);
inu ~ normal(Nmupr, Nsigpr);
lam ~ normal(1, 0.5);
fac_dist_helper ~ normal(0, 1);
to_vector(kap) ~ double_exponential(0, 1 ./ sqrt(to_vector(alpha_i)));
to_vector(lom) ~ double_exponential(0, 1 ./ sqrt(to_vector(alpha_l)));

mu_imp ~ normal(Nmupr, Nsigpr);
phi_imp ~ normal(Nmupr, Nsigpr);

// The likelihood
Y ~ bernoulli_logit(to_vector(mu));
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
x <- as.matrix(cbind(data11[,2], sd(data11[,2])*scale(data11[,1])))

N=dim(y)[1]
J=dim(y)[2]
P=dim(x)[2]


##model input and starting values
set.seed(10)
fa.data <- list(N = N, J = J, Y = as.vector(t(y)), P = P, X = x, Nmupr = 0, Nsigpr = 2, alpha_a = 10, alpha_b = 1)

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

#stan_ssp <- sampling(stan_m, data = fa.data, pars =c("inu","kap","lam","lom", "mu_imp", "phi_imp"), chains = 3, iter = 5000, warmup = 2500, init = init_model, cores = 3)
stan_ssp <- sampling(stan_m, data = fa.data, pars =c("inu","kap","lam","lom", "mu_imp", "phi_imp"), chains = 3, iter = 6000, warmup = 3000, init = init_model, cores = 3, 
           control = list(adapt_delta=0.95, max_treedepth=12))

aa <- summary(stan_ssp, probs = c(0.025, 0.50, 0.975), pars = c("inu","kap","lam","lom", "mu_imp", "phi_imp") )

et<-proc.time()
print((et-bt)[3])


get_num_divergent(stan_ssp)
sum(aa$summary[,c(8)]>=1.05)  
sum(aa$summary[,c(7)]<1000)  

round(aa$summary[,c(1)],3)

save.image(paste("DataNURstan922",".RData",sep=""))



