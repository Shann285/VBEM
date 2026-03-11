# VBEM
An efficient regularized variational Bayesian estimation for assessing differential item functioning in moderated nonlinear factor analysis (MNLFA)

MNLFA has emerged as a significant and flexible psychometric model for examining measurement invariance (MI) and differential item functioning (DIF). 
Since MNLFA models can accommodate various types of response variables, with model parameters moderated by a range of exogenous covariates, model estimation becomes a critical issue. 
Recent research has explored the use of Markov chain Monte Carlo (MCMC) estimation for MNLFA. While MCMC
estimation provides statistical and practical advantages over likelihood-based metheods, it is computationally demanding and time-consuming. 
This paper introduces an efficient regularized variational Bayesian expectation-maximization (VBEM) algorithm to accelerate the estimation of MNLFA models.

The Data file folder includes five files. Nrstan26.R implements the MCMC method for the nonuniform model in the real data. 
NUvb.R implements our VBEM method for the nonuniform model in the real data. Urstan26.R implements the MCMC method for the uniform model in the real data. Uvb.R implements our VBEM method for the uniform model in the real data. b21d.csv is the final data file.

The rstanNN1500LS file folder includes the regularized MCMC estimation procedure for the nonuniform model under the condition of N=1500, 20% DIF and small DIF.

The rstanUN1500LS file folder includes the regularized MCMC estimation procedure for the nuniform model under the condition of N=1500, 20% DIF and small DIF.

The VBEMNN1500LS file folder includes the regularized VBEM estimation procedure for the nonuniform model under the condition of N=1500, 20% DIF and small DIF.

The VBEMUN1500LS file folder includes the regularized VBEM estimation procedure for the uniform model under the condition of N=1500, 20% DIF and small DIF.
