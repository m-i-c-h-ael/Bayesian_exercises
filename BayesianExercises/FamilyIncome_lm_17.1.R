#17.1 Compare models for family income as a function of family size via linear model
 #w/o and with quadratic term
  #based on median income from different states

rm(list=ls())
cat('\014')
graphics.off()

library('rjags')

modStr=' #transform data

model{
  for(i in 1:Ntotal){
    zy[i] ~ dt(mu[i],1/sig^2,nu)   #common sig and nu for all states
    mu[i]= zbeta0[s[i]] + zbeta1[s[i]] * zx[i]
  }
  
  for(j in 1:Nstates){
    zbeta0[j] ~ dnorm(mu0,1/sig0^2)
    zbeta1[j] ~ dnorm(mu1,1/sig1^2)
  }
  
  mu0 ~ dnorm(M,1/S^2)
  sig0 ~ dunif(L,H)
  mu1 ~ dnorm(M,1/S^2)
  sig1 ~ dunif(L,H)
  
  nuMin1 ~ dexp(1/29.0)
  nu= 1 + nuMin1
  sig ~ dunif(L,H)
  
  M= 0
  S= 10
  L= 1/1000
  H= 1000
  
  #Retransform
}'