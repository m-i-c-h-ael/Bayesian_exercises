# Model with categorical predictor and metric response (corresponding to ANOVA)

rm(list=ls())
cat('\014')
graphics.off()

library('rjags')
source('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/DBDA2E-utilities.R')

modStr='
data{
  yMean= mean(y)
  ySD= sd(y)
}
model{
  for(i in 1:Ntotal) { y[i] ~ dnorm(alpha0 + alpha[group[i]], 1/sig_y^2)  }
  sig_y ~ dunif(ySD/100,ySD*10)
  alpha0 ~ dnorm(yMean,1/(5*ySD)^2)
  
  for(j in 1:Ngroups){
    alpha[j] ~ dnorm(0.0,1/sig_beta^2)
  }
  sig_beta ~ dgamma(shape,scale)
  shape= gammaPar[1]
  scale= gammaPar[2]
  
  # Ensure total deviance is 0
  for(j in 1:Ngroups){
    mu[j]= alpha0 + alpha[j]
  }
  beta0= mean(mu[1:Ngroups])
  for(j in 1:Ngroups){
    beta[j]= mu[j]-beta0
  }
  
  # meanDev= sum(alpha[1:Ngroups])/Ngroups
  # beta0= alpha0 - meanDev
  # for(j in 1:Ngroups){
  #   beta[j]= alpha[j]-meanDev
  # }
}
'
data=read.csv('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/FruitflyDataReduced.csv')
head(data)
gammaPar= gammaShRaFromModeSD(mode= sd(data$Longevity)/2, sd=2*sd(data$Longevity))
dataList= list(
  y= data$Longevity,
  group= as.numeric(as.factor(data$CompanionNumber)),
  Ntotal= dim(data)[1],
  Ngroups= length(unique(data$CompanionNumber)),
  gammaPar= gammaPar
)

jagMod= jags.model(file=textConnection(modStr),data=dataList,n.chains=3)
cS= coda.samples(model=jagMod,variable.names=c('beta0','beta','sig_y','sig_beta'),n.iter=11000)
DF1= data.frame(cS[[1]])
plot(density(DF1$beta0))
