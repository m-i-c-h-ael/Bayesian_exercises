#17.1 Compare models for family income as a function of family size via linear model
 #w/o and with quadratic term
  #based on median income from different states

rm(list=ls())
cat('\014')
graphics.off()

library('rjags')
library('tidyverse')

modStr=' 
data{
  #transform data
  mean_x= mean(x)
  sd_x= sd(x)
  mean_y= mean(y)
  sd_y= sd(y)
  for(i in 1:Ntotal){
    zx[i]= (x[i]-mean_x)/sd_x
    zy[i]= (y[i]-mean_y)/sd_y
  }
}

model{
  for(i in 1:Ntotal){
    zy[i] ~ dt(mu[i],1/zsig^2,nu)   #common sig and nu for all states
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
  zsig ~ dunif(L,H)
  
  M= 0
  S= 10
  L= 1/1000
  H= 1000
  
  #Retransform
  beta0= zbeta0*sd_y + mean_y - zbeta1*mean_x*sd_y/sd_x
  beta1= zbeta1*sd_y/sd_x
  sig= zsig*sd_y
}'

csv= read.csv('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/IncomeFamszState.csv')
csv$s= as.numeric(as.factor(csv$State))
head(csv)
dataList= list(
  y= csv$Income,
  x= csv$Famsz,
  s= csv$s,
  Ntotal= dim(csv)[1],
  Nstates= length(unique(csv$State))
)

jagMod= jags.model(file= textConnection(modStr),data=dataList,n.chains=3)
cS= coda.samples(model= jagMod,variable.names= c('beta0','beta1','sig','nu'),n.iter=10000)

chain1= data.frame(cS[[1]])
chain2= data.frame(cS[[2]])
chain3= data.frame(cS[[3]])
head(chain1)

#plot state 20
state20= levels(as.factor(csv$State))[20]
set.seed(11082023)
rand_idx= sample(1:dim(chain1)[1],size=50)
ggplot(data= csv[csv$State==state20,],aes(x=Famsz,y=Income))+
  geom_point()+
  coord_cartesian(ylim=c(0,max(csv$Income[csv$State==state20])))+
  geom_abline(intercept=chain1$beta0.20.[rand_idx],slope=chain1$beta1.20.[rand_idx],col='blue')

#generate posterior predictive of credible values (variation not considered)
chain1_50= chain1[rand_idx,]
head(csv)
postPredMtx= matrix(NA,nrow=dim(csv)[1],ncol=length(rand_idx))  #50 simulations for each datapoint
for(i in 1:dim(csv)[1]){
  postPredMtx[i,]= chain1_50[,csv$s[i]] + chain1_50[,csv$s[i]+52] * csv$Famsz[i]
}
head(postPredMtx)

#calculate mean difference of each datapoint from prediction (absolute difference!)
csv$meanDiff= rep(NA,dim(csv)[1])
for(i in 1:dim(csv)[1]){
  csv$meanDiff[i]= mean(csv$Income[i]-postPredMtx[i,])
}

ggplot(data=csv,aes(x=State,y=meanDiff,group=as.factor(Famsz),fill=as.factor(Famsz)))+
  geom_bar(stat='identity',position='dodge')

#the data show that the true value is below the predicted value for small and big families
 #the true value is above predicted value for intermediate-size families

####### ///////////////////////////////////// ###########