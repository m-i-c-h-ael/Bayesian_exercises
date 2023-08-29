# Metric response, two predictors

rm(list=ls())
graphics.off()
cat('\014')

setwd('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms')
set.seed(29082023)

library('rjags')
library('reshape2')
library('tidyverse')
source('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/DBDA2E-utilities.R')
source('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/Jags-Ymet-Xnom2fac-MnormalHom.R')

modStr=
'data{
  mean_y= mean(y)
  sd_y= sd(y)
}
model{
  for(i in 1:Ntotal){
    y[i] ~ dnorm(mu[i], 1/sig_y^2)
    mu[i]= a0 + a1[x1[i]] + a2[x2[i]] + a1a2[x1[i],x2[i]]
  }
  sig_y ~ dunif(sd_y/100,sd_y*10)
  a0 ~ dnorm(mean_y,1/(sd_y*5)^2)
  for(i in 1:x1Lvls){
    a1[i] ~ dnorm(0.0, 1/sig_b1^2)
  }
  for (j in 1:x2Lvls){
    a2[j] ~ dnorm(0, 1/sig_b2^2)
    for(i in 1:x1Lvls){
      a1a2[i,j] ~ dnorm(0.0, 1/sig_b1b2^2)
    }
  }
  sig_b1 ~ dgamma(sh[1],sh[2])
  sig_b2 ~ dgamma(sh[1],sh[2])
  sig_b1b2 ~ dgamma(sh[1],sh[2])
  
  # 0-sum game
  for(i in 1:x1Lvls){
    for(j in 1:x2Lvls){
      m[i,j]= a0 + a1[i] + a2[j] + a1a2[i,j]
    }
  }
  b0= mean(m[1:x1Lvls,1:x2Lvls])
  for(i in 1:x1Lvls){
    b1[i]= mean(m[i,x2Lvls])-b0
  }
  for(j in 1:x2Lvls){
    b2[j]= mean(m[1:x1Lvls,j])-b0
    for(i in 1:x1Lvls){
      b1b2[i,j]= m[i,j]-(b0+b1[i]+b2[j])
    }
  }
}'

modStr_datIn=    #mean and sd in input instead of calculated in data-block
'model{
  for(i in 1:Ntotal){
    y[i] ~ dnorm(mu[i], 1/sig_y^2)
    mu[i]= a0 + a1[x1[i]] + a2[x2[i]] + a1a2[x1[i],x2[i]]
  }
  sig_y ~ dunif(sd_y/100,sd_y*10)
  a0 ~ dnorm(mean_y,1/(sd_y*5)^2)
  for(i in 1:x1Lvls){
    a1[i] ~ dnorm(0.0, 1/sig_b1^2)
  }
  for (j in 1:x2Lvls){
    a2[j] ~ dnorm(0, 1/sig_b2^2)
    for(i in 1:x1Lvls){
      a1a2[i,j] ~ dnorm(0.0, 1/sig_b1b2^2)
    }
  }
  sig_b1 ~ dgamma(sh[1],sh[2])
  sig_b2 ~ dgamma(sh[1],sh[2])
  sig_b1b2 ~ dgamma(sh[1],sh[2])
  
  # 0-sum game
  for(i in 1:x1Lvls){
    for(j in 1:x2Lvls){
      m[i,j]= a0 + a1[i] + a2[j] + a1a2[i,j]
    }
  }
  b0= mean(m[1:x1Lvls,1:x2Lvls])
  for(i in 1:x1Lvls){
    b1[i]= mean(m[i,x2Lvls])-b0
  }
  for(j in 1:x2Lvls){
    b2[j]= mean(m[1:x1Lvls,j])-b0
    for(i in 1:x1Lvls){
      b1b2[i,j]= m[i,j]-(b0+b1[i]+b2[j])
    }
  }
}'

modelstring = "  #from book
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dnorm( mu[i] , 1/ySigma^2 )
      mu[i] <- a0 + a1[x1[i]] + a2[x2[i]] + a1a2[x1[i],x2[i]]
    }
    ySigma ~ dunif( ySD/100 , ySD*10 )
    a0 ~ dnorm( yMean , 1/(ySD*5)^2 ) 
    #
    for ( j1 in 1:Nx1Lvl ) { a1[j1] ~ dnorm( 0.0 , 1/a1SD^2 ) }
    a1SD ~ dgamma(agammaShRa[1],agammaShRa[2]) # or try a folded t (Cauchy)
    #
    for ( j2 in 1:Nx2Lvl ) { a2[j2] ~ dnorm( 0.0 , 1/a2SD^2 ) }
    a2SD ~ dgamma(agammaShRa[1],agammaShRa[2]) # or try a folded t (Cauchy)
    #
    for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
      a1a2[j1,j2] ~ dnorm( 0.0 , 1/a1a2SD^2 )
    } }
    a1a2SD ~ dgamma(agammaShRa[1],agammaShRa[2]) # or try a folded t (Cauchy)
    # Convert a0,a1[],a2[],a1a2[,] to sum-to-zero b0,b1[],b2[],b1b2[,] :
    for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
      m[j1,j2] <- a0 + a1[j1] + a2[j2] + a1a2[j1,j2] # cell means 
    } }
    b0 <- mean( m[1:Nx1Lvl,1:Nx2Lvl] )
    for ( j1 in 1:Nx1Lvl ) { b1[j1] <- mean( m[j1,1:Nx2Lvl] ) - b0 }
    for ( j2 in 1:Nx2Lvl ) { b2[j2] <- mean( m[1:Nx1Lvl,j2] ) - b0 }
    for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
      b1b2[j1,j2] <- m[j1,j2] - ( b0 + b1[j1] + b2[j2] )  
    } }
  }
  " # close quote for modelstring

modStr_bookNames=    #my code, but with variable names from book
  'model{
  for(i in 1:Ntotal){
    y[i] ~ dnorm(mu[i], 1/ySigma^2)
    mu[i]= a0 + a1[x1[i]] + a2[x2[i]] + a1a2[x1[i],x2[i]]
  }
  ySigma ~ dunif(ySD/100,ySD*10)
  a0 ~ dnorm(yMean,1/(ySD*5)^2)
  for(i in 1:Nx1Lvl){
    a1[i] ~ dnorm(0.0, 1/a1SD^2)
  }
  for (j in 1:Nx2Lvl){
    a2[j] ~ dnorm(0, 1/a2SD^2)
    for(i in 1:Nx1Lvl){
      a1a2[i,j] ~ dnorm(0.0, 1/a1a2SD^2)
    }
  }
  a1SD ~ dgamma(agammaShRa[1],agammaShRa[2])
  a2SD ~ dgamma(agammaShRa[1],agammaShRa[2])
  a1a2SD ~ dgamma(agammaShRa[1],agammaShRa[2])
  
  # 0-sum game
  for(i in 1:Nx1Lvl){
    for(j in 1:Nx2Lvl){
      m[i,j]= a0 + a1[i] + a2[j] + a1a2[i,j]
    }
  }
  b0= mean(m[1:Nx1Lvl,1:Nx2Lvl])
  for(i in 1:Nx1Lvl){
    b1[i]= mean(m[i,Nx2Lvl])-b0
  }
  for(j in 1:Nx2Lvl){
    b2[j]= mean(m[1:Nx1Lvl,j])-b0
    for(i in 1:Nx1Lvl){
      b1b2[i,j]= m[i,j]-(b0+b1[i]+b2[j])
    }
  }
}'

data= read.csv('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/SeaweedData.csv')
head(data)
data$ZoneNum= as.numeric(as.factor(data$Zone))
data$GrazerNum= as.numeric(as.factor(data$Grazer))

dataList= list(
  y= data$SeaweedAmt,
  x1= data$ZoneNum,
  x2= data$GrazerNum,
  x1Lvls= length(unique(data$ZoneNum)),
  x2Lvls= length(unique(data$GrazerNum)),
  Ntotal= dim(data)[1],
  sh= unlist( gammaShRaFromModeSD(mode= sd(data$SeaweedAmt)/2, sd= sd(data$SeaweedAmt)*2) )
)

dataList2 = list(   #from book, but names adjusted
  y = data$SeaweedAmt ,
  x1 = data$ZoneNum ,
  x2 = data$GrazerNum ,
  Ntotal = dim(data)[1] ,
  Nx1Lvl = length(unique(data$ZoneNum)) ,
  Nx2Lvl = length(unique(data$GrazerNum)) ,
  # data properties for scaling the prior:
  yMean = mean(data$SeaweedAmt) ,
  ySD = sd(data$SeaweedAmt) ,
  agammaShRa = unlist( gammaShRaFromModeSD(mode= sd(data$SeaweedAmt)/2, sd= sd(data$SeaweedAmt)*2) ) 
)

dataList_datIn= list(
  y= data$SeaweedAmt,
  x1= data$ZoneNum,
  x2= data$GrazerNum,
  x1Lvls= length(unique(data$ZoneNum)),
  x2Lvls= length(unique(data$GrazerNum)),
  Ntotal= dim(data)[1],
  mean_y = mean(data$SeaweedAmt) ,
  sd_y = sd(data$SeaweedAmt) ,
  sh= unlist( gammaShRaFromModeSD(mode= sd(data$SeaweedAmt)/2, sd= sd(data$SeaweedAmt)*2) )
)

# jagMod= jags.model(file=textConnection(modStr),data=dataList,n.chains=3)   #my model
# cS= coda.samples(model=jagMod, variable.names=c('b0','b1','b2','b1b2','sig_y',
#                                                 'sig_b1','sig_b2','sig_b1b2'),n.iter=10000)

# jagMod= jags.model(file=textConnection(modelstring),data=dataList2,n.chains=3)   #script from book
# cS= coda.samples(model=jagMod, variable.names=c('b0','b1','b2','b1b2','ySigma',
#           'a1SD','a2SD','a1a2SD'),n.iter=10000)

# jagMod= jags.model(file=textConnection(modStr_datIn),data=dataList_datIn,n.chains=3)   #my model with mean and SD in input
# cS= coda.samples(model=jagMod, variable.names=c('b0','b1','b2','b1b2','sig_y',
#                                                 'sig_b1','sig_b2','sig_b1b2'),n.iter=10000)

jagMod= jags.model(file=textConnection(modStr_bookNames),data=dataList2,n.chains=3)   #my model, but book names
cS= coda.samples(model=jagMod, variable.names=c('b0','b1','b2','b1b2','ySigma',
                                                'a1SD','a2SD','a1a2SD'),n.iter=10000)

#parNames= colnames(data.frame(cS[[1]]))
# for(i in c(1:12,57:67)){
#   diagMCMC(cS, parName = varnames(cS)[i])
# }

mcmcList= coda::as.mcmc.list(cS)
summaryInfo = smryMCMC( mcmcList)

Zones= data$Zone[match(1:max(data$ZoneNum), data$ZoneNum)]
Grazers= data$Grazer[match(1:max(data$GrazerNum), data$GrazerNum)]

colnames(cS[[1]]) [match( paste('b2[',1:length(Grazers),']',sep=''), colnames(cS[[1]]))]= Grazers
colnames(cS[[1]]) [match( paste('b1[',1:length(Zones),']',sep=''), colnames(cS[[1]]))]= Zones
colnames(cS[[1]]) [match(
    paste('b1b2[',
      rep(1:length(Zones),times=length(Grazers)),',',
      rep(1:length(Grazers),each=length(Zones)),
    ']',sep=''),
  colnames(cS[[1]]))]= 
  paste( rep(Zones,each=length(Grazers)),rep(Grazers,times=length(Zones)),sep='x' )
colnames(cS[[1]])

means= apply(DF1,2,mean)
names(means)= c('basel',
                  Zones,
                  paste( rep(Zones,each=length(Grazers)),rep(Grazers,times=length(Zones)),sep='x' ),
                  Grazers,
                  'sig_b1','sig_b1b2','sig_b2','sig_y')
#melt(DF1a,measure.vars=2:9,value.name='zoneMean')

# ZoneDF= DF1a[ DF1a$names %in% Zones, ]
# GrazerDF= DF1a[ DF1a$names %in% Grazers, ]
# ZoneVec= as.numeric(ZoneDF$means,names=ZoneDF$names)

# gridName= c( rep(names(means)[2:9],times=length(58:63)),
#               rep(names(means)[58:63],each=length(2:9))
# )
gridNames= expand.grid(Zone=names(means)[2:9],Grazer=names(means)[58:63])
grid= expand.grid(SeawZone= means[2:9], SeawGrazers= means[58:63])
grid= cbind.data.frame(gridNames,basel=rep(means['basel'],dim(grid)[1]),grid)

#rearrange factor levels
levels(grid$Grazer)
#ord= match(levels(grid$Grazer),c("None","f","fF","L","Lf","LfF"))
grid$Grazer= factor(grid$Grazer,levels=c("None","f","fF","L","Lf","LfF"))
data$Grazer= factor(data$Grazer,levels=c("None","f","fF","L","Lf","LfF"))

ggplot(data=grid,aes(group=Zone))+   #
  geom_point(data=data,aes(x=Grazer,y=SeaweedAmt,group=Zone))+
  geom_point(data=grid,aes(x=Grazer,y=basel+SeawGrazers+SeawZone),col='blue')+
  geom_line(data=grid,aes(x=Grazer,y=basel,col='red'))+
  geom_line(data=grid,aes(x=Grazer,y=basel+SeawZone),col='blue',lty=2)+
  theme(legend.position = 'none')+
  facet_wrap(~Zone)

plot(density(DF1$b0))

#B, Effect of small fish
DF1a= DF1
colnames(DF1a)= c('basel',
  Zones,
  paste( rep(Zones,each=length(Grazers)),rep(Grazers,times=length(Zones)),sep='x' ),
  Grazers,
  'sig_b1','sig_b1b2','sig_b2','sig_y')
head(DF1a)
 
 #small fish vs. none 
diff_f.None= DF1a$f - DF1a$None
plot(density(diff_f.None))

 #limpets+small fish vs. limpets
diff_Lf.L= DF1a$Lf-DF1a$L
plot(density(diff_Lf.L))

 #avg(f+Lf) vs. avg(none+limpets)
diff_fLf.noneL= (DF1a$f+DF1a$Lf)/2 - (DF1a$None+DF1a$L)/2
plot(density(diff_fLf.noneL))

#C, Effect of limpets
diffLimp= (DF1a$L+DF1a$Lf+DF1a$LfF)/3 - (DF1a$None+DF1a$f+DF1a$fF)/3
plot(density(diffLimp))

#D, Difference zone A vs. D
diff_D_A= DF1a$D - DF1a$A
plot(density(diff_D_A))
 
#E, Limpet effect in A vs. D
Limp_DvsA= (DF1a$DxL+DF1a$DxLf+DF1a$DxLfF)/3 - (DF1a$AxL+DF1a$AxLf+DF1a$AxLfF)/3
plot(density(Limp_DvsA))

#F: Heterogenous variance model: Variances around group means are drawn from gamma 
 # distribution

modStr_hetVar=
  'data{
  mean_y= mean(y)
  sd_y= sd(y)
}
model{
  for(i in 1:Ntotal){
    y[i] ~ dnorm(mu[i], 1/sig_y[x1[i],x2[i]]^2)
    mu[i]= a0 + a1[x1[i]] + a2[x2[i]] + a1a2[x1[i],x2[i]]
  }
  a0 ~ dnorm(mean_y,1/(sd_y*5)^2)
  for(i in 1:x1Lvls){
    a1[i] ~ dnorm(0.0, 1/sig_b1^2)
  }
  for (j in 1:x2Lvls){
    a2[j] ~ dnorm(0, 1/sig_b2^2)
    for(i in 1:x1Lvls){
      a1a2[i,j] ~ dnorm(0.0, 1/sig_b1b2^2)
      sig_y[i,j] ~ dgamma(sh[1],sh[2])
    }
  }
  sig_b1 ~ dgamma(sh[1],sh[2])
  sig_b2 ~ dgamma(sh[1],sh[2])
  sig_b1b2 ~ dgamma(sh[1],sh[2])
  
  # 0-sum game
  for(i in 1:x1Lvls){
    for(j in 1:x2Lvls){
      m[i,j]= a0 + a1[i] + a2[j] + a1a2[i,j]
    }
  }
  b0= mean(m[1:x1Lvls,1:x2Lvls])
  for(i in 1:x1Lvls){
    b1[i]= mean(m[i,x2Lvls])-b0
  }
  for(j in 1:x2Lvls){
    b2[j]= mean(m[1:x1Lvls,j])-b0
    for(i in 1:x1Lvls){
      b1b2[i,j]= m[i,j]-(b0+b1[i]+b2[j])
    }
  }
}'

jagMod_hetVar= jags.model(file=textConnection(modStr_hetVar),data=dataList,n.chains=3)
cS_hetVar= coda.samples(model=jagMod_hetVar, variable.names=c('b0','b1','b2','b1b2','sig_y',
                                                'sig_b1','sig_b2','sig_b1b2'),n.iter=10000)

mcmcList_hetVar= coda::as.mcmc.list(cS_hetVar)
summaryInfo_hetVar = smryMCMC( mcmcList_hetVar)

DF1_hetVar= data.frame(cS_hetVar[[1]])
head(DF1_hetVar)


#DF1a= data.frame( matrix(apply(DF1,2,mean),nrow=1 ))
means_hetVar= apply(DF1,2,mean)
replaceNames= c('basel',
                Zones,
                paste( rep(Zones,each=length(Grazers)),rep(Grazers,times=length(Zones)),sep='x' ),
                Grazers,
                'sig_b1','sig_b1b2','sig_b2')
names(means_hetVar)[1:length(replaceNames)]= replaceNames

gridNames_hetVar= expand.grid(Zone=names(means_hetVar)[2:9],Grazer=names(means_hetVar)[58:63])
grid_hetVar= expand.grid(SeawZone= means_hetVar[2:9], SeawGrazers= means_hetVar[58:63])
grid_hetVar= cbind.data.frame(gridNames_hetVar,basel=rep(means_hetVar['basel'],dim(grid_hetVar)[1]),grid_hetVar)

ggplot(data=grid_hetVar,aes(group=Zone))+   #
  geom_point(data=data,aes(x=Grazer,y=SeaweedAmt,group=Zone))+
  geom_point(data=grid_hetVar,aes(x=Grazer,y=basel+SeawGrazers+SeawZone),col='blue')+
  geom_line(data=grid_hetVar,aes(x=Grazer,y=basel,col='red'))+
  geom_line(data=grid_hetVar,aes(x=Grazer,y=basel+SeawZone),col='blue',lty=2)+
  theme(legend.position = 'none')+
  facet_wrap(~Zone)
