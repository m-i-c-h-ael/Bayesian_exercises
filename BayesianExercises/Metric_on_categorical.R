# Groups of categorical variable with metric measurements (similar to 
 #ANOVA / ANCOVA)

rm(list=ls())
graphics.off()
cat('\014')

library('rjags')
library('tidyverse')
source('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/DBDA2E-utilities.R')

gammaShRaFromModeSD = function( mode , sd ) {
  if ( mode <=0 ) stop("mode must be > 0")
  if ( sd <=0 ) stop("sd must be > 0")
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

modStr= 'model{
  #sd_y= sd(y[1:Ntot])
  #mean_y= mean(y[1:Ntot])

  for(i in 1:Ntot){ y[i] ~ dnorm(a0+a[j[i]], 1/sig_y^2) }
  sig_y ~ dunif(sd_y/100, sd_y*10)
  a0 ~ dnorm(mean_y,1/(5*sd_y)^2)
  for(j in 1:NxLvl) { a[j] ~ dnorm(0,1/sd_y^2) }

  mu= mean(a[1:NxLvl])
  b0= a0 + mu
  b[1:NxLvl]= a[1:NxLvl]-mu
}
'

#alternative model: deviance SD is estimated
modStr2= 'model{
  #sd_y= sd(y[1:Ntot])
  #mean_y= mean(y[1:Ntot])

  for(i in 1:Ntot){ y[i] ~ dnorm(a0+a[j[i]], 1/sig_y^2) }
  sig_y ~ dunif(sd_y/100, sd_y*10)
  a0 ~ dnorm(mean_y,1/(5*sd_y)^2)
  for(j in 1:NxLvl) { a[j] ~ dnorm(0.0, 1/sig_a^2) }
  sig_a ~ dgamma(gamma_sh, gamma_ra)

  mu= mean(a[1:NxLvl])
  b0= a0 + mu
  b[1:NxLvl]= a[1:NxLvl]-mu
}
'

# simulate data
means= c(3,3,3,5)
NperGroup= 7
set.seed(22082023)
y= c( rnorm(n= NperGroup,mean= means[1],sd=1),rnorm(n= NperGroup,mean= means[2],sd=1),
        rnorm(n= NperGroup,mean= means[3],sd=1),rnorm(n= NperGroup,mean= means[4],sd=1) )
dataList= list(
  y= y,
  j= c( rep(1,NperGroup),rep(2,NperGroup),rep(3,NperGroup),rep(4,NperGroup) ),
  Ntot= length(means)*NperGroup,
  NxLvl= length(means),
  mean_y= mean(y),
  sd_y= sd(y),
  gamma_sh= gammaShRaFromModeSD(mode= sd(y)/2, sd= sd(y)*2)[1],
  gamma_ra= gammaShRaFromModeSD(mode= sd(y)/2, sd= sd(y)*2)[2]
)

# visualize
data= data.frame(dataList[c("y","j")])
data$j= as.factor(data$j)
ggplot(data=data,aes(x=j,y=y))+
  geom_jitter(width=.1)+
  coord_cartesian(ylim=c(0,max(data$y)))

# run MCMC
# jagMod= jags.model(file=textConnection(modStr),data=dataList,n.chains=3)
# cS= coda.samples(model=jagMod,variable.names=c('b0','b','sig_y'),n.iter=10000)

jagMod= jags.model(file=textConnection(modStr2),data=dataList,n.chains=3)
cS= coda.samples(model=jagMod,variable.names=c('b0','b','sig_y','sig_a'),n.iter=10000)

DF1= data.frame(cS[[1]])
DF2= data.frame(cS[[2]])
DF3= data.frame(cS[[3]])

for(i in 1:dim(DF1)[2]){
  plot(DF1[,i],type='l',ylab=colnames(DF1)[i],
       main=colnames(DF1)[i])
  lines(DF2[,i],col='blue')
  lines(DF3[,i],col='green')
  
  plot(density(DF1[,i]),main=colnames(DF1)[i])
  lines(density(DF2[,i]),col='blue')
  lines(density(DF3[,i]),col='green')
}

#####################
# Exercise 19.1
rm(list=ls())
library('rjags')
library('tidyverse')
source('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/DBDA2E-utilities.R')
source('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/BayesianExercises/getHDI.R')

dat19.1= read.csv('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/AnovaShrinkageData.csv')
head(dat19.1)
(length(unique(dat19.1$Group)))   #number of groups: A-U

#alternative model: deviance SD is estimated
modStr2= 'model{
  #sd_y= sd(y[1:Ntot])
  #mean_y= mean(y[1:Ntot])

  for(i in 1:Ntot){ y[i] ~ dnorm(a0+a[j[i]], 1/sig_y^2) }
  sig_y ~ dunif(sd_y/100, sd_y*10)
  a0 ~ dnorm(mean_y,1/(5*sd_y)^2)
  for(j in 1:NxLvl) { a[j] ~ dnorm(0.0, 1/sig_a^2) }
  sig_a ~ dgamma(gamma_sh, gamma_ra)

  mu= mean(a[1:NxLvl])
  b0= a0 + mu
  b[1:NxLvl]= a[1:NxLvl]-mu
}
'

#A, Looking at the data
grouped19.1= dat19.1 %>%
  group_by(Group) %>%
  #count()
  summarise(count=length(Y), mean=mean(Y), median=median(Y))   # 4 observations each
  
#B, Is there shrinkage?
j= dat19.1$Group
j_fac= as.factor(j)
j_num= as.numeric(j_fac)
j_levels= levels(j_fac)

dataList19.1= list(
  y= dat19.1$Y,
  j= j_num,  #as numeric
  Ntot= dim(dat19.1)[1],
  NxLvl= length(unique(dat19.1$Group)),
  mean_y= mean(dat19.1$Y),
  sd_y= sd(dat19.1$Y),
  gamma_sh= gammaShRaFromModeSD(mode= sd(dat19.1$Y)/2, sd= sd(dat19.1$Y)*2)[1],
  gamma_ra= gammaShRaFromModeSD(mode= sd(dat19.1$Y)/2, sd= sd(dat19.1$Y)*2)[2]
)


jagMod19.1= jags.model(file=textConnection(modStr2),data=dataList19.1,n.chains=3)
cS19.1= coda.samples(model=jagMod19.1,variable.names = c('b0','b','sig_a','sig_y','a0','a'),n.iter=10000)  #

for(parName in c("sig_y","b0","b[1]","sig_a")){      
  diagMCMC( codaObject= cS19.1 , parName=parName,     
                     saveName=NULL , saveType="jpg" )
}

for(i in 1:length(varnames(cS19.1))){
  diagMCMC( codaObject= cS19.1 , parName=varnames(cS19.1)[i] ,
            saveName=NULL , saveType="jpg" )
}

cS19.1_1= data.frame(cS19.1[[1]])
head(cS19.1_1)
meds19.1= apply(cS19.1_1,2,median)
meds19.1_DFa= cbind.data.frame(j_levels,
                              group_median= meds19.1[22]+meds19.1[1:21],
                              group_min= meds19.1[22]+meds19.1[1:21]-meds19.1[46],
                              group_max= meds19.1[22]+meds19.1[1:21]+meds19.1[46])

ggplot()+
  geom_point(data=dat19.1,aes(x=Group,y=Y))+
  geom_hline(yintercept=meds19.1[22],lty=2,col='green')+
  geom_point(data=meds19.1_DFa,aes(x=j_levels,y=group_median),shape=5,color='red')+
  geom_linerange(data=meds19.1_DFa,aes(x=j_levels,ymin=group_min,ymax=group_max),col='red')+
  ggtitle('Before zero-sum')+
  theme_minimal()

meds19.1_DFb= cbind.data.frame(j_levels,
                               group_median= meds19.1[44]+meds19.1[23:43],
                               group_min= meds19.1[44]+meds19.1[23:43]-meds19.1[46],
                               group_max= meds19.1[44]+meds19.1[23:43]+meds19.1[46])

ggplot()+
  geom_point(data=dat19.1,aes(x=Group,y=Y))+
  geom_hline(yintercept=meds19.1[44],lty=2,col='green')+
  geom_point(data=meds19.1_DFb,aes(x=j_levels,y=group_median),shape=5,color='red')+
  geom_linerange(data=meds19.1_DFb,aes(x=j_levels,ymin=group_min,ymax=group_max),col='red')+
  ggtitle('After zero-sum')+
  theme_minimal()

cS19.1_1_copy= cS19.1_1[,23:46]
colnames(cS19.1_1_copy)[1:21]= j_levels

#calculate contrasts
U_A= cS19.1_1_copy$U - cS19.1_1_copy$A
M_A= cS19.1_1_copy$M - cS19.1_1_copy$A
G_A= cS19.1_1_copy$G - cS19.1_1_copy$A

graphics.off()
HDI_U_A= getHDI(U_A,.95)
hist(U_A,breaks =20,col='skyblue')
lines(x=HDI_U_A,y=c(0,0),lwd=5)

HDI_M_A= getHDI(M_A,.95)
hist(M_A,breaks =20,col='skyblue')
lines(x=HDI_M_A,y=c(0,0),lwd=5)

HDI_G_A= getHDI(G_A,.95)
hist(G_A,breaks =20,col='skyblue')
lines(x=HDI_G_A,y=c(0,0),lwd=5)
 #none is credibly different from 0
 #all three contrasts have a inverted funnel shape

# (Median) estimated deviances vs. (median) group differences 
head(grouped19.1)

(median(U_A))
(grouped19.1$median[grouped19.1$Group=='U'] - grouped19.1$median[grouped19.1$Group=='A'])

(median(M_A))
(grouped19.1$median[grouped19.1$Group=='M'] - grouped19.1$median[grouped19.1$Group=='A'])

(median(G_A))
(grouped19.1$median[grouped19.1$Group=='G'] - grouped19.1$median[grouped19.1$Group=='A'])


# C: Fix SD of normal distribution of deviances to 10*sd(y)
modStr_19.1c= 'model{
  #sd_y= sd(y[1:Ntot])
  #mean_y= mean(y[1:Ntot])

  for(i in 1:Ntot){ y[i] ~ dnorm(a0+a[j[i]], 1/sig_y^2) }
  sig_y ~ dunif(sd_y/100, sd_y*10)
  a0 ~ dnorm(mean_y,1/(5*sd_y)^2)
  for(j in 1:NxLvl) { a[j] ~ dnorm(0,1/(10*sd_y)^2) }

  mu= mean(a[1:NxLvl])
  b0= a0 + mu
  b[1:NxLvl]= a[1:NxLvl]-mu
}
'

jagMod19.1c= jags.model(file=textConnection(modStr_19.1c),data=dataList19.1,n.chains=3)
cS19.1c= coda.samples(model=jagMod19.1c,variable.names = c('b0','b','sig_y'),n.iter=10000)  #

for(parName in c("sig_y","b0","b[1]")){      
  diagMCMC( codaObject= cS19.1c , parName=parName,     
            saveName=NULL , saveType="jpg" )
}

cS19.1c_1= data.frame(cS19.1c[[1]])
head(cS19.1c_1)
meds19.1c= apply(cS19.1c_1,2,median)
meds19.1c_DF= cbind.data.frame(j_levels,
                               group_median= meds19.1c[22]+meds19.1c[1:21],
                               group_min= meds19.1c[22]+meds19.1c[1:21]-meds19.1c[23],
                               group_max= meds19.1c[22]+meds19.1c[1:21]+meds19.1c[23])

graphics.off()
ggplot()+
  geom_point(data=dat19.1,aes(x=Group,y=Y))+
  geom_hline(yintercept=meds19.1c[22],lty=2,col='green')+
  geom_point(data=meds19.1c_DF,aes(x=j_levels,y=group_median),shape=5,color='red')+
  geom_linerange(data=meds19.1c_DF,aes(x=j_levels,ymin=group_min,ymax=group_max),col='red')+
  theme_minimal()

cS19.1c_1_copy= cS19.1c_1
colnames(cS19.1c_1_copy)[1:21]= j_levels

#calculate contrasts
U_Ac= cS19.1c_1_copy$U - cS19.1c_1_copy$A
M_Ac= cS19.1c_1_copy$M - cS19.1c_1_copy$A
G_Ac= cS19.1c_1_copy$G - cS19.1c_1_copy$A

graphics.off()
HDI_U_Ac= getHDI(U_Ac,.95)
hist(U_Ac,breaks =20,col='skyblue')
lines(x=HDI_U_Ac,y=c(0,0),lwd=5)

HDI_M_Ac= getHDI(M_Ac,.95)
hist(M_Ac,breaks =20,col='skyblue')
lines(x=HDI_M_Ac,y=c(0,0),lwd=5)

HDI_G_Ac= getHDI(G_Ac,.95)
hist(G_Ac,breaks =20,col='skyblue')
lines(x=HDI_G_Ac,y=c(0,0),lwd=5)
 #contrast U_A excludes 0

# (Median) estimated deviances vs. (median) group differences 
head(grouped19.1)

(median(U_Ac))
(grouped19.1$median[grouped19.1$Group=='U'] - grouped19.1$median[grouped19.1$Group=='A'])

(median(M_Ac))
(grouped19.1$median[grouped19.1$Group=='M'] - grouped19.1$median[grouped19.1$Group=='A'])

(median(G_Ac))
(grouped19.1$median[grouped19.1$Group=='G'] - grouped19.1$median[grouped19.1$Group=='A'])

# D
 #The median of the estimated sig_a was 3.9 in case of the cross-informed deviances, but
 # 106 in case of fixed sd on prior of deviances
 #The shrinkage is high, probably because there are relatively many groups with relatively
 # few observations and the variation within each group is relatively large

#################
# 19.4. Stan model
rm(list=ls())
graphics.off()
cat('\014')

setwd('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/BayesianExercises')
library('rstan')
source('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/DBDA2E-utilities.R')


modStr= 
'  data{
    real gammaShRa[2];
    int<lower=1> Ntotal;
    int<lower=1> Ngroup;
    real y[Ntotal];
    int<lower=1> x[Ntotal];
  }
  transformed data{
    real<lower=0> sd_y;
    real mean_y;
    sd_y= sd(y);
    mean_y= mean(y);
  }
  parameters{
    real a0;
    real a[Ngroup];
    real<lower=0> y_sigma;
    real<lower=0> b_sigma;
  }
  transformed parameters{
    real mu;
    real b0;
    real b[Ngroup];
    mu= mean(a);
    b0= a0-mu;
    for(j in 1:Ngroup){ b[j]= a[j]+mu; }
  }
  model{
    b_sigma ~ gamma(gammaShRa[1],gammaShRa[2]);
    for(j in 1:Ngroup){
      a[j] ~ normal(0,b_sigma);
    }
    a0 ~ normal(mean_y,sd_y*5);
    y_sigma ~ uniform(sd_y/100,sd_y*10);
    for(i in 1:Ntotal){
      y[i] ~ normal( a0+a[x[i]], y_sigma);
    }
  }
'

data= read.csv('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/FruitflyDataReduced.csv')
head(data)
gammaShRa= gammaShRaFromModeSD(mode= 1/2*sd(data$Longevity),sd=2*sd(data$Longevity))

dataList= list(y=data$Longevity,x=as.numeric(as.factor(data$CompanionNumber)),
               gammaShRa=unlist(gammaShRa), Ntotal= dim(data)[1], 
               Ngroup= length(unique(data$CompanionNumber)))

start= Sys.time()
stMod= stan_model(model_code= modStr)   #call only upon changes above
sp= sampling(object= stMod, data=dataList, chains=3, iter=10000, warmup=500, thin=1)
end= Sys.time()
(end-start)

codaSamp= mcmc.list( lapply( 1:ncol(sp) , 
                   function(x) { mcmc(as.array(sp)[,x,]) } ) )
for(i in 7:dim(sp1)[2]){
  diagMCMC( codaObject= codaSamp , parName=varnames(codaSamp)[i])
}

sp1= as.array(sp)[,1,]
for(i in 7:dim(sp1)[2]){
  hist(sp1[,i],col='skyblue',main=paste(colnames(sp1)[i]),breaks=50)
}

# 10000 iterations take about 14s; ESS on b_sigma ~7500
 #rjags takes a bit longer for the same ESS, but the difference is not huge