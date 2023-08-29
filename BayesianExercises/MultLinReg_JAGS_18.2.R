# Multiple linear regression using JAGS

rm(list=ls())
cat('\014')
graphics.off()

library('rjags')

modStr= '
#Transform data
data{
  ymean= mean(y)
  ysd= sd(y)
  for(i in 1:Ntotal){
    zy[i]= (y[i]-ymean)/ysd
  }
  for(j in 1:Npred){
    xmean[j]= mean(x[,j])
    xsd[j]= sd(x[,j])
    for(i in 1:Ntotal){
        zx[i,j]= (x[i,j]-xmean[j])/xsd[j]
    }
  }
}

#Construct model
model{
  for(i in 1:Ntotal){
    zy[i] ~ dt( zbeta0 + sum(zbeta[1:Npred]*zx[i,1:Npred]), 1/zsigma^2,nu )
  }
  
  for(j in 1:Npred){
    zbeta[j] ~ dnorm(0,1/2^2)
  }
  zbeta0 ~ dnorm(0,1/2^2)
  zsigma ~ dunif(1/1000, 1000)
  nuMin1 ~ dexp(1/29.0)
  nu= nuMin1+1
  
  #Retransform
  beta0= zbeta0*ysd + ymean - sum (zbeta[1:Npred] * xmean[1:Npred] / xsd[1:Npred]) * ysd
  for(j in 1:Npred){
    beta[j]= zbeta[j]/xsd[j]*ysd
  }
  sigma= zsigma*ysd
}'

## 18.2.
data= read.csv('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/MultLinRegrPlotUnif.csv')
head(data) 
dataList= list(
  x=data[,1:2],
  y=data[,3],
  Ntotal= dim(data)[1],
  Npred= 2
)

# First analyse with frequentist stats
# freqMod= lm(y~x1+x2,data= data)
# summary(freqMod)

jagMod= jags.model(file= textConnection(modStr),data=dataList,n.chains=3)
cS= coda.samples(model=jagMod,variable.names=c('beta0','beta','sigma','nu'),n.iter=10000)
DF1= data.frame(cS[[1]])
DF2= data.frame(cS[[2]])
DF3= data.frame(cS[[3]])

plot(DF1$beta0,type='l')
lines(DF2$beta0,col='blue')
lines(DF3$beta0,col='green')

plot(DF1$beta.1.,type='l')
lines(DF2$beta.1.,col='blue')
lines(DF3$beta.1.,col='green')

plot(DF1$beta.2.,type='l')
lines(DF2$beta.2.,col='blue')
lines(DF3$beta.2.,col='green')

plot(density(DF1$beta0))
plot(density(DF1$beta.1.))
plot(density(DF1$beta.2.))
plot(density(DF1$sigma))

cor(data$x1,data$x2)

#calculate R-squared
 # get predictions
#SSQ_res= 0
#SSQ_mod= 0
y_pred= matrix(NA,nrow=dim(data)[1],ncol=1000)
for (i in 1:dim(data)[1]){
  y_pred[i,]= rnorm(1000,
                    mean= DF1$beta0[9001:10000] + 
                          DF1$beta.1.[9001:10000]*data$x1[i] + 
                          DF1$beta.2.[9001:10000]*data$x2[i],
                    sd= DF1$sigma[9001:10000]
                    ) #based on 1000 predictions
  #SSQ_res= SSQ_res + sum( (data$y[i] - y_pred)^2 )
  #SSQ_mod= SSQ_mod + sum((y_mean - y_pred)^2 * 1000)
}
#(R_sq= SSQ_mod/(SSQ_res+SSQ_mod) )

y_mean= mean(data$y)
SSQres_vec= apply(y_pred,2,function(x){ sum((x-data$y)^2)})
SSQmod_vec= apply(y_pred,2,function(x){ sum((x-y_mean)^2)})
Rsq_vec= SSQmod_vec/(SSQres_vec+SSQmod_vec)
hist(Rsq_vec,col='skyblue')   #not the same as in the solution..


########
# B,
 # use only predictor x1
head(data)

dataList_x1= list(
  x=data.frame(data[,1]),   #needs to be data.frame, even if only one column
  y=data[,3],
  Ntotal= dim(data)[1],
  Npred= 1
)

jagMod_x1= jags.model(file= textConnection(modStr),data=dataList_x1,n.chains=3)
cS_x1= coda.samples(model=jagMod_x1,variable.names=c('beta0','beta','sigma','nu'),n.iter=10000)
DF1_x1= data.frame(cS_x1[[1]])
DF2_x1= data.frame(cS_x1[[2]])
DF3_x1= data.frame(cS_x1[[3]])

for(i in 1:dim(DF1_x1)[2]){
  plot(DF1_x1[,i],type='l',ylab=colnames(DF1_x1)[i],main=colnames(DF1_x1)[i])
  lines(DF2_x1[,i],col='blue')
  lines(DF3_x1[,i],col='green')
}

# The y values at x1=0 are distributed in a range 10 to 30 (x2 is in range 0 to 10); therefore
 #the intercept is 20 if not explicitly considering x2
# As the variation around the mean now includes variation due to y2, sigma has increased to 6

# C
dim(data)
data_101to150= data[101:150,]
head(data_101to150)
summary(data)
summary(data_101to150) #selection seems random

# Both predictors
dataList_101to150= list(
  x= data_101to150[,1:2],
  y= data_101to150[,3],
  Ntotal= dim(data_101to150)[1],
  Npred= dim(dataList_101to150$x)[2]
)

jagMod_101to150= jags.model(file= textConnection(modStr),data=dataList_101to150,n.chains=3)
cS_101to150= coda.samples(model=jagMod_101to150,variable.names=c('beta0','beta','sigma','nu'),
                          n.iter=10000)

DF1_101to150= data.frame(cS_101to150[[1]])
DF2_101to150= data.frame(cS_101to150[[2]])
DF3_101to150= data.frame(cS_101to150[[3]])

for(i in 1:dim(DF1_101to150)[2]){
  plot(DF1_101to150[,i],type='l',ylab=colnames(DF1_101to150)[i],main=colnames(DF1_101to150)[i])
  lines(DF2_101to150[,i],col='blue')
  lines(DF3_101to150[,i],col='green')
  
  plot(density(DF1_101to150[,i]),main=colnames(DF1_101to150)[i])
  lines(density(DF2_101to150[,i]),col='blue')
  lines(density(DF3_101to150[,i]),col='green')
}
 #The mean of beta1 is now a bit smaller, ~0.9; the menan of beta2 a bit higher, ~2.1

# Only x1 as predictor
dataList_101to150_x1= list(
  x= data.frame(data_101to150[,1]),
  y= data_101to150[,3],
  Ntotal= dim(data_101to150)[1],
  Npred= dim(dataList_101to150_x1$x)[2]
)

jagMod_101to150_x1= jags.model(file= textConnection(modStr),data=dataList_101to150_x1,n.chains=3)
cS_101to150_x1= coda.samples(model=jagMod_101to150_x1,variable.names=c('beta0','beta','sigma','nu'),
                          n.iter=10000)

DF1_101to150_x1= data.frame(cS_101to150_x1[[1]])
DF2_101to150_x1= data.frame(cS_101to150_x1[[2]])
DF3_101to150_x1= data.frame(cS_101to150_x1[[3]])

for(i in 1:dim(DF1_101to150_x1)[2]){
  plot(DF1_101to150_x1[,i],type='l',ylab=colnames(DF1_101to150_x1)[i],
       main=colnames(DF1_101to150_x1)[i])
  lines(DF2_101to150_x1[,i],col='blue')
  lines(DF3_101to150_x1[,i],col='green')
  
  plot(density(DF1_101to150_x1[,i]),main=colnames(DF1_101to150_x1)[i])
  lines(density(DF2_101to150_x1[,i]),col='blue')
  lines(density(DF3_101to150_x1[,i]),col='green')
}

#beta0 is a bit higher at ~22, beta lower at ~0.7
 #When x2 is included, beta1 is interpreted as the rate of change/slope of y in the x1-
 #direction for each value of x2
  #When x2 is not included, beta1 is the overall rate of change in the x1-direction