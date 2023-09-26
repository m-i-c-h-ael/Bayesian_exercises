# 16.2. Caloric restriction extends the lifespan of rats; but does it also increase the
 #variability in lifespan?

#Two-group analysis describing the two groups with student t-distributions with
 #individual shape parameters, but common normality parameter

rm(list=ls())
cat('\014')
graphics.off()

library('rstan')
source('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/BayesianExercises/getHDI.R')

#specify model
stanCode= '
data {
  int<lower=1> Ntotal;       //number of rats
  real<lower=0> y[Ntotal];  //rat lifespan
  int x[Ntotal];            //category: ad libidum or CR
  real<lower=0> y_mean;
  real<lower=0> y_sd;
}
transformed data{
  real<lower=0> sd_norm;
  real<lower=0> unif_lo;
  real<lower=1> unif_hi;
  sd_norm= y_sd*100;       //sd for prior on mu  
  unif_lo= y_sd/1000;
  unif_hi= y_sd*1000;
}
parameters{
  real<lower=0> nuMin1;
  real<lower=0> mu[2];    //mean for both groups
  real<lower=0> sigma[2];
}
transformed parameters{
  real<lower=0> lambda;
  real<lower=1> nu;
  lambda= 1/29.0;
  nu= nuMin1+1;
}
model{
  nuMin1 ~ exponential(lambda);
  mu ~ normal(y_mean, sd_norm);
  sigma ~ uniform(unif_lo, unif_hi);
  for (i in 1:Ntotal) {
    y[i] ~ student_t(nu, mu[x[i]], sigma[x[i]]);
  }
}
'
# prepare data
data= read.csv('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/RatLives.csv')
head(data)
dataList= list(
  x= as.numeric(as.factor(data$Group)),
  y= data$DaysLive,
  Ntotal= dim(data)[1],
  y_mean= mean(data$DaysLive),
  y_sd= sd(data$DaysLive)
)

#run MCMC
if(! exists(DSO)){  #creation of DSO takes long; only create if not yet
  DSO= stan_model(model_code= stanCode)
}
stanOUT= sampling(DSO, data=dataList, chains=3, iter=1000, warmup=500)

# isolate chains
chain1= as.array(stanOUT)[,1,]
chain2= as.array(stanOUT)[,2,]
chain3= as.array(stanOUT)[,3,]
head(chain1)

#check runs
for(i in 1:dim(chain1)[2]){
  plot(chain1[,i],main=paste(colnames(chain1)[i]),type='l')
  lines(chain2[,i],col='red')
  lines(chain3[,i],col='green')
}

#visualize posterior sample
for(i in 1:dim(chain1)[2]){
  #calculate HDI
  HDI= getHDI(chain1[,i],HDI_width= 0.95)
  
  plot(density(chain1[,i]),main=paste(colnames(chain1)[i]),type='l')
  lines(density(chain2[,i]),col='red')
  lines(density(chain3[,i]),col='green')
  
  lines(x=HDI,y=c(0,0),lwd=3)  #HDI based on chain1
}

#effect size: (mu2-mu1) /sqrt((sig1^2+sig2^2)/2)
effSize1= (chain1[,3]-chain1[,2])/ sqrt( (chain1[,4]^2+chain1[,5]^2)/2 )
plot(density(effSize1),main="Effect Size")
effSize2= (chain2[,3]-chain2[,2])/ sqrt( (chain2[,4]^2+chain2[,5]^2)/2 )
lines(density(effSize2),col='red')
effSize3= (chain3[,3]-chain3[,2])/ sqrt( (chain3[,4]^2+chain3[,5]^2)/2 )
lines(density(effSize3),col='green')

HDI_effSize1= getHDI(effSize1,HDI_width= 0.95)
lines(x=HDI_effSize1,y=c(0,0),lwd=3)  #HDI based on chain1

ROPE= c(-.1,.1)  #arbitrary ROPE
lines(x=ROPE,y=c(0.01,0.01),col='blue',lwd=3)  #HDI based on chain1

#differences in scale
diffSig1= chain1[,5]-chain1[,4]
plot(density(diffSig1), main="Diff. in sigma: CR vs. ad lib.")
diffSig2= chain2[,5]-chain2[,4]
lines(density(diffSig2), col='red')
diffSig3= chain3[,5]-chain3[,4]
lines(density(diffSig3), col='green')

# calculate posterior predictive values - overlay on original data
par(mfrow=c(2,1))
pred_grp1_ch1= rnorm(dim(chain1)[1],mean=chain1[,2],sd=chain1[,4])
pred_grp1_ch2= rnorm(dim(chain2)[1],mean=chain2[,2],sd=chain2[,4])
pred_grp1_ch3= rnorm(dim(chain3)[1],mean=chain3[,2],sd=chain3[,4])
plot(density(pred_grp1_ch1),main="Post. pred. ad lib.",xlim=c(300,1700))
lines(density(pred_grp1_ch2),col='red')
lines(density(pred_grp1_ch3),col='green')
lines(density(data$DaysLive[data$Group=='Adlib']),col='blue')
#hist(data$DaysLive[data$Group=='Adlib'],freq=FALSE,add=TRUE,alpha=.7,col='skyblue')

pred_grp2_ch1= rnorm(dim(chain1)[1],mean=chain1[,3],sd=chain1[,5])
pred_grp2_ch2= rnorm(dim(chain2)[1],mean=chain2[,3],sd=chain2[,5])
pred_grp2_ch3= rnorm(dim(chain3)[1],mean=chain3[,3],sd=chain3[,5])
plot(density(pred_grp2_ch1),main="Post. restrict",xlim=c(300,1700))
lines(density(pred_grp2_ch2),col='red')
lines(density(pred_grp2_ch3),col='green')
lines(density(data$DaysLive[data$Group=='Restrict']),col='blue')

#The groups are largely different in their central tendencies
#The groups largely differ in their scale parameter (higher variation with CR)
#The normality parameter is only around 3; this suggests sample is heavily influenced
 #by outliers and far from normality

# B
 #does transformation through squaring improve models in presence of left-tailing?
head(data)
data$DaysLiveSq= data$DaysLive^2

dataList_sq= list(
  x= as.numeric(as.factor(data$Group)),
  y= data$DaysLiveSq,
  Ntotal= dim(data)[1],
  y_mean= mean(data$DaysLiveSq),
  y_sd= sd(data$DaysLiveSq)
)

hist(dataList_sq$y)

if (! exists('DSO_sq')){
  DSO_sq= stan_model(model_code= stanCode)
}
stanOUT_sq= sampling(DSO_sq, data=dataList_sq, chains=3, iter=1000, warmup=500,
                     init=list( list(nuMin1=3,mu=c(70000,70000),sigma=c(150000,150000)),
                                list(nuMin1=3,mu=c(750000,750000),sigma=c(200000,200000)),
                                list(nuMin1=3,mu=c(800000,800000),sigma=c(250000,250000)) ) 
                     )

# isolate chains
chain1_sq= as.array(stanOUT_sq)[,1,]
chain2_sq= as.array(stanOUT_sq)[,2,]
chain3_sq= as.array(stanOUT_sq)[,3,]
head(chain1_sq)

#check runs
for(i in 1:dim(chain1_sq)[2]){
  plot(chain1_sq[,i],main=paste(colnames(chain1_sq)[i]),type='l')
  lines(chain2_sq[,i],col='red')
  lines(chain3_sq[,i],col='green')
}

#visualize posterior sample
for(i in 1:dim(chain1_sq)[2]){
  #calculate HDI
  HDI= getHDI(chain1_sq[,i],HDI_width= 0.95)
  
  plot(density(chain1_sq[,i]),main=paste(colnames(chain1_sq)[i]),type='l')
  lines(density(chain2_sq[,i]),col='red')
  lines(density(chain3_sq[,i]),col='green')
  
  lines(x=HDI,y=c(0,0),lwd=3)  #HDI based on chain1
}

#effect size: (mu2-mu1) /sqrt((sig1^2+sig2^2)/2)
effSize1_sq= (chain1_sq[,3]-chain1_sq[,2])/ sqrt( (chain1_sq[,4]^2+chain1_sq[,5]^2)/2 )
plot(density(effSize1_sq),main="Effect Size")
effSize2_sq= (chain2_sq[,3]-chain2_sq[,2])/ sqrt( (chain2_sq[,4]^2+chain2_sq[,5]^2)/2 )
lines(density(effSize2_sq),col='red')
effSize3_sq= (chain3_sq[,3]-chain3_sq[,2])/ sqrt( (chain3_sq[,4]^2+chain3_sq[,5]^2)/2 )
lines(density(effSize3_sq),col='green')

#Groups clearly differ in central tendencies
#Scale parameters clearly differ between groups
#Normality parameter covers relatively wide range, but its mean is higher
#Effect size remains similar to un-transformed data