# 16.2. Caloric restriction extends the lifespan of rats; but does it also increase the
 #variability in lifespan?

#Two-group analysis describing the two groups with student t-distributions with
 #individual shape parameters, but common normality parameter

stanCode= '
data{
  int<lower=1> Ntotal;       //number of rats
  real<lower=0> y[Ntotal]:  //rat lifespan
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

