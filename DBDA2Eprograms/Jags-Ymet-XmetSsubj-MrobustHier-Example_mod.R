# Example for Jags-Ymet-XmetSsubj-MrobustHier.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# Load data file and specity column names of x (predictor) and y (predicted):

# myData = read.csv( file="HierLinRegressData.csv" )
# xName = "X" ; yName = "Y" ; sName="Subj"
# fileNameRoot = "HierLinRegressData-Jags-" 

myData = read.csv( file="IncomeFamszState.csv" )
xName = "Famsz" ; yName = "Income" ; sName="State"
fileNameRoot = "IncomeFamszState-Lin-Jags-"

# myData = read.csv( file="BugsRatsData.csv" )
# xName = "Day" ; yName = "Weight" ; sName="Subj"
# fileNameRoot = "BugsRatsData-Jags-" 

graphFileType = "pdf" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-XmetSsubj-MrobustHier.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
#startTime = proc.time()
mcmcCoda = genMCMC( data=myData , xName=xName , yName=yName , sName=sName ,
                    numSavedSteps=20000 , thinSteps=15 , saveName=fileNameRoot )
#stopTime = proc.time()
#duration = stopTime - startTime
#show(duration)

myData$stateNo= as.numeric(as.factor(myData$State))
chain1= data.frame(mcmcCoda[[1]])

# from each datapoint, simulate posterior predictive distribution for the given state
 #calculate "p-value" of actual income within predictive distribution
set.seed(14082023)

pVal= rep(NA,length(myData))
for (i in 1:dim(myData)[1]){
  stateNo= myData$stateNo[i]
  state_mus= chain1[,stateNo] + chain1[,52+stateNo]*myData$Famsz[i]
  state_devis= rt(n=dim(chain1)[1],df=chain1$nu,ncp=chain1$sigma)
  stateDist= state_mus+state_devis
  
  #probability of actual value being more extreme than observed
  pLower= sum(myData$Income[i] < stateDist)/length(stateDist)
  pHigher= sum(myData$Income[i] > stateDist)/length(stateDist)
  pVal[i]= min(pLower,pHigher)
}
plot(density(pVal))

# #------------------------------------------------------------------------------- 
# # Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in c("beta0mu","beta1mu","nu","sigma","beta0[1]","beta1[1]") ) {
 diagMCMC( codaObject=mcmcCoda , parName=parName , 
           saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , data=myData , xName=xName , yName=yName , sName=sName ,
          compValBeta1=0.0 , ropeBeta1=c(-0.5,0.5) ,
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
