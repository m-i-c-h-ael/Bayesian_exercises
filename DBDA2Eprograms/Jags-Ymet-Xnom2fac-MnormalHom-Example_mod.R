# Example for Jags-Ymet-Xnom2fac-MnormalHom.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
#Load The data file 

fileNameRoot = "SeaweedNormalHom-" 
graphFileType = "eps" 
myDataFrame = read.csv( file="SeaweedData.csv" )
# Re-label and re-order the Grazer factor:
myDataFrame$Grazer = factor( myDataFrame$Grazer , 
                          levels=c("None","f","fF","L","Lf","LfF") , 
                          ordered=TRUE , 
                          labels=c("None","f","fF","L","Lf","LfF") )
# Specify the column names in the data file relevant to the analysis:
yName="SeaweedAmt" 
# x1 should be factor with fewer levels, to plot in single pane:
x1Name="Zone" 
x2Name="Grazer" 
# Specify desired contrasts.
# Each main-effect contrast is a list of 2 vectors of level names, 
# a comparison value (typically 0.0), and a ROPE (which could be NULL):
x1contrasts = list( 
  list( c("Full") , c("Assoc") , compVal=0.0 , ROPE=c(-1000,1000) ) ,
  list( c("Assoc") , c("Assis") , compVal=0.0 , ROPE=c(-1000,1000) ) 
)
x2contrasts = list( 
  list( c("CHEM") , c("ENG") , compVal=0.0 , ROPE=c(-1000,1000) ) ,
  list( c("CHEM") , c("PSY") , compVal=0.0 , ROPE=c(-1000,1000) ) ,
  list( c("BFIN") , c("PSY","CHEM","ENG") , compVal=0.0 , ROPE=c(-1000,1000) ) 
)
# Each interaction contrast is a list of 2 lists of 2 vectors of level names, 
# a comparison value (typically 0.0), and a ROPE (which could be NULL)::
x1x2contrasts = list( 
  list( list( c("Full") , c("Assis") ) ,
        list( c("CHEM") , c("ENG") ) ,
        compVal=0.0 , ROPE=c(-1000,1000) ) ,
  list( list( c("Full") , c("Assis") ) ,
        list( c("CHEM") , c("PSY") ) ,
        compVal=0.0 , ROPE=c(-1000,1000) ) ,
  list( list( c("Full") , c("Assoc","Assis") ) ,
        list( c("BFIN") , c("PSY","CHEM","ENG") ) , 
        compVal=0.0 , ROPE=c(-1000,1000) )
) 

# fileNameRoot = "SplitPlotAgriData-NOSUBJ-" 
# graphFileType = "eps" 
# myDataFrame = read.csv( "SplitPlotAgriData.csv" )
# # Specify the column names in the data file relevant to the analysis:
# yName="Yield" 
# x1Name="Fert" 
# x2Name="Till" 
# #xSubjectName="Field"
# x1contrasts = list( 
#   list( c("Deep","Surface") , c("Broad") , compVal=0.0 , ROPE=c(-5,5) ) 
# )
# x2contrasts = list( 
#   list( c("Moldbrd") , c("Ridge") , compVal=0.0 , ROPE=c(-5,5) ) ,
#   list( c("Moldbrd","Ridge") , c("Chisel") , compVal=0.0 , ROPE=c(-5,5) ) 
# )
# x1x2contrasts = list(
#   list( list(  c("Broad") , c("Deep","Surface") ) ,
#         list( c("Chisel","Moldbrd") , c("Ridge") ) ,
#         compVal=0.0 , ROPE=c(-5,5) ) 
# )

#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-Xnom2fac-MnormalHom.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , 
                    yName=yName , x1Name=x1Name , x2Name=x2Name ,
                    numSavedSteps=15000 , thinSteps=5 , 
                    saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in c("b0","b1[1]","b2[1]","b1b2[1,1]","ySigma") ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}




mcmcList= coda::as.mcmc.list(mcmcCoda)
summaryInfo = smryMCMC( mcmcList)

DF1= data.frame(mcmcCoda[[1]])
head(DF1)

data=myDataFrame
data$ZoneNum= as.numeric(as.factor(data$Zone))
data$GrazerNum= as.numeric(as.factor(data$Grazer))
Zones= data$Zone[match(1:max(data$ZoneNum), data$ZoneNum)]
Grazers= as.character(data$Grazer[match(1:max(data$GrazerNum), data$GrazerNum)])

#DF1a= data.frame( matrix(apply(DF1,2,mean),nrow=1 ))
means= apply(DF1,2,mean)
names(means)= c('basel',
                Zones,
                Grazers,
                paste( rep(Zones,each=length(Grazers)),rep(Grazers,times=length(Zones)),sep='x' ),
                paste('m', rep(Zones,each=length(Grazers)),rep(Grazers,times=length(Zones)),sep='x' ),
                'sig_y','sig_b1','sig_b1b2','sig_b2')
#melt(DF1a,measure.vars=2:9,value.name='zoneMean')

# ZoneDF= DF1a[ DF1a$names %in% Zones, ]
# GrazerDF= DF1a[ DF1a$names %in% Grazers, ]
# ZoneVec= as.numeric(ZoneDF$means,names=ZoneDF$names)

# gridName= c( rep(names(means)[2:9],times=length(58:63)),
#               rep(names(means)[58:63],each=length(2:9))
# )
gridNames= expand.grid(Zone=names(means)[2:9],Grazer=names(means)[10:15])
grid= expand.grid(SeawZone= means[2:9], SeawGrazers= means[10:15])
grid= cbind.data.frame(gridNames,basel=rep(means['basel'],dim(grid)[1]),grid)

ggplot(data=grid,aes(group=Zone))+   #
  geom_point(data=data,aes(x=Grazer,y=SeaweedAmt,group=Zone))+
  geom_point(data=grid,aes(x=Grazer,y=basel+SeawGrazers+SeawZone),col='blue')+
  geom_line(data=grid,aes(x=Grazer,y=basel,col='red'))+
  geom_line(data=grid,aes(x=Grazer,y=basel+SeawZone),col='blue',lty=2)+
  theme(legend.position = 'none')+
  facet_wrap(~Zone)




#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , 
                        datFrm=myDataFrame , x1Name=x1Name , x2Name=x2Name ,
                        x1contrasts=x1contrasts , 
                        x2contrasts=x2contrasts , 
                        x1x2contrasts=x1x2contrasts ,
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , 
          datFrm=myDataFrame , yName=yName , x1Name=x1Name , x2Name=x2Name ,
          x1contrasts=x1contrasts , 
          x2contrasts=x2contrasts , 
          x1x2contrasts=x1x2contrasts ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
# Other specific comparisons of cells:
if ( fileNameRoot == "SalaryNormalHom-" ) {
  # THIS x1level minus THAT x1level at AT x2level:
  THISx1 = "Full"
  THATx1 = "Assis"
  ATx2 = "CHEM"
  THISidx = which(levels(myDataFrame[,x1Name])==THISx1)
  THATidx = which(levels(myDataFrame[,x1Name])==THATx1)
  ATidx   = which(levels(myDataFrame[,x2Name])==ATx2)
  openGraph(height=4,width=4)
  compInfo = plotPost( 
    as.matrix(mcmcCoda)[,paste("m[",THISidx,",",ATidx,"]",sep="")] -
      as.matrix(mcmcCoda)[,paste("m[",THATidx,",",ATidx,"]",sep="")] , 
    main=paste(THISx1,"-",THATx1,"@",ATx2) , 
    xlab=paste("Difference in",yName) , 
    compVal=0 ,ROPE=c(-1000,1000) )
  show(compInfo)
  saveGraph(file=paste(fileNameRoot,THISx1,"-",THATx1,"At",ATx2,sep=""),
            type=graphFileType)
  # THIS x1level minus THAT x1level at AT x2level:
  THISx1 = "Full"
  THATx1 = "Assis"
  ATx2 = "PSY"
  THISidx = which(levels(myDataFrame[,x1Name])==THISx1)
  THATidx = which(levels(myDataFrame[,x1Name])==THATx1)
  ATidx   = which(levels(myDataFrame[,x2Name])==ATx2)
  openGraph(height=4,width=4)
  compInfo = plotPost( 
    as.matrix(mcmcCoda)[,paste("m[",THISidx,",",ATidx,"]",sep="")] -
      as.matrix(mcmcCoda)[,paste("m[",THATidx,",",ATidx,"]",sep="")] , 
    main=paste(THISx1,"-",THATx1,"@",ATx2) , 
    xlab=paste("Difference in",yName) , 
    compVal=0 ,ROPE=c(-1000,1000) )
  show(compInfo)
  saveGraph(file=paste(fileNameRoot,THISx1,"-",THATx1,"At",ATx2,sep=""),
            type=graphFileType)
  # THIS x2level minus THAT x2level at AT x1level:
  THISx2 = "PSY"
  THATx2 = "ENG"
  ATx1 = "Full"
  THISidx = which(levels(myDataFrame[,x2Name])==THISx2)
  THATidx = which(levels(myDataFrame[,x2Name])==THATx2)
  ATidx   = which(levels(myDataFrame[,x1Name])==ATx1)
  openGraph(height=4,width=4)
  compInfo = plotPost( 
    as.matrix(mcmcCoda)[,paste("m[",ATidx,",",THISidx,"]",sep="")] -
      as.matrix(mcmcCoda)[,paste("m[",ATidx,",",THATidx,"]",sep="")] , 
    main=paste(THISx2,"-",THATx2,"@",ATx1) , 
    xlab=paste("Difference in",yName) , 
    compVal=0 ,ROPE=c(-1000,1000) )
  show(compInfo)
  saveGraph(file=paste(fileNameRoot,THISx2,"-",THATx2,"At",ATx1,sep=""),
            type=graphFileType)
}
#------------------------------------------------------------------------------- 
