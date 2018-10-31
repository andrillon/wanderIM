library('lme4')
library('RWiener')

t<-read.delim("/Users/tand0009/Data/WanderIM/behav/WanderIM_TestResults.txt",sep = ",")
t$SubID<-as.factor(t$SubID)
t$TrCat<-as.factor(t$TrCat)
t$StimID<-as.factor(t$StimID)
t$Task<-as.factor(t$Task)

GO<-subset(t,TrCat=="0")
NOGO<-subset(t,TrCat=="1")

rt<-GO$RT
a<-0.7
v<-3
t0<-0.2

#generate synthetic data
parms <- c(2.00, #boundary separation
           0.30, #non-decision time
           0.50, #initial bias
           0.50) #drift rate
startParms <- c(1, 0.1, 0.1, 1) #boundary, non-decision, initial bias, & drift

tempData <- rwiener(n=100, alpha=parms[1], tau=parms[2], beta=parms[3],
                    delta=parms[4])
#now fit the model to this synthetic data
GO2<-subset(GO,Corr==1)
tempData2<-data.frame(GO$RT,GO$TrCat)
names(tempData2) <- c("q", "resp")
levels(tempData2$resp)<-c("upper","lower")
tempData2$resp<-droplevels(tempData2$resp)
findParms <- optim(startParms, wiener_deviance, dat=tempData2, method="Nelder-Mead")