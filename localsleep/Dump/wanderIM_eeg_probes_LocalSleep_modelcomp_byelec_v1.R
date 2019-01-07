library('lme4')
#library('lmerTest')

t<-read.delim("/Users/tand0009/Data/WanderIM/preproc_eeg/wanderIM_table_localsleep.txt",sep = ",",na.strings = "NaN")
t$SubID<-as.factor(t$SubID)
t$State<-as.factor(t$State)
t$Task<-as.factor(t$Task)

res <- data.frame(Chi=double(63),
                 pVal=double(63),
                 Chan=double(63),
                 MWvsON=double(63),
                 MBvsON=double(63),
                 Task_MWvsON=double(63),
                 Task_MBvsON=double(63),
                 pMWvsON=double(63),
                 pMBvsON=double(63),
                 pTask_MWvsON=double(63),
                 pTask_MBvsON=double(63))
count<-0
for (chan in c(1:63)) {
  print(chan)
  subt<-subset(t,Chan==chan)
  mdl1 <- glmer(nWave ~ nBlock + Task+ (1| SubID), data=subt,family = Gamma(link = "inverse"))
  mdl2 <- glmer(nWave ~ nBlock + State * Task + (1| SubID), data=subt,family = Gamma(link = "inverse"))
  comp <- anova(mdl1,mdl2)
  smdl2<-summary(mdl2)
  
  count<-count+1
  res$Chi[count]<-comp$Chisq[2]
  res$pVal[count]<-comp$`Pr(>Chisq)`[2]
  res$Chan[count]<-chan

  res$MWvsON[count]<-smdl2$coefficients[3,3]
  res$MBvsON[count]<-smdl2$coefficients[4,3]
  res$Task_MWvsON[count]<-smdl2$coefficients[6,3]
  res$Task_MBvsON[count]<-smdl2$coefficients[7,3]
  res$pMWvsON[count]<-smdl2$coefficients[3,4]
  res$pMBvsON[count]<-smdl2$coefficients[4,4]
  res$pTask_MWvsON[count]<-smdl2$coefficients[6,4]
  res$pTask_MBvsON[count]<-smdl2$coefficients[7,4]
}
#res$Cond<-as.factor(res$Cond)
#res$Chan<-as.factor(res$Chan)
write.table (res,file="/Users/tand0009/Data/WanderIM/preproc_eeg/wanderIM_table_localsleep_results.txt",sep = ",",na = "NaN")

res <- data.frame(Chi=double(63),
                  pVal=double(63),
                  Chan=double(63),
                  MWvsON=double(63),
                  MBvsON=double(63),
                  Task_MWvsON=double(63),
                  Task_MBvsON=double(63),
                  pMWvsON=double(63),
                  pMBvsON=double(63),
                  pTask_MWvsON=double(63),
                  pTask_MBvsON=double(63))
count<-0
for (chan in c(1:63)) {
  print(chan)
  subt<-subset(t,Chan==chan)
  mdl1 <- glmer(pWave ~ nBlock + Task+ (1| SubID), data=subt,family = binomial)
  mdl2 <- glmer(pWave ~ nBlock + State * Task + (1| SubID), data=subt,family = binomial)
  comp <- anova(mdl1,mdl2)
  smdl2<-summary(mdl2)
  
  count<-count+1
  res$Chi[count]<-comp$Chisq[2]
  res$pVal[count]<-comp$`Pr(>Chisq)`[2]
  res$Chan[count]<-chan
  
  res$MWvsON[count]<-smdl2$coefficients[3,3]
  res$MBvsON[count]<-smdl2$coefficients[4,3]
  res$Task_MWvsON[count]<-smdl2$coefficients[6,3]
  res$Task_MBvsON[count]<-smdl2$coefficients[7,3]
  res$pMWvsON[count]<-smdl2$coefficients[3,4]
  res$pMBvsON[count]<-smdl2$coefficients[4,4]
  res$pTask_MWvsON[count]<-smdl2$coefficients[6,4]
  res$pTask_MBvsON[count]<-smdl2$coefficients[7,4]
}
#res$Cond<-as.factor(res$Cond)
#res$Chan<-as.factor(res$Chan)
write.table (res,file="/Users/tand0009/Data/WanderIM/preproc_eeg/wanderIM_table_localsleep_results3.txt",sep = ",",na = "NaN")


res <- data.frame(Chi=double(63),
                  pVal=double(63),
                  Chan=double(63),
                  MWvsON=double(63),
                  MBvsON=double(63),
                  Task_MWvsON=double(63),
                  Task_MBvsON=double(63),
                  pMWvsON=double(63),
                  pMBvsON=double(63),
                  pTask_MWvsON=double(63),
                  pTask_MBvsON=double(63))
count<-0
for (chan in c(1:63)) {
  print(chan)
  subt<-subset(t,Chan==chan)
  mdl1 <- lmer(DownSlope ~ nBlock + Task+ (1| SubID), data=subt)
  mdl2 <- lmer(DownSlope ~ nBlock + State * Task + (1| SubID), data=subt)
  comp <- anova(mdl1,mdl2)
  smdl2<-summary(mdl2)
  
  count<-count+1
  res$Chi[count]<-comp$Chisq[2]
  res$pVal[count]<-comp$`Pr(>Chisq)`[2]
  res$Chan[count]<-chan
  
  res$MWvsON[count]<-smdl2$coefficients[3,4]
  res$MBvsON[count]<-smdl2$coefficients[4,4]
  res$Task_MWvsON[count]<-smdl2$coefficients[6,4]
  res$Task_MBvsON[count]<-smdl2$coefficients[7,4]
  res$pMWvsON[count]<-smdl2$coefficients[3,5]
  res$pMBvsON[count]<-smdl2$coefficients[4,5]
  res$pTask_MWvsON[count]<-smdl2$coefficients[6,5]
  res$pTask_MBvsON[count]<-smdl2$coefficients[7,5]
}
#res$Cond<-as.factor(res$Cond)
#res$Chan<-as.factor(res$Chan)
write.table (res,file="/Users/tand0009/Data/WanderIM/preproc_eeg/wanderIM_table_localsleep_results2.txt",sep = ",",na = "NaN")




