library('lme4')
#library('lmerTest')

t<-read.delim("/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults_EEG_allCh3.txt",sep = ",",na.strings = "NaN")
t$SubID<-as.factor(t$SubID)
t$State<-as.factor(t$State)
t$Task<-as.factor(t$Task)
t$Look<-as.factor(t$Look)
t$Orig<-as.factor(t$Orig)
t$Freq<-as.factor(t$Freq)

res <- data.frame(Chi=double(63*3),
                 pVal=double(63*3),
                 Chan=double(63*3),
                 Cond=double(63*3))
count<-0
for (chan in c(1:63)) {
  print(chan)
  subt_F<-subset(t,Chan==chan & (Freq=="1" | Freq=="2"))
  mdl1_F <- lmer(SNR ~ nBlock + Task * Freq + (1| SubID), data=subt_F,REML=FALSE)
  mdl2_F <- lmer(SNR ~ nBlock + State * Task * Freq+ (1| SubID), data=subt_F,REML=FALSE)
  comp_F <- anova(mdl1_F,mdl2_F)
  
  count<-count+1
  res$Chi[count]<-comp_F$Chisq[2]
  res$pVal[count]<-comp_F$`Pr(>Chisq)`[2]
  res$Chan[count]<-chan
  res$Cond[count]<-1
  
  subt_2F<-subset(t,Chan==chan & (Freq=="3" | Freq=="4"))
  mdl1_2F <- lmer(SNR ~ nBlock + Task * Freq + (1| SubID), data=subt_2F,REML=FALSE)
  mdl2_2F <- lmer(SNR ~ nBlock + State * Task * Freq+ (1| SubID), data=subt_2F,REML=FALSE)
  comp_2F <- anova(mdl1_2F,mdl2_2F)

  count<-count+1
  res$Chi[count]<-comp_2F$Chisq[2]
  res$pVal[count]<-comp_2F$`Pr(>Chisq)`[2]
  res$Chan[count]<-chan
  res$Cond[count]<-2
  
  
  subt_IM<-subset(t,Chan==chan & (Freq=="5"))
  mdl1_IM <- lmer(SNR ~ nBlock + Task + (1| SubID), data=subt_IM,REML=FALSE)
  mdl2_IM <- lmer(SNR ~ nBlock + State * Task + (1| SubID), data=subt_IM,REML=FALSE)
  comp_IM <- anova(mdl1_IM,mdl2_IM)
  
  count<-count+1
  res$Chi[count]<-comp_IM$Chisq[2]
  res$pVal[count]<-comp_IM$`Pr(>Chisq)`[2]
  res$Chan[count]<-chan
  res$Cond[count]<-3
}
#res$Cond<-as.factor(res$Cond)
#res$Chan<-as.factor(res$Chan)
write.table (res,file="/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults_EEG_allCh3_results.txt",sep = ",",na = "NaN")

resD <- data.frame(Chi=double(63*3),
                  pVal=double(63*3),
                  Chan=double(63*3),
                  Cond=double(63*3))
count<-0
for (chan in c(1:63)) {
  print(chan)
  subt_F<-subset(t,Chan==chan & (Freq=="1" | Freq=="2") & Task="2")
  mdl1_F <- lmer(SNR ~ nBlock + Task * Freq + (1| SubID), data=subt_F,REML=FALSE)
  mdl2_F <- lmer(SNR ~ nBlock + State * Task * Freq+ (1| SubID), data=subt_F,REML=FALSE)
  comp_F <- anova(mdl1_F,mdl2_F)
  
  count<-count+1
  resD$Chi[count]<-comp_F$Chisq[2]
  resD$pVal[count]<-comp_F$`Pr(>Chisq)`[2]
  resD$Chan[count]<-chan
  resD$Cond[count]<-1
  
  subt_2F<-subset(t,Chan==chan & (Freq=="3" | Freq=="4") & Task="2")
  mdl1_2F <- lmer(SNR ~ nBlock + Task * Freq + (1| SubID), data=subt_2F,REML=FALSE)
  mdl2_2F <- lmer(SNR ~ nBlock + State * Task * Freq+ (1| SubID), data=subt_2F,REML=FALSE)
  comp_2F <- anova(mdl1_2F,mdl2_2F)
  
  count<-count+1
  resD$Chi[count]<-comp_2F$Chisq[2]
  resD$pVal[count]<-comp_2F$`Pr(>Chisq)`[2]
  resD$Chan[count]<-chan
  resD$Cond[count]<-2
  
  
  subt_IM<-subset(t,Chan==chan & (Freq=="5") & Task="2")
  mdl1_IM <- lmer(SNR ~ nBlock + Task + (1| SubID), data=subt_IM,REML=FALSE)
  mdl2_IM <- lmer(SNR ~ nBlock + State * Task + (1| SubID), data=subt_IM,REML=FALSE)
  comp_IM <- anova(mdl1_IM,mdl2_IM)
  
  count<-count+1
  resD$Chi[count]<-comp_IM$Chisq[2]
  resD$pVal[count]<-comp_IM$`Pr(>Chisq)`[2]
  resD$Chan[count]<-chan
  resD$Cond[count]<-3
  }
#res$Cond<-as.factor(res$Cond)
#res$Chan<-as.factor(res$Chan)
write.table (resD,file="/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults_EEG_allCh3_results_D.txt",sep = ",",na = "NaN")

resF <- data.frame(Chi=double(63*3),
                  pVal=double(63*3),
                  Chan=double(63*3),
                  Cond=double(63*3))
count<-0
for (chan in c(1:63)) {
  print(chan)
  subt_F<-subset(t,Chan==chan & (Freq=="1" | Freq=="2") & Task="1")
  mdl1_F <- lmer(SNR ~ nBlock + Task * Freq + (1| SubID), data=subt_F,REML=FALSE)
  mdl2_F <- lmer(SNR ~ nBlock + State * Task * Freq+ (1| SubID), data=subt_F,REML=FALSE)
  comp_F <- anova(mdl1_F,mdl2_F)
  
  count<-count+1
  resF$Chi[count]<-comp_F$Chisq[2]
  resF$pVal[count]<-comp_F$`Pr(>Chisq)`[2]
  resF$Chan[count]<-chan
  resF$Cond[count]<-1
  
  subt_2F<-subset(t,Chan==chan & (Freq=="3" | Freq=="4") & Task="1")
  mdl1_2F <- lmer(SNR ~ nBlock + Task * Freq + (1| SubID), data=subt_2F,REML=FALSE)
  mdl2_2F <- lmer(SNR ~ nBlock + State * Task * Freq+ (1| SubID), data=subt_2F,REML=FALSE)
  comp_2F <- anova(mdl1_2F,mdl2_2F)
  
  count<-count+1
  resF$Chi[count]<-comp_2F$Chisq[2]
  resF$pVal[count]<-comp_2F$`Pr(>Chisq)`[2]
  resF$Chan[count]<-chan
  resF$Cond[count]<-2
  
  
  subt_IM<-subset(t,Chan==chan & (Freq=="5") & Task="1")
  mdl1_IM <- lmer(SNR ~ nBlock + Task + (1| SubID), data=subt_IM,REML=FALSE)
  mdl2_IM <- lmer(SNR ~ nBlock + State * Task + (1| SubID), data=subt_IM,REML=FALSE)
  comp_IM <- anova(mdl1_IM,mdl2_IM)
  
  count<-count+1
  resF$Chi[count]<-comp_IM$Chisq[2]
  resF$pVal[count]<-comp_IM$`Pr(>Chisq)`[2]
  resF$Chan[count]<-chan
  resF$Cond[count]<-3
  
  
}
#res$Cond<-as.factor(res$Cond)
#res$Chan<-as.factor(res$Chan)
write.table (resF,file="/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults_EEG_allCh3_results_F.txt",sep = ",",na = "NaN")


