library('lme4')

t<-read.delim("/Users/tand0009/Data/WanderIM/behav/WanderIM_TestResults.txt",sep = ",")
t$SubID<-as.factor(t$SubID)
t$TrCat<-as.factor(t$TrCat)
t$StimID<-as.factor(t$StimID)
t$Task<-as.factor(t$Task)


mdl_0 <- lmer(Corr ~ 1 + (1| SubID), data=t)
mdl_1 <- lmer(Corr ~ TrCat + (1| SubID), data=t)
mdl_2 <- lmer(Corr ~ TrCat+Task + (1| SubID), data=t)
mdl_3 <- lmer(Corr ~ TrCat*Task + (1| SubID), data=t)

#mdl_4 <- lmer(Corr ~ TrCat*Task+StimID + (1| SubID), data=t)

p<-read.delim("/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults.txt",sep = ",")
#SubID,nBlock,Task,nTrial,Look,State,Orig,Awa,Int,Eng,Perf,Vig,Corr,RT,TrCat,DistProbe
p$SubID<-as.factor(p$SubID)
p$State<-as.factor(p$State)
p$Orig<-as.factor(p$Orig)
p$Task<-as.factor(p$Task)
p$TrCat<-as.factor(p$TrCat)
p<-subset(p,DistProbe>-21)
p$State<-factor(p$State,c("2","1","3"))

mdl_0 <- lmer(Corr ~ 1 + (1| SubID), data=p)
mdl_1 <- lmer(Corr ~ TrCat + (1| SubID), data=p)
mdl_2 <- lmer(Corr ~ TrCat+Task + (1| SubID), data=p)
mdl_3 <- lmer(Corr ~ TrCat*Task + (1| SubID), data=p)
mdl_4 <- lmer(Corr ~ TrCat*Task+DistProbe + (1| SubID), data=p)

mdl_5 <- lmer(Corr ~ TrCat*Task+DistProbe+State + (1| SubID), data=p)
mdl_6 <- lmer(Corr ~ TrCat*Task*State+DistProbe + (1| SubID), data=p)
