library('lme4')
#library('lmerTest')

t<-read.delim("/Users/tand0009/Data/WanderIM/behav/WanderIM_TestResults.txt",sep = ",")
t$SubID<-as.factor(t$SubID)
t$TrCat<-as.factor(t$TrCat)
t$StimID<-as.factor(t$StimID)
t$Task<-as.factor(t$Task)

GO<-subset(t,TrCat=="0")
NOGO<-subset(t,TrCat=="1")


mdl_0 <- lmer(Corr ~ 1 + (1| SubID), data=t)
mdl_1 <- lmer(Corr ~ TrCat + (1| SubID), data=t)
mdl_2 <- lmer(Corr ~ TrCat+Task + (1| SubID), data=t)
mdl_3 <- lmer(Corr ~ TrCat+Task + (1 + TrCat + Task| SubID), data=t)
mdl_4 <- lmer(Corr ~ TrCat*Task + (1 + TrCat + Task| SubID), data=t)


mdl2_0 <- lmer(Corr ~ 1 + (1| SubID), data=GO)
mdl2_1 <- lmer(Corr ~ Task + (1| SubID), data=GO)
mdl2_2 <- lmer(Corr ~ Task + (1+Task| SubID), data=GO)

mdl3_0 <- lmer(Corr ~ 1 + (1| SubID), data=NOGO)
mdl3_1 <- lmer(Corr ~ Task + (1| SubID), data=NOGO)
mdl3_2 <- lmer(Corr ~ Task + (1+Task| SubID), data=NOGO)


#mdl_4 <- lmer(Corr ~ TrCat*Task+StimID + (1| SubID), data=t)

p<-read.delim("/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults.txt",sep = ",")
#SubID,nBlock,Task,nTrial,Look,State,Orig,Awa,Int,Eng,Perf,Vig,Corr,RT,TrCat,DistProbe
p$SubID<-as.factor(p$SubID)
p$State<-as.factor(p$State)
p$Orig<-as.factor(p$Orig)
p$Task<-as.factor(p$Task)
p$TrCat<-as.factor(p$TrCat)
p<-subset(p,(DistProbe>-19 & TrCat=="0") | (DistProbe>-3 & TrCat=="1"))
p<-subset(p,State != "4")
p$State<-droplevels(p$State)
p$State<-factor(p$State,c("1","2","3"))

pGO<-subset(p,TrCat=="0")
pNOGO<-subset(p,TrCat=="1")

mdl4_0 <- lmer(Corr ~ 1 + (1| SubID), data=pGO)
mdl4_1 <- lmer(Corr ~ Task + (1| SubID), data=pGO)
mdl4_2 <- lmer(Corr ~ Task + (1+Task| SubID), data=pGO)
mdl4_3 <- lmer(Corr ~ Task + State + (1+Task| SubID), data=pGO)
mdl4_4 <- lmer(Corr ~ Task + State + (1+Task+State| SubID), data=pGO)
mdl4_5 <- lmer(Corr ~ Task*State + (1+Task+State| SubID), data=pGO)
mdl4_6 <- lmer(Corr ~ Task*State*DistProbe + (1+Task+State| SubID), data=pGO)

mdl5_0 <- lmer(Corr ~ 1 + (1| SubID), data=pNOGO)
mdl5_1 <- lmer(Corr ~ 1 + (1+Task| SubID), data=pNOGO)
mdl5_2 <- lmer(Corr ~ State + (1+Task| SubID), data=pNOGO)
mdl5_3 <- lmer(Corr ~ Task*State + (1+Task| SubID), data=pNOGO)
#mdl5_4 <- lmer(Corr ~ Task*State + DistProbe + (1+Task| SubID), data=pNOGO)
#mdl5_5 <- lmer(Corr ~ Task*State*DistProbe + (1+Task| SubID), data=pNOGO)
#mdl5_5 <- lmer(Corr ~ Task*State + DistProbe + (1+Task+DistProbe| SubID), data=pNOGO)

p2<-p
p2$State<-droplevels(p2$State)
p2$State<-factor(p2$State,c("2","1","3"))
p2GO<-subset(p2,TrCat=="0")
p2NOGO<-subset(p2,TrCat=="1")

mdl6_0 <- lmer(Corr ~ 1 + (1| SubID), data=p2GO)
mdl6_1 <- lmer(Corr ~ Task + (1| SubID), data=p2GO)
mdl6_2 <- lmer(Corr ~ Task + (1+Task| SubID), data=p2GO)
mdl6_3 <- lmer(Corr ~ Task + State + (1+Task| SubID), data=p2GO)
mdl6_4 <- lmer(Corr ~ Task + State + (1+Task+State| SubID), data=p2GO)
mdl6_5 <- lmer(Corr ~ Task*State + (1+Task+State| SubID), data=p2GO)
mdl6_6 <- lmer(Corr ~ Task*State*DistProbe + (1+Task+State| SubID), data=p2GO)

mdl7_0 <- lmer(Corr ~ 1 + (1| SubID), data=p2NOGO)
mdl7_1 <- lmer(Corr ~ 1 + (1+Task| SubID), data=p2NOGO)
mdl7_2 <- lmer(Corr ~ State + (1+Task| SubID), data=p2NOGO)
mdl7_3 <- lmer(Corr ~ Task*State + (1+Task| SubID), data=p2NOGO)
