library('lme4')
#library('lmerTest')

p<-read.delim("/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults_EEG.txt",sep = ",")
#SubID,nBlock,Task,nTrial,Look,State,Orig,Awa,Int,Eng,Perf,Vig,Corr,RT,TrCat,DistProbe
p$SubID<-as.factor(p$SubID)
p$State<-as.factor(p$State)
p$Orig<-as.factor(p$Orig)
p$Task<-as.factor(p$Task)
p$L_freq<-as.factor(p$L_freq)
p$R_freq<-as.factor(p$R_freq)
p<-subset(p,State != "4")
p$State<-droplevels(p$State)
p$State<-factor(p$State,c("1","2","3"))

mdl4_0 <- lmer(pOZ_f1 ~ 1 + (1| SubID), data=p)
mdl4_1 <- lmer(pOZ_f1 ~ Task + (1+Task| SubID), data=p)
mdl4_1 <- lmer(pOZ_f1 ~ Task + (1+Task| SubID), data=p)
