library('lme4')
#library('lmerTest')

p<-read.delim("/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults_EEG_allCh.txt",sep = ",")
#SubID,nBlock,Task,nTrial,Look,State,Orig,Awa,Int,Eng,Perf,Vig,Corr,RT,TrCat,DistProbe
p$SubID<-as.factor(p$SubID)
p$State<-as.factor(p$State)
p$Orig<-as.factor(p$Orig)
p$Task<-as.factor(p$Task)
p$L_freq<-as.factor(p$L_freq)
p$R_freq<-as.factor(p$R_freq)
p$Chan<-as.factor(p$Chan)
p$Freq<-as.factor(p$Freq)

p<-subset(p,State != "4")
p$State<-droplevels(p$State)
p$State<-factor(p$State,c("2","1","3"))


mdl1_1 <- lmer(SNR ~ 1 + (1| SubID), data=p)
mdl1_2 <- lmer(SNR ~ Freq + (1| SubID), data=p)
mdl1_3 <- lmer(SNR ~ Freq + Chan + (1| SubID), data=p)
mdl1_4 <- lmer(SNR ~ Freq + Chan + State + (1| SubID), data=p)
mdl1_5 <- lmer(SNR ~ Freq*State + Chan*State + (1| SubID), data=p)
mdl1_6 <- lmer(SNR ~ Freq*State*Chan + (1| SubID), data=p)

pF<-subset(p,Freq == 1 | Freq == 2)
p2F<-subset(p,Freq == 3 | Freq == 4)
pIM<-subset(p,Freq == 5)

mdl2_1 <- lmer(SNR ~ 1 + (1| SubID), data=pF)
mdl2_2 <- lmer(SNR ~ Chan + (1| SubID), data=pF)
mdl2_3 <- lmer(SNR ~ Freq + State + (1| SubID), data=pF)
mdl2_4 <- lmer(SNR ~ Freq + Chan + State + (1| SubID), data=pF)
mdl2_5 <- lmer(SNR ~ Freq*State + Chan*State + (1| SubID), data=pF)
mdl2_6 <- lmer(SNR ~ Freq*State*Chan + (1| SubID), data=pF)

mdl3_1 <- lmer(SNR ~ 1 + (1| SubID), data=pIM)
mdl3_2 <- lmer(SNR ~ Chan + (1| SubID), data=pIM)
mdl3_3 <- lmer(SNR ~ Chan + State + (1| SubID), data=pIM)
mdl3_4 <- lmer(SNR ~ Chan * State + (1| SubID), data=pIM)

