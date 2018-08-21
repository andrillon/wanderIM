library('lme4')

t<-read.table("/Users/tand0009/Data/WanderIM/behav/WanderIM_ProbeResults_HB.txt", header = TRUE,
              sep = ",", quote = "\"'",
              dec = ".",
              na.strings = "<undefined>",
              colClasses=c("factor","numeric","factor","numeric","factor","factor","factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

#SubID,nBlock,Task,nTrial,Look,State,Orig,Awa,Int,Eng,Perf,Vig,Go,NoGo,RT,HB,HB2
mdl_0 <- lmer(HB ~ 1 + (1| SubID), data=t)
mdl_1 <- lmer(Go ~ Orig + (1| SubID), data=t)


