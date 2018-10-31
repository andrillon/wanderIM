library('lme4')
#library('lmerTest')

t<-read.delim("/Users/tand0009/Data/WanderIM/behav/wanderIM_binary_table_ONOFF_MWMB.txt",sep = ",",na.strings = "NaN")
t$SubID<-as.factor(t$SubID)
t$State<-as.factor(t$State)
t$Task<-as.factor(t$Task)



mdl_0 <- glmer(ONOFF~1+(1|SubID), data=t, family=binomial)
mdl_1 <- glmer(ONOFF~nBlock+(1|SubID), data=t, family=binomial)
mdl_2 <- glmer(ONOFF~nBlock+nProbe+(1|SubID), data=t, family=binomial)
mdl_3 <- glmer(ONOFF~nBlock+nProbe+CorrGo+(1|SubID), data=t, family=binomial)
mdl_4 <- glmer(ONOFF~nBlock+nProbe+CorrGo+CorrNoGo+(1|SubID), data=t, family=binomial)
#mdl_5 <- glmer(ONOFF~nBlock+nProbe+CorrGo+CorrNoGo+dp+(1|SubID), data=t, family=binomial)
mdl_5 <- glmer(ONOFF~nBlock+nProbe+(CorrGo+CorrNoGo+RTGO)+(1|SubID), data=t, family=binomial)

mdl2_0 <- glmer(MWMB~nBlock+nProbe+(CorrGo+CorrNoGo+RTGO)+Task+(1|SubID), data=t, family=binomial)
mdl2_1 <- glmer(MWMB~Task*(CorrNoGo+RTGO)+Task+(1|SubID), data=t, family=binomial)
mdl2_2 <- glmer(MWMB~nBlock+nProbe+Task*CorrNoGo+RTGO+(1|SubID), data=t, family=binomial)
mdl2_3 <- glmer(MWMB~nBlock+nProbe+(RTGO)+(1|SubID), data=t, family=binomial)

mdl2_2 <- glmer(MWMB~nBlock+nProbe+(1|SubID), data=t, family=binomial)
mdl2_3 <- glmer(MWMB~nBlock+nProbe+CorrGo+(1|SubID), data=t, family=binomial)
mdl2_4 <- glmer(MWMB~nBlock+nProbe+CorrGo+CorrNoGo+(1|SubID), data=t, family=binomial)
#mdl_5 <- glmer(ONOFF~nBlock+nProbe+CorrGo+CorrNoGo+dp+(1|SubID), data=t, family=binomial)
mdl2_5 <- glmer(MWMB~nBlock+nProbe+(CorrGo+CorrNoGo+RTGO)+(1|SubID), data=t, family=binomial)


library('yarrr')
mycol.col <- c(
  rgb(252, 141, 98, maxColorValue = 255),    # ON
  rgb(102, 194, 165, maxColorValue = 255),    # MW
  rgb(141, 160, 203, maxColorValue = 255))   # MB
  
  
dat=t
dat$ONOFF<-as.factor(t$ONOFF)
dat$MWMB<-as.factor(t$MWMB)
levels(dat$ONOFF)<-c("OFF","ON")
levels(dat$MWMB)<-c("MB","MW")
levels(dat$State)<-c("ON","MW","MB")
levels(dat$Task)<-c("Face","Digit")

dat_F<-subset(dat,Task == "Face")
dat_D<-subset(dat,Task == "Digit")
par(mfrow=c(1,5))
pirateplot(formula = CorrGo ~ State,
           data = dat_F,
           theme = 1,
           ,
           main = "GO")
pirateplot(formula = CorrNoGo ~ State,
           data = dat_F,
           theme = 2,
           main = "NOGO")
pirateplot(formula = dp ~ State,
           data = dat_F,
           theme = 2,
           main = "d-prime")
pirateplot(formula = crit ~ State,
           data = dat_F,
           theme = 2,
           main = "crit")
pirateplot(formula = RTGO ~ State,
           data = dat_F,
           theme = 2,
           main = "RT(go)")

pirateplot(formula = RTGO ~ State* Task,
           data = dat,
           theme = 2,
           main = "RTGO")