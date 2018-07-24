#### Import Data ####
WanderIM_TestResults <- read.csv("~/Google Drive/GitHub/data/behav/WanderIM_ProbeResults.txt")
View(WanderIM_TestResults)

WanderIM_ProbeResults <- read.csv("~/Google Drive/GitHub/data/behav/WanderIM_ProbeResults.txt")
View(WanderIM_ProbeResults)

NAB.df <- read.delim("~/Google Drive/GitHub/data/NAB_workplace/Shortened Participant Output.txt")
View(NAB.df)

Survey.df <- read.csv("~/Google Drive/GitHub/wanderIM/survey/data/output/MWI_Computed_Scales.csv")
View(Survey.df)

#Recode performance variable from nab.df
library(stringr)
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 

NAB.df$EOY2017Rating_num <- numextract(NAB.df$EOY2017Rating)

NAB.df$EOY2017Rating_int[NAB.df$EOY2017Rating_num=="1"] <- 1
NAB.df$EOY2017Rating_int[NAB.df$EOY2017Rating_num=="2"] <- 2
NAB.df$EOY2017Rating_int[NAB.df$EOY2017Rating_num=="3"] <- 3
NAB.df$EOY2017Rating_int[NAB.df$EOY2017Rating_num=="4"] <- 4
NAB.df$EOY2017Rating_int[NAB.df$EOY2017Rating_num=="5"] <- 5

as.numeric(NAB.df$EOY2017Rating_int)
is.numeric(NAB.df$EOY2017Rating_int)

library("psych")


# EXPLORE BEHAVIOURAL DATA

#### OVERALL CORRECTNESS ACROSS TRIALS by NBLOCK ####

#Reshape data to view correctness by individual participant across trials
tbl_perf_by_block_ind <- table(WanderIM_ProbeResults$nBlock, WanderIM_ProbeResults$SubID, WanderIM_ProbeResults$Corr)
wrongs <- tbl_perf_by_block_ind[1:6, 1:20, 1]
rights <- tbl_perf_by_block_ind[1:6, 1:20, 2]
tbl_perf_by_block_ind <- rights / (rights+wrongs)
tbl_perf_by_block_ind <- as.data.frame.matrix(tbl_perf_by_block_ind)
tbl_perf_by_block_ind$block <- c(1, 2, 3, 4, 5, 6)
library(reshape2)
tbl_perf_by_block_ind <- melt(tbl_perf_by_block_ind, id.vars = c("block"))
names(tbl_perf_by_block_ind)[names(tbl_perf_by_block_ind)=="variable"] <- "participant"
names(tbl_perf_by_block_ind)[names(tbl_perf_by_block_ind)=="value"] <- "perc_correct"

#plot correctness by participant across trials

#Violin plot
library(ggplot2)
tbl_perf_by_block_ind$block <- as.factor(tbl_perf_by_block_ind$block)
perf_block_ind_violin <- ggplot(tbl_perf_by_block_ind, aes(x=block, y=perc_correct)) + 
  geom_violin()
perf_block_ind_violin + stat_summary(fun.y=mean, geom="point", shape=23, size=2)

#boxplot
boxplot(perc_correct~block,data=tbl_perf_by_block_ind, main="Percent Correct by Block Number", 
        xlab="Block Number", ylab="Percent Correct")

#Calculate mean scores across trials
aggregate(perc_correct~block, data=tbl_perf_by_block_ind, FUN="mean")

# Summary statistics
describeBy(tbl_perf_by_block_ind, group = tbl_perf_by_block_ind$block, mat=FALSE)

#Statistical analysis - effect of block
aov_effect_block <- aov(perc_correct~block, data=tbl_perf_by_block_ind)
summary(aov_effect_block)
#No effect of block on correctness



#### CORRECTNESS ACROSS NOGO TRIALS by NBLOCK ####

#Reshape data to view correctness by individual participant across trials
tbl_perf_by_block_ind_nogo <- subset(WanderIM_ProbeResults, TrCat == "1")

tbl_perf_by_block_ind_nogo <- table(tbl_perf_by_block_ind_nogo$nBlock, tbl_perf_by_block_ind_nogo$SubID, tbl_perf_by_block_ind_nogo$Corr)
wrongs <- tbl_perf_by_block_ind_nogo[1:6, 1:20, 1]
rights <- tbl_perf_by_block_ind_nogo[1:6, 1:20, 2]
tbl_perf_by_block_ind_nogo <- rights / (rights+wrongs)
tbl_perf_by_block_ind_nogo <- as.data.frame.matrix(tbl_perf_by_block_ind_nogo)
tbl_perf_by_block_ind_nogo$block <- c(1, 2, 3, 4, 5, 6)
library(reshape2)
tbl_perf_by_block_ind_nogo <- melt(tbl_perf_by_block_ind_nogo, id.vars = c("block"))
names(tbl_perf_by_block_ind_nogo)[names(tbl_perf_by_block_ind_nogo)=="variable"] <- "participant"
names(tbl_perf_by_block_ind_nogo)[names(tbl_perf_by_block_ind_nogo)=="value"] <- "perc_correct"

#plot correctness by participant across trials

#Violin plot
library(ggplot2)
tbl_perf_by_block_ind_nogo$block <- as.factor(tbl_perf_by_block_ind_nogo$block)
perf_block_ind_violin <- ggplot(tbl_perf_by_block_ind_nogo, aes(x=block, y=perc_correct)) + 
  geom_violin()
perf_block_ind_violin + stat_summary(fun.y=mean, geom="point", shape=23, size=2)

#boxplot
boxplot(perc_correct~block,data=tbl_perf_by_block_ind_nogo, main="Percent Correct by Block Number", 
        xlab="Block Number", ylab="Percent Correct")

#Calculate mean scores across trials
aggregate(perc_correct~block, data=tbl_perf_by_block_ind_nogo, FUN="mean")


# Summary statistics
describeBy(tbl_perf_by_block_ind_nogo, group = tbl_perf_by_block_ind_nogo$block, mat=FALSE)

#Statistical analysis - effect of block
aov_effect_block <- aov(perc_correct~block, data=tbl_perf_by_block_ind_nogo)
summary(aov_effect_block)
#No effect of block on correctness


#### CORRECTNESS ACROSS GO TRIALS by NBLOCK ####

#Reshape data to view correctness by individual participant across trials
tbl_perf_by_block_ind_go <- subset(WanderIM_ProbeResults, TrCat == "0")

tbl_perf_by_block_ind_go <- table(tbl_perf_by_block_ind_go$nBlock, tbl_perf_by_block_ind_go$SubID, tbl_perf_by_block_ind_go$Corr)
wrongs <- tbl_perf_by_block_ind_go[1:6, 1:20, 1]
rights <- tbl_perf_by_block_ind_go[1:6, 1:20, 2]
tbl_perf_by_block_ind_go <- rights / (rights+wrongs)
tbl_perf_by_block_ind_go <- as.data.frame.matrix(tbl_perf_by_block_ind_go)
tbl_perf_by_block_ind_go$block <- c(1, 2, 3, 4, 5, 6)
library(reshape2)
tbl_perf_by_block_ind_go <- melt(tbl_perf_by_block_ind_go, id.vars = c("block"))
names(tbl_perf_by_block_ind_go)[names(tbl_perf_by_block_ind_go)=="variable"] <- "participant"
names(tbl_perf_by_block_ind_go)[names(tbl_perf_by_block_ind_go)=="value"] <- "perc_correct"

#plot correctness by participant across trials

#Violin plot
library(ggplot2)
tbl_perf_by_block_ind_go$block <- as.factor(tbl_perf_by_block_ind_go$block)
perf_block_ind_violin <- ggplot(tbl_perf_by_block_ind_go, aes(x=block, y=perc_correct)) + 
  geom_violin()
perf_block_ind_violin + stat_summary(fun.y=mean, geom="point", shape=23, size=2)

#boxplot
boxplot(perc_correct~block,data=tbl_perf_by_block_ind_go, main="Percent Correct by Block Number", 
        xlab="Block Number", ylab="Percent Correct")

#Calculate mean scores across trials
aggregate(perc_correct~block, data=tbl_perf_by_block_ind_go, FUN="mean")

# Summary statistics
describeBy(tbl_perf_by_block_ind_go, group = tbl_perf_by_block_ind_go$block, mat=FALSE)

#Statistical analysis - effect of block
aov_effect_block <- aov(perc_correct~block, data=tbl_perf_by_block_ind_go)
summary(aov_effect_block)
#No effect of block on correctness





#### REACTION TIME ACROSS GO TRIALS by NBLOCK ####

#Reshape data to view response time by individual participant across trials
tbl_RT_by_block_ind_go <- subset(WanderIM_ProbeResults, TrCat == "0")

tbl_RT_by_block_ind_go <- aggregate(tbl_RT_by_block_ind_go$RT~tbl_RT_by_block_ind_go$nBlock+tbl_RT_by_block_ind_go$SubID, FUN="mean")

#plot RT by participant across trials

#Violin plot
library(ggplot2)
tbl_RT_by_block_ind_go$`tbl_RT_by_block_ind_go$nBlock` <- as.factor(tbl_RT_by_block_ind_go$`tbl_RT_by_block_ind_go$nBlock`)

perf_block_ind_violin <- ggplot(tbl_RT_by_block_ind_go, aes(x=`tbl_RT_by_block_ind_go$nBlock`, y=`tbl_RT_by_block_ind_go$RT`)) + 
  geom_violin()
perf_block_ind_violin + stat_summary(fun.y=mean, geom="point", shape=23, size=2)

#boxplot
boxplot(`tbl_RT_by_block_ind_go$RT`~`tbl_RT_by_block_ind_go$nBlock`,data=tbl_RT_by_block_ind_go, main="Reaction time by Block Number", 
        xlab="Block Number", ylab="Reaction Time")

#Calculate mean scores across trials
aggregate(`tbl_RT_by_block_ind_go$RT`~`tbl_RT_by_block_ind_go$nBlock`, data=tbl_RT_by_block_ind_go, FUN="mean")

# Summary statistics
describeBy(tbl_RT_by_block_ind_go, group = tbl_RT_by_block_ind_go$block, mat=FALSE)

#Statistical analysis - effect of block
aov_effect_block <- aov(`tbl_RT_by_block_ind_go$RT`~`tbl_RT_by_block_ind_go$nBlock`, data=tbl_RT_by_block_ind_go)
summary(aov_effect_block)
#No effect of block on reaction time







#### CORRECTNESS BY TASK ####

#Reshape data to view correctness by individual participant across trials
#Task coding: 1=faces, 2=squares - taken from Thomas' comments in matlab script

tbl_perf_task_ind <- table(WanderIM_ProbeResults$Task, WanderIM_ProbeResults$SubID, WanderIM_ProbeResults$Corr)
wrongs <- tbl_perf_task_ind[1:2, 1:20, 1]
rights <- tbl_perf_task_ind[1:2, 1:20, 2]
tbl_perf_task_ind <- rights / (rights+wrongs)
tbl_perf_task_ind <- as.data.frame.matrix(tbl_perf_task_ind)
tbl_perf_task_ind$Task <- c(1, 2)
library(reshape2)
tbl_perf_task_ind <- melt(tbl_perf_task_ind, id.vars = c("Task"))
names(tbl_perf_task_ind)[names(tbl_perf_task_ind)=="variable"] <- "participant"
names(tbl_perf_task_ind)[names(tbl_perf_task_ind)=="value"] <- "perc_correct"

#plot correctness by participant across tasks

#Violin plot
library(ggplot2)
tbl_perf_task_ind$Task <- as.factor(tbl_perf_task_ind$Task)
perf_block_ind_violin <- ggplot(tbl_perf_task_ind, aes(x=Task, y=perc_correct)) + 
  geom_violin()
perf_block_ind_violin + stat_summary(fun.y=mean, geom="point", shape=23, size=2)

#boxplot
boxplot(perc_correct~Task,data=tbl_perf_task_ind, main="Percent Correct by Block Number", 
        xlab="Task (1=faces, 2=Squares)", ylab="Percent Correct")

#Calculate mean scores across trials
aggregate(perc_correct~Task, data=tbl_perf_task_ind, FUN="mean")

# Summary statistics
describeBy(tbl_perf_task_ind, group = tbl_perf_task_ind$Task, mat=FALSE)

#Statistical analysis - effect of block
aov_effect_Task <- aov(perc_correct~Task, data=tbl_perf_task_ind)
summary(aov_effect_Task)
TukeyHSD(aov_effect_Task)
#No effect of block on correctness





#### REACTION TIME BY TASK ####
#Task coding: 1=faces, 2=squares - taken from Thomas' comments in matlab script


#Reshape data to view response time by individual participant across trials
tbl_RT_by_task_ind_go <- subset(WanderIM_ProbeResults, TrCat == "0")

tbl_RT_by_task_ind_go <- aggregate(tbl_RT_by_task_ind_go$RT~tbl_RT_by_task_ind_go$Task+tbl_RT_by_task_ind_go$SubID, FUN="mean")


#plot RT by participant across trials

#Violin plot
library(ggplot2)
tbl_RT_by_task_ind_go$`tbl_RT_by_task_ind_go$Task` <- as.factor(tbl_RT_by_task_ind_go$`tbl_RT_by_task_ind_go$Task`)

tbl_RT_by_task_ind_go_violin <- ggplot(tbl_RT_by_task_ind_go, aes(x=`tbl_RT_by_task_ind_go$Task`, y=`tbl_RT_by_task_ind_go$RT`)) + 
  geom_violin()
tbl_RT_by_task_ind_go_violin + stat_summary(fun.y=mean, geom="point", shape=23, size=2)

#boxplot
boxplot(`tbl_RT_by_task_ind_go$RT`~`tbl_RT_by_task_ind_go$Task`,data=tbl_RT_by_task_ind_go, main="Reaction time by Task", 
        xlab="Task (1=faces, 2=Squares)", ylab="Reaction Time")

#Calculate mean scores across trials
aggregate(`tbl_RT_by_task_ind_go$RT`~`tbl_RT_by_task_ind_go$Task`, data=tbl_RT_by_task_ind_go, FUN="mean")

# Summary statistics
describeBy(tbl_RT_by_task_ind_go, group = tbl_RT_by_task_ind_go$`tbl_RT_by_task_ind_go$Task`, mat=FALSE)

#Statistical analysis - effect of block
aov_RTeffect_Task <- aov(`tbl_RT_by_task_ind_go$RT`~`tbl_RT_by_task_ind_go$Task`, data=tbl_RT_by_task_ind_go)
summary(aov_RTeffect_Task)
#No effect of block on reaction time




#### OVERALL CORRECTNESS BY TRIAL N WITHIN BLOCK - issue viewing this - too few instances to make sense ####

tbl_perf_by_trial_ind <- table(WanderIM_TestResults$nTrial, WanderIM_TestResults$Corr)
wrongs <- tbl_perf_by_trial_ind[1:474, 1]
rights <- tbl_perf_by_trial_ind[1:474, 2]
tbl_perf_by_trial_ind <- rights / (rights+wrongs)

remove(tbl_perf_by_trial_ind)

tbl_perf_by_trial_ind <- melt(tbl_perf_by_trial_ind, id.vars = c("nTrial"))
names(tbl_perf_by_trial_ind)[names(tbl_perf_by_trial_ind)=="value"] <- "perc_correct"

plot(tbl_perf_by_trial_ind)


#### REACTION TIME BY TRIAL N WITHIN BLOCK ####
#Task coding: 1=faces, 2=squares - taken from Thomas' comments in matlab script

#Reshape data to view response time by individual participant across trials
tbl_RT_by_trial_ind_go <- subset(WanderIM_ProbeResults, TrCat == "0")

tbl_RT_by_trial_ind_go <- aggregate(tbl_RT_by_trial_ind_go$RT~tbl_RT_by_trial_ind_go$nTrial+tbl_RT_by_trial_ind_go$SubID, FUN="mean")

plot(tbl_RT_by_trial_ind_go$`tbl_RT_by_trial_ind_go$nTrial`, tbl_RT_by_trial_ind_go$`tbl_RT_by_trial_ind_go$RT`)
abline(lm(tbl_RT_by_trial_ind_go$`tbl_RT_by_trial_ind_go$RT`~tbl_RT_by_trial_ind_go$`tbl_RT_by_trial_ind_go$nTrial`), col="red") # regression line (y~x) 
lines(lowess(tbl_RT_by_trial_ind_go$`tbl_RT_by_trial_ind_go$nTrial`, tbl_RT_by_trial_ind_go$`tbl_RT_by_trial_ind_go$RT`), col="blue") # lowess line (x,y)



#### OVERALL CORRECTNESS BY PROBE (10 TRIALS BEFORE PROBE) - finish lme analysis ####

#Remove all trials >-10 away from probe
tbl_perf_by_probe <- subset(WanderIM_ProbeResults, DistProbe >-11)

# Aggregate to correctness, including relevant variables
attach(tbl_perf_by_probe)
tbl_perf_by_probe <- aggregate(Corr~SubID+nBlock+Task+nTrial+Look+State+Orig+Awa+Int+Eng+Perf+Vig, FUN="mean")
detach(tbl_perf_by_probe)

# Analysis Across Probe Type
agg1_perf_probe <- aggregate(Corr~State, FUN="mean")
agg1_perf_probe
plot(agg1_perf_probe)

agg2_perf_probe <- aggregate(Corr~State+SubID, FUN="mean")
agg2_perf_probe
plot(agg2_perf_probe)

#Violin plot
library(ggplot2)
agg2_perf_probe$State <- as.factor(agg2_perf_probe$State)

agg2_perf_probe_violin <- ggplot(agg2_perf_probe, aes(x=agg2_perf_probe$State, y=agg2_perf_probe$Corr)) + 
  geom_violin()
agg2_perf_probe_violin + stat_summary(fun.y=mean, geom="point", shape=23, size=2)

#boxplot
boxplot(agg2_perf_probe$Corr~agg2_perf_probe$State,data=agg2_perf_probe, main="Correctness by probe", 
        xlab="Probe Type (1=ON, 2=MW, 3=MB, 4=UN)", ylab="Correctness %")

#Statistical analysis - effect of STATE (mixed models analysis)



#### REACTION TIME BY PROBE (10 TRIALS BEFORE PROBE) - finish lme analysis ####

#Remove all trials >-10 away from probe
tbl_RT_by_probe <- subset(WanderIM_ProbeResults, DistProbe >-11)

# Aggregate to correctness, including relevant variables
attach(tbl_RT_by_probe)
tbl_RT_by_probe <- aggregate(RT~SubID+nBlock+Task+nTrial+Look+State+Orig+Awa+Int+Eng+Perf+Vig, FUN="mean")
detach(tbl_RT_by_probe)

# Analysis Across Probe Type
agg1_RT_probe <- aggregate(RT~State, FUN="mean")
agg1_RT_probe
plot(agg1_RT_probe)

agg2_RT_probe <- aggregate(RT~State+SubID, FUN="mean")
agg2_RT_probe
plot(agg2_RT_probe)

#Violin plot
library(ggplot2)
agg2_RT_probe$State <- as.factor(agg2_RT_probe$State)

agg2_RT_probe_violin <- ggplot(agg2_RT_probe, aes(x=agg2_RT_probe$State, y=agg2_RT_probe$RT)) + 
  geom_violin()
agg2_RT_probe_violin + stat_summary(fun.y=mean, geom="point", shape=23, size=2)

#boxplot
boxplot(RT~State,data=agg2_RT_probe, main="Reaction Time by probe", 
        xlab="Probe Type (1=ON, 2=MW, 3=MB, 4=UN)", ylab="Reaction Time (ms)")

#Statistical analysis - effect of STATE (mixed models analysis)




#### MODELLING IMPACT OF TASK PARAMETERS - finish lme analysis ####

install.packages("nlme")
library("nlme")
lme0 <- lme(WanderIM_ProbeResults$Corr~1, (1|WanderIM_ProbeResults$SubID))


## From Thomas' matlab script:
#%% Examples of models (to do in R)
#lme_0= fitlme(tbl_probe,'Corr~1+(1|SubID)');
#lme_1= fitlme(tbl_probe,'Corr~1+(1|SubID)');

#lme_full= fitlme(tbl_probe,'Corr~Task+nBlock+nTrial+TrCat+DistProbe+State+Task*State+(1|SubID)');
#lme_full2= fitlme(tbl_probe,'Corr~Task+nBlock+nTrial+TrCat+DistProbe+State+DistProbe*State+Task*State+(1|SubID)');

#[v] = compare(lme_full,lme_full2);
#if v.pValue(end)<0.05
#fprintf('More complex model (%s) wins (pV=%1.6f)\n',lme_full2.Formula,v.pValue(end))
#else
#  fprintf('Less complex model (%s) wins (pV=%1.6f)\n',lme_full.Formula,v.pValue(end))
#end






#### SART PERFORMANCE AND RELATIONSHIP TO QUESTIONNAIRES ####

# Is overall performance on the SART impacted by the the variables we've measured?

#Reshape data to view correctness by individual participant

# Aggregate to correctness by participant ID
attach(WanderIM_ProbeResults)
tbl_perf_ind <- aggregate(Corr~SubID, FUN="mean")
tbl_RT_ind <- aggregate(RT~SubID, FUN="mean")
detach(WanderIM_ProbeResults)

#Join survey and outcome variables
tbl_perf_ind_merged <- merge(tbl_perf_ind,Survey.df, by.x="SubID", by.y="Participant_No")
tbl_perf_ind_merged <- merge(tbl_perf_ind_merged,NAB.df, by.x="SubID", by.y="Participant_No")
tbl_perf_ind_merged <- merge(tbl_perf_ind_merged,tbl_RT_ind, by.x="SubID", by.y="SubID")

#myvars = c(tbl_perf_ind_merged$Corr, tbl_perf_ind_merged$D2...Age, tbl_perf_ind_merged$DASS_Depression, tbl_perf_ind_merged$DASS_Anxiety, tbl_perf_ind_merged$DASS_Stress,
#                                          tbl_perf_ind_merged$ESS, tbl_perf_ind_merged$MWQ_mean, tbl_perf_ind_merged$MiniIPIP_Agreeableness, tbl_perf_ind_merged$MiniIPIP_Openness, tbl_perf_ind_merged$MiniIPIP_Extraversion, tbl_perf_ind_merged$MiniIPIP_Conscientiousness, tbl_perf_ind_merged$MiniIPIP_Neuroticism,
#                                          tbl_perf_ind_merged$PANAS_Positive_Affect, tbl_perf_ind_merged$PANAS_Negative_Affect, tbl_perf_ind_merged$ASRS_ptA_num, tbl_perf_ind_merged$EOY2017Rating, tbl_perf_ind_merged$Employee.Group, tbl_perf_ind_merged$Org.Tenure, tbl_perf_ind_merged$FTE, tbl_perf_ind_merged$Promotions_per_year)

library("dplyr")
tbl_numeric_perf_ind_merged <- select_if(tbl_perf_ind_merged, is.numeric)

#install.packages("PerformanceAnalytics")
#library("PerformanceAnalytics")
#chart.Correlation(tbl_numeric_perf_ind_merged, histogram=TRUE, pch=19)

## Create correlation matrix

cormatrix_numeric_perf_ind_merged <- round(cor(tbl_numeric_perf_ind_merged, method="spearman"),2)

library("Hmisc")
cormatrix_numeric_perf_ind_merged <- rcorr(as.matrix(cormatrix_numeric_perf_ind_merged, type ="spearman"))

# Extract the correlation coefficients
cormatrix_numeric_perf_ind_merged$r
# Extract p-values
cormatrix_numeric_perf_ind_merged$P

write.table(cormatrix_numeric_perf_ind_merged$r, file = "./data/output/cormatrix_numeric_perf_ind_merged_r.csv",row.names=FALSE, na="",col.names=TRUE, sep=",")
write.table(cormatrix_numeric_perf_ind_merged$P, file = "./data/output/cormatrix_numeric_perf_ind_merged_p.csv",row.names=FALSE, na="",col.names=TRUE, sep=",")

#Population too small to interpret promotions per year. It is not correlated with performance and therefore not a good measure of performance outcomes.


#### PROBE RESULTS AND RELATIONSHIP TO QUESTIONNAIRES ####

#### end ####





Hmisc::summarize()
aggregate(WanderIM_ProbeResults$Corr, by=WanderIM_ProbeResults$nBlock, FUN="mean") 

plot(WanderIM_ProbeResults$Corr, WanderIM_ProbeResults$nBlock, main="Correctness by block", 
     xlab="Correct", ylab="Block number", pch=19)

