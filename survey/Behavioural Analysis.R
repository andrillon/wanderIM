#### Import Data ####
WanderIM_TestResults <- read.csv("~/Google Drive/GitHub/data/behav/WanderIM_ProbeResults.txt")
View(WanderIM_TestResults)
WanderIM_TestResults$SubID <- as.factor(WanderIM_TestResults$SubID)
WanderIM_TestResults$Task <- as.factor(WanderIM_TestResults$Task)
WanderIM_TestResults$Look <- as.factor(WanderIM_TestResults$Look)
WanderIM_TestResults$State <- as.factor(WanderIM_TestResults$State)
WanderIM_TestResults$Awa <- as.ordered(WanderIM_TestResults$Awa)
WanderIM_TestResults$Eng <- as.ordered(WanderIM_TestResults$Eng)
WanderIM_TestResults$Perf <- as.ordered(WanderIM_TestResults$Perf)
WanderIM_TestResults$Vig <- as.ordered(WanderIM_TestResults$Vig)
WanderIM_TestResults$TrCat <- as.factor(WanderIM_TestResults$TrCat)

WanderIM_ProbeResults <- read.csv("~/Google Drive/GitHub/data/behav/WanderIM_ProbeResults.txt")
View(WanderIM_ProbeResults)
WanderIM_ProbeResults$SubID <- as.factor(WanderIM_ProbeResults$SubID)
WanderIM_ProbeResults$Task <- as.factor(WanderIM_ProbeResults$Task)
WanderIM_ProbeResults$Look <- as.factor(WanderIM_ProbeResults$Look)
WanderIM_ProbeResults$State <- as.factor(WanderIM_ProbeResults$State)
WanderIM_ProbeResults$Awa <- as.numeric(WanderIM_ProbeResults$Awa)
WanderIM_ProbeResults$Int <- as.numeric(WanderIM_ProbeResults$Int)
WanderIM_ProbeResults$Eng <- as.numeric(WanderIM_ProbeResults$Eng)
WanderIM_ProbeResults$Perf <- as.numeric(WanderIM_ProbeResults$Perf)
WanderIM_ProbeResults$Vig <- as.numeric(WanderIM_ProbeResults$Vig)
WanderIM_ProbeResults$TrCat <- as.factor(WanderIM_ProbeResults$TrCat)


NAB.df <- read.delim("~/Google Drive/GitHub/data/NAB_workplace/Shortened Participant Output.txt")
View(NAB.df)
as.factor(NAB.df$Participant_No)


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

# Import Manually scored Guilford's AUT data
library(readr)
Guilford_AUT <- read_csv("data/input/AUT Input R.csv")

# Import OSPAN
OSPAN.dat <- read_csv("data/input/automatedospan_summary_18_07_09.csv")

# Import behavioural response metrics exported from MATLAB
behav_resp_DV <- read_csv("data/output/all_behav_mat2.csv")

# EXPLORE BEHAVIOURAL DATA

#### Plots more logical to do in R ####
# Analysis of spider plot

#Remove all trials >-1 away from probe
temp_probe_results <- subset(WanderIM_ProbeResults, DistProbe >-2)

## Frequency of attention states
# Aggregate to Frequency of MB/MW/ON

attach(temp_probe_results)
temp_probe_results1 <- aggregate(nTrial~SubID+State, FUN = "length")
temp_probe_results2 <- aggregate(Awa~SubID+State, FUN = "mean")
temp_probe_results3 <- aggregate(Int~SubID+State, FUN = "mean")
temp_probe_results4 <- aggregate(Eng~SubID+State, FUN = "mean")
temp_probe_results5 <- aggregate(Perf~SubID+State, FUN = "mean")
temp_probe_results6 <- aggregate(Vig~SubID+State, FUN = "mean")

detach(temp_probe_results)
#Pivot so attention each attention state is own variable
library(reshape2)
#temp_probe_results1 <- dcast(temp_probe_results1, SubID~State)
#temp_probe_results2 <- dcast(temp_probe_results2, SubID~State)
#temp_probe_results3 <- dcast(temp_probe_results3, SubID~State)
#temp_probe_results4 <- dcast(temp_probe_results4, SubID~State)
#temp_probe_results5 <- dcast(temp_probe_results5, SubID~State)
#temp_probe_results6 <- dcast(temp_probe_results6, SubID~State)

#Rename variables
library(plyr)
temp_probe_results1 <- plyr::rename(temp_probe_results1, c("1"="ON_Freq", "2"="MW_Freq", "3"="MB_Freq", "4"="DR"))
temp_probe_results2 <- plyr::rename(temp_probe_results2, c("1"="ON_Awa", "2"="MW_Awa", "3"="MB_Awa", "4"="DR"))
temp_probe_results3 <- plyr::rename(temp_probe_results3, c("1"="ON_Int", "2"="MW_Int", "3"="MB_Int", "4"="DR"))
temp_probe_results4 <- plyr::rename(temp_probe_results4, c("1"="ON_Eng", "2"="MW_Eng", "3"="MB_Eng", "4"="DR"))
temp_probe_results5 <- plyr::rename(temp_probe_results5, c("1"="ON_Perf", "2"="MW_Perf", "3"="MB_Perf", "4"="DR"))
temp_probe_results6 <- plyr::rename(temp_probe_results6, c("1"="ON_Vig", "2"="MW_Vig", "3"="MB_Vig", "4"="DR"))

# Drop don't remember probes
#library(dplyr)
#temp_probe_results1 <- select(temp_probe_results1, -DR)
#temp_probe_results2 <- select(temp_probe_results2, -DR)
#temp_probe_results3 <- select(temp_probe_results3, -DR)
#temp_probe_results4 <- select(temp_probe_results4, -DR)
#temp_probe_results5 <- select(temp_probe_results5, -DR)
#temp_probe_results6 <- select(temp_probe_results6, -DR)

# Merge tables
Probe_Results_OtQs <- merge(temp_probe_results2, temp_probe_results3, by=c("SubID","State"))
Probe_Results_OtQs <- merge(Probe_Results_OtQs, temp_probe_results4, by=c("SubID","State"))
Probe_Results_OtQs <- merge(Probe_Results_OtQs, temp_probe_results5, by=c("SubID","State"))
Probe_Results_OtQs <- merge(Probe_Results_OtQs, temp_probe_results6, by=c("SubID","State"))

#Remove Don't know responses
Probe_Results_OtQs$State <- mapvalues(Probe_Results_OtQs$State, from = c(1,2,3,4), to = c("ON","MW","MB","DR"))
Probe_Results_OtQs <- subset(Probe_Results_OtQs$State, select=c("ON","MW","MB"))
Probe_Results_OtQs <- subset(Probe_Results_OtQs, State!="DR")

krusk_Awa <- kruskal.test(Awa~State, data=Probe_Results_OtQs)
krusk_Awa
dunnTest(Awa~State, data=Probe_Results_OtQs)

krusk_Int <- kruskal.test(Int~State, data=Probe_Results_OtQs)
krusk_Int
dunnTest(Int~State, data=Probe_Results_OtQs)

krusk_Eng <- kruskal.test(Eng~State, data=Probe_Results_OtQs)
krusk_Eng
dunnTest(Eng~State, data=Probe_Results_OtQs)

krusk_Perf <- kruskal.test(Perf~State, data=Probe_Results_OtQs)
krusk_Perf
dunnTest(Perf~State, data=Probe_Results_OtQs)

krusk_Vig <- kruskal.test(Vig~State, data=Probe_Results_OtQs)
krusk_Vig
dunnTest(Vig~State, data=Probe_Results_OtQs)


#### ARCHIVED / NOT USED ####

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



# CORRECTNESS ACROSS NOGO TRIALS by NBLOCK

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


# CORRECTNESS ACROSS GO TRIALS by NBLOCK #

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





# REACTION TIME ACROSS GO TRIALS by NBLOCK

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







# CORRECTNESS BY TASK #

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





# REACTION TIME BY TASK #
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




# OVERALL CORRECTNESS BY TRIAL N WITHIN BLOCK - issue viewing this - too few instances to make sense #

tbl_perf_by_trial_ind <- table(WanderIM_TestResults$nTrial, WanderIM_TestResults$Corr)
wrongs <- tbl_perf_by_trial_ind[1:474, 1]
rights <- tbl_perf_by_trial_ind[1:474, 2]
tbl_perf_by_trial_ind <- rights / (rights+wrongs)

remove(tbl_perf_by_trial_ind)

tbl_perf_by_trial_ind <- melt(tbl_perf_by_trial_ind, id.vars = c("nTrial"))
names(tbl_perf_by_trial_ind)[names(tbl_perf_by_trial_ind)=="value"] <- "perc_correct"

plot(tbl_perf_by_trial_ind)


# REACTION TIME BY TRIAL N WITHIN BLOCK #
#Task coding: 1=faces, 2=squares - taken from Thomas' comments in matlab script

#Reshape data to view response time by individual participant across trials
tbl_RT_by_trial_ind_go <- subset(WanderIM_ProbeResults, TrCat == "0")

tbl_RT_by_trial_ind_go <- aggregate(tbl_RT_by_trial_ind_go$RT~tbl_RT_by_trial_ind_go$nTrial+tbl_RT_by_trial_ind_go$SubID, FUN="mean")

plot(tbl_RT_by_trial_ind_go$`tbl_RT_by_trial_ind_go$nTrial`, tbl_RT_by_trial_ind_go$`tbl_RT_by_trial_ind_go$RT`)
abline(lm(tbl_RT_by_trial_ind_go$`tbl_RT_by_trial_ind_go$RT`~tbl_RT_by_trial_ind_go$`tbl_RT_by_trial_ind_go$nTrial`), col="red") # regression line (y~x) 
lines(lowess(tbl_RT_by_trial_ind_go$`tbl_RT_by_trial_ind_go$nTrial`, tbl_RT_by_trial_ind_go$`tbl_RT_by_trial_ind_go$RT`), col="blue") # lowess line (x,y)



# OVERALL CORRECTNESS MODELLING - ALL GO/NOGO TRIALS COMBINED (COMPLEX) #
# Modelling the task parameters
library(lme4)

## Model: correctness ~ 100% (baseline)
mdl_0 <- lmer(Corr ~ 1 + (1| SubID), data=WanderIM_ProbeResults)
## Model: correctness ~ Trial category (Go / No-go)
mdl_1 <- lmer(Corr ~ TrCat + (1| SubID), data=WanderIM_ProbeResults)

    #### Model comparison: Mdl_1 wins
    anova(mdl_0,mdl_1)

## Model: correctness ~ trial category and task (face/digit)
mdl_2 <- lmer(Corr ~ TrCat+Task + (1| SubID), data=WanderIM_ProbeResults)

    #### Model comparison: Mdl_2 wins
    anova(mdl_1,mdl_2)

## Model: correctness ~ trial category and interaction effect with task (face/digit)
mdl_3 <- lmer(Corr ~ TrCat*Task + (1| SubID), data=WanderIM_ProbeResults)

    #### Model comparison: Mdl_3 wins
    anova(mdl_2,mdl_3)

## Model: correctness ~ trial category and interaction effect with task (face/digit) and block number
mdl_4 <- lmer(Corr ~ TrCat*Task+nBlock + (1| SubID), data=WanderIM_ProbeResults)

    #### Model comparison: Mdl_4 wins
    anova(mdl_3,mdl_4)

## Model: correctness ~ trial category and interaction effect with task (face/digit) and interaction effect with block number
mdl_5 <- lmer(Corr ~ TrCat*Task*nBlock + (1| SubID), data=WanderIM_ProbeResults)

    #### Model comparison: Mdl_5 wins
    anova(mdl_4,mdl_5)

## Model: correctness ~ trial category and interaction effect with task (face/digit) and interaction effect with block number with effect of state
mdl_6 <- lmer(Corr ~ TrCat*Task*nBlock+Look + (1| SubID), data=WanderIM_ProbeResults)

    #### Model comparison: Mdl_6 wins
    anova(mdl_5,mdl_6)    

## Model: correctness ~ trial category and interaction effect with task (face/digit) and interaction effect with block number and interaction effect of state
mdl_7 <- lmer(Corr ~ TrCat*Task*nBlock*Look + (1| SubID), data=WanderIM_ProbeResults)

    #### Model comparison: Mdl_7 wins
    anova(mdl_6,mdl_7)    

## Model: 
mdl_8 <- lmer(Corr ~ TrCat*Task*nBlock*Look+Perf + (1| SubID), data=WanderIM_ProbeResults)
    
    #### Model comparison: Mdl_8 wins
    anova(mdl_7,mdl_8) 
       
## Model: 
mdl_9 <- lmer(Corr ~ TrCat*Task*nBlock*Look*Perf + (1| SubID), data=WanderIM_ProbeResults)
    
    #### Model comparison: Mdl_9 wins
    anova(mdl_8,mdl_9) 
    
## Model: 
mdl_10 <- lmer(Corr ~ TrCat*Task*nBlock*Look*Perf+Vig + (1| SubID), data=WanderIM_ProbeResults)
    
    #### Model comparison: Mdl_9 wins
    anova(mdl_9,mdl_10) 
    
## Model: 
mdl_11 <- lmer(Corr ~ TrCat*Task*nBlock*Look*Perf*Vig + (1| SubID), data=WanderIM_ProbeResults)
    
    #### Model comparison: Mdl_11 wins
    anova(mdl_9,mdl_11) 
    
## Model: 
mdl_12 <- lmer(Corr ~ TrCat*Task*nBlock*Look*Perf*Vig+State + (1| SubID), data=WanderIM_ProbeResults)
    
    #### Model comparison: Mdl_12 wins
    anova(mdl_11,mdl_12) 
    
## Model: 
mdl_13 <- lmer(Corr ~ TrCat*Task*nBlock*Look*Perf*Vig*State + (1| SubID), data=WanderIM_ProbeResults)
    
    #### Model comparison: Mdl_13 wins
    anova(mdl_12,mdl_13) 
    
#inspect winning model
summary(mdl_13)
anova(mdl_13)

# Create SEM Plot -- not working. Can we do SEM/Path Analysis on a ML model?
model <- WanderIM_ProbeResults$Corr ~ WanderIM_ProbeResults$TrCat*WanderIM_ProbeResults$Task*WanderIM_ProbeResults$nBlock*WanderIM_ProbeResults$State*WanderIM_ProbeResults$Awa*WanderIM_ProbeResults$Int*WanderIM_ProbeResults$Eng*WanderIM_ProbeResults$Vig + (1| WanderIM_ProbeResults$SubID)

fit <- cfa(model, data = WanderIM_ProbeResults)

lavaan.diagram()
lavaan(model = mdl_15, data = WanderIM_ProbeResults)

# REACTION TIME MODELLING #
# Modelling the task parameters

## Model: Reaction Time ~ 100% (baseline)
RTmdl_0 <- lmer(RT ~ 1 + (1| SubID), data=WanderIM_ProbeResults)
## Model: Reaction Time ~ Trial category (Go / No-go)
RTmdl_1 <- lmer(RT ~ TrCat*Task*nBlock*Look*Perf*Vig*State + (1| SubID), data=WanderIM_ProbeResults)

  #### Model comparison: Mdl_1 wins
  anova(RTmdl_0,RTmdl_1)

## Model: Reaction Time ~ 
RTmdl_2 <- lmer(RT ~ TrCat+Task + (1| SubID), data=WanderIM_ProbeResults)
  
  #### Model comparison: Mdl_2 wins
  anova(RTmdl_1,RTmdl_2)

## Model: Reaction Time ~ 
RTmdl_3 <- lmer(RT ~ TrCat*Task + (1| SubID), data=WanderIM_ProbeResults)
  
  #### Model comparison: Mdl_3 wins
  anova(RTmdl_2,RTmdl_3)

## Model: Reaction Time ~ 
RTmdl_4 <- lmer(RT ~ TrCat*Task+nBlock + (1| SubID), data=WanderIM_ProbeResults)
  
  #### Model comparison: Mdl_4 wins
  anova(RTmdl_3,RTmdl_4)

## Model: Reaction Time ~ 
RTmdl_5 <- lmer(RT ~ TrCat*Task*nBlock + (1| SubID), data=WanderIM_ProbeResults)
  
  #### Model comparison: Mdl_5 wins
  anova(RTmdl_4,RTmdl_5)

## Model: Reaction Time ~ 
RTmdl_6 <- lmer(RT ~ TrCat*Task*nBlock+Look + (1| SubID), data=WanderIM_ProbeResults)
  
  #### Model comparison: Mdl_6 wins
  anova(RTmdl_5,RTmdl_6)

## Model: Reaction Time ~ 
RTmdl_7 <- lmer(RT ~ TrCat*Task*nBlock*Look + (1| SubID), data=WanderIM_ProbeResults)

  #### Model comparison: Mdl_7 wins
  anova(RTmdl_6,RTmdl_7)

## Model: Reaction Time ~ 
RTmdl_8 <- lmer(RT ~ TrCat*Task*nBlock*Look+Perf + (1| SubID), data=WanderIM_ProbeResults)

## Model: Reaction Time ~ 
RTmdl_9 <- lmer(RT ~ TrCat*Task*nBlock*Look*Perf + (1| SubID), data=WanderIM_ProbeResults)

## Model: Reaction Time ~ 
RTmdl_10 <- lmer(RT ~ TrCat*Task*nBlock*Look*Perf+Vig + (1| SubID), data=WanderIM_ProbeResults)

## Model: Reaction Time ~ 
RTmdl_11 <- lmer(RT ~ TrCat*Task*nBlock*Look*Perf*Vig + (1| SubID), data=WanderIM_ProbeResults)

## Model: Reaction Time ~ 
RTmdl_12 <- lmer(RT ~ TrCat*Task*nBlock*Look*Perf*Vig+State + (1| SubID), data=WanderIM_ProbeResults)

## Model: Reaction Time ~ 
RTmdl_13 <- lmer(RT ~ TrCat*Task*nBlock*Look*Perf*Vig*State + (1| SubID), data=WanderIM_ProbeResults)











#### OVERALL CORRECTNESS MODELLING - NOGO TRIALS ONLY ####
# Modelling the task parameters
library(lme4)

WIM_NOGO <- subset(WanderIM_ProbeResults, TrCat == 0)

## Model: correctness ~ 100% (baseline)
nogo_mdl_0 <- lmer(Corr ~ 1 + (1| SubID), data=WIM_NOGO)
## Model: correctness ~ Trial category (Go / No-go)
nogo_mdl_1 <- lmer(Corr ~ Task + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogo_mdl_1 wins
anova(nogo_mdl_0,nogo_mdl_1)

## Model: correctness ~ trial category and task (face/digit)
nogo_mdl_2 <- lmer(Corr ~ Task+nBlock + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogo_mdl_2 wins
anova(nogo_mdl_1,nogo_mdl_2)

## Model: correctness ~ trial category and interaction effect with task (face/digit)
nogo_mdl_3 <- lmer(Corr ~ Task*nBlock + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogo_mdl_3 wins
anova(nogo_mdl_2,nogo_mdl_3)

## Model: correctness ~ trial category and interaction effect with task (face/digit) and block number
nogo_mdl_4 <- lmer(Corr ~ Task*nBlock+Perf + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogo_mdl_4 wins
anova(nogo_mdl_3,nogo_mdl_4)

## Model: correctness ~ trial category and interaction effect with task (face/digit) and interaction effect with block number
nogo_mdl_5 <- lmer(Corr ~ Task*nBlock*Perf + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogo_mdl_5 wins
anova(nogo_mdl_4,nogo_mdl_5)

## Model: correctness ~ trial category and interaction effect with task (face/digit) and interaction effect with block number with effect of state
nogo_mdl_6 <- lmer(Corr ~ Task*nBlock*Perf+Vig + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogo_mdl_6 wins
anova(nogo_mdl_5,nogo_mdl_6)    

## Model: correctness ~ trial category and interaction effect with task (face/digit) and interaction effect with block number and interaction effect of state
nogo_mdl_7 <- lmer(Corr ~ Task*nBlock*Perf*Vig + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogo_mdl_7 wins
anova(nogo_mdl_6,nogo_mdl_7)    

## Model: 
nogo_mdl_8 <- lmer(Corr ~ Task*nBlock*Perf*Vig+State + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogo_mdl_8 wins
anova(nogo_mdl_7,nogo_mdl_8) 

## Model: 
nogo_mdl_9 <- lmer(Corr ~ Task*nBlock*Perf*Vig*State + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogo_mdl_9 wins
anova(nogo_mdl_8,nogo_mdl_9) 


#inspect winning model
summary(nogo_mdl_9)
anova(nogo_mdl_9)






#### OVERALL CORRECTNESS MODELLING - GO TRIALS ONLY ####
# Modelling the task parameters
library(lme4)

WIM_GO <- subset(WanderIM_ProbeResults, TrCat == 1)

## Model: correctness ~ 100% (baseline)
go_mdl_0 <- lmer(Corr ~ 1 + (1| SubID), data=WIM_GO)
## Model: correctness ~ Trial category (Go / No-go)
go_mdl_1 <- lmer(Corr ~ Task + (1| SubID), data=WIM_GO)

#### Model comparison: go_mdl_1 wins
anova(go_mdl_0,go_mdl_1)

## Model: correctness ~ trial category and task (face/digit)
go_mdl_2 <- lmer(Corr ~ Task+nBlock + (1| SubID), data=WIM_GO)

#### Model comparison: go_mdl_2 wins
anova(go_mdl_1,go_mdl_2)

## Model: correctness ~ trial category and interaction effect with task (face/digit)
go_mdl_3 <- lmer(Corr ~ Task*nBlock + (1| SubID), data=WIM_GO)

#### Model comparison: go_mdl_2 wins
anova(go_mdl_2,go_mdl_3)

## Model: correctness ~ trial category and interaction effect with task (face/digit) and block number
go_mdl_4 <- lmer(Corr ~ Task*nBlock+Perf + (1| SubID), data=WIM_GO)

#### Model comparison: go_mdl_4 wins
anova(go_mdl_2,go_mdl_4)

## Model: correctness ~ trial category and interaction effect with task (face/digit) and interaction effect with block number
go_mdl_5 <- lmer(Corr ~ Task*nBlock*Perf + (1| SubID), data=WIM_GO)

#### Model comparison: go_mdl_5 wins
anova(go_mdl_4,go_mdl_5)

## Model: correctness ~ trial category and interaction effect with task (face/digit) and interaction effect with block number with effect of state
go_mdl_6 <- lmer(Corr ~ Task*nBlock*Perf+Vig + (1| SubID), data=WIM_GO)

#### Model comparison: go_mdl_5 wins
anova(go_mdl_5,go_mdl_6)    

## Model: correctness ~ trial category and interaction effect with task (face/digit) and interaction effect with block number and interaction effect of state
go_mdl_7 <- lmer(Corr ~ Task*nBlock*Perf*Vig + (1| SubID), data=WIM_GO)

#### Model comparison: go_mdl_5 wins
anova(go_mdl_5,go_mdl_7)    

## Model: 
go_mdl_8 <- lmer(Corr ~ Task*nBlock*Perf+State + (1| SubID), data=WIM_GO)

#### Model comparison: go_mdl_8 wins
anova(go_mdl_5,go_mdl_8) 

## Model: 
go_mdl_9 <- lmer(Corr ~ Task*nBlock*Perf*State + (1| SubID), data=WIM_GO)

#### Model comparison: go_mdl_8 wins
anova(go_mdl_8,go_mdl_9) 


#inspect winning model
summary(go_mdl_8)
anova(go_mdl_8)








#### REACTION TIME MODELLING - NOGO TRIALS ONLY ####
# Modelling the task parameters
library(lme4)

WIM_NOGO <- subset(WanderIM_ProbeResults, TrCat == 0)

## Model: RT ~ 100% (baseline)
nogoRT_mdl_0 <- lmer(RT ~ 1 + (1| SubID), data=WIM_NOGO)
## Model: 
nogoRT_mdl_1 <- lmer(RT ~ Task + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogoRT_mdl_1 wins
anova(nogoRT_mdl_0,nogoRT_mdl_1)

## Model: 
nogoRT_mdl_2 <- lmer(RT ~ Task+nBlock + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogoRT_mdl_1 wins
anova(nogoRT_mdl_1,nogoRT_mdl_2)

## Model: 
nogoRT_mdl_3 <- lmer(RT ~ Task*nBlock + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogoRT_mdl_3 wins
anova(nogoRT_mdl_1,nogoRT_mdl_3)

## Model: 
nogoRT_mdl_4 <- lmer(RT ~ Task*nBlock+Perf + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogoRT_mdl_4 wins
anova(nogoRT_mdl_3,nogoRT_mdl_4)

## Model: 
nogoRT_mdl_5 <- lmer(RT ~ Task*nBlock*Perf + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogoRT_mdl_4 wins
anova(nogoRT_mdl_4,nogoRT_mdl_5)

## Model: 
nogoRT_mdl_6 <- lmer(RT ~ Task*nBlock+Perf+Vig + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogoRT_mdl_6 wins
anova(nogoRT_mdl_4,nogoRT_mdl_6)    

## Model: 
nogoRT_mdl_7 <- lmer(RT ~ Task*nBlock+Perf*Vig + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogoRT_mdl_7 wins
anova(nogoRT_mdl_6,nogoRT_mdl_7)    

## Model: 
nogoRT_mdl_8 <- lmer(RT ~ Task*nBlock+Perf*Vig+State + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogoRT_mdl_8 wins
anova(nogoRT_mdl_7,nogoRT_mdl_8) 

## Model: 
nogoRT_mdl_9 <- lmer(RT ~ Task*nBlock+Perf*Vig*State + (1| SubID), data=WIM_NOGO)

#### Model comparison: nogoRT_mdl_9 wins
anova(nogoRT_mdl_8,nogoRT_mdl_9) 


#inspect winning model
summary(nogoRT_mdl_9)
anova(nogoRT_mdl_9)






#### REACTION TIME MODELLING - GO TRIALS ONLY ####
# Modelling the task parameters
library(lme4)

WIM_GO <- subset(WanderIM_ProbeResults, TrCat == 1)

## Model: RT ~ 100% (baseline)
goRT_mdl_0 <- lmer(RT ~ 1 + (1| SubID), data=WIM_GO)
## Model: 
goRT_mdl_1 <- lmer(RT ~ Task + (1| SubID), data=WIM_GO)

#### Model comparison: goRT_mdl_1 wins
anova(goRT_mdl_0,goRT_mdl_1)

## Model: 
goRT_mdl_2 <- lmer(RT ~ Task+nBlock + (1| SubID), data=WIM_GO)

#### Model comparison: goRT_mdl_2 wins
anova(goRT_mdl_1,goRT_mdl_2)

## Model: 
goRT_mdl_3 <- lmer(RT ~ Task*nBlock + (1| SubID), data=WIM_GO)

#### Model comparison: goRT_mdl_2 wins
anova(goRT_mdl_2,goRT_mdl_3)

## Model: 
goRT_mdl_4 <- lmer(RT ~ Task+nBlock+Perf + (1| SubID), data=WIM_GO)

#### Model comparison: goRT_mdl_4 wins
anova(goRT_mdl_2,goRT_mdl_4)

## Model: 
goRT_mdl_5 <- lmer(RT ~ Task+nBlock*Perf + (1| SubID), data=WIM_GO)

#### Model comparison: goRT_mdl_4 wins
anova(goRT_mdl_4,goRT_mdl_5)

## Model: 
goRT_mdl_6 <- lmer(RT ~ Task+nBlock+Perf+Vig + (1| SubID), data=WIM_GO)

#### Model comparison: goRT_mdl_4 wins
anova(goRT_mdl_4,goRT_mdl_6)    

## Model: 
goRT_mdl_7 <- lmer(RT ~ Task+nBlock+Perf*Vig + (1| SubID), data=WIM_GO)

#### Model comparison: goRT_mdl_4 wins
anova(goRT_mdl_4,goRT_mdl_7)    

## Model: 
goRT_mdl_8 <- lmer(RT ~ Task+nBlock+Perf+State + (1| SubID), data=WIM_GO)

#### Model comparison: goRT_mdl_4 wins
anova(goRT_mdl_4,goRT_mdl_8) 

## Model: 
goRT_mdl_9 <- lmer(RT ~ Task+nBlock+Perf*State + (1| SubID), data=WIM_GO)

#### Model comparison: goRT_mdl_4 wins
anova(goRT_mdl_4,goRT_mdl_9) 


#inspect winning model
summary(goRT_mdl_4)
anova(goRT_mdl_4)











#### OVERALL CORRECTNESS BY PROBE (20 TRIALS BEFORE PROBE) - finish lme analysis ####

#Remove all trials >-20 away from probe
tbl_perf_by_probe <- subset(WanderIM_ProbeResults, DistProbe >-21)

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



#### REACTION TIME BY PROBE (20 TRIALS BEFORE PROBE) - finish lme analysis ####

#Remove all trials >-20 away from probe
tbl_RT_by_probe <- subset(WanderIM_ProbeResults, DistProbe >-21)

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

#### ANOVA OF SURVEY RESPONSES AGAINST MW/MB/ON FREQUENCY & Performance - PREP DATA ####

## BEHAVIOURAL EXPERIMENT DATA
#Remove all trials >-20 away from probe
tbl_perf_by_probe <- subset(WanderIM_ProbeResults, DistProbe >-21)

# Aggregate to correctness & RT, including relevant variables
attach(tbl_perf_by_probe)
tbl_Correctness_ID_State <- aggregate(Corr~SubID+State, FUN="mean")
tbl_RT_ID_State <- aggregate(RT~SubID+State, FUN="mean")
detach(tbl_perf_by_probe)

#Reshape so attention state is own variable
library(reshape2)
tbl_Correctness_ID_State <- dcast(tbl_Correctness_ID_State, SubID~State)
tbl_RT_ID_State <- dcast(tbl_RT_ID_State, SubID~State)
#Rename variables
library(plyr)
tbl_Correctness_ID_State <- plyr::rename(tbl_Correctness_ID_State, c("1"="ON_Correctness", "2"="MW_Correctness", "3"="MB_Correctness", "4"="DR"))
tbl_RT_ID_State <- plyr::rename(tbl_RT_ID_State, c("1"="ON_RT", "2"="MW_RT", "3"="MB_RT", "4"="DR"))
# Drop don't remember probes
library(dplyr)
tbl_Correctness_ID_State <- select(tbl_Correctness_ID_State, -DR)
tbl_RT_ID_State <- select(tbl_RT_ID_State, -DR)



## Frequency of attention states
# Aggregate to Frequency of MB/MW/ON
tbl_freq_probe <- subset(WanderIM_ProbeResults, DistProbe >-2)
attach(tbl_freq_probe)
tbl_freq_probe <- aggregate(nTrial~SubID+State, FUN = "length")
detach(tbl_freq_probe)
#Pivot so attention each attention state is own variable
library(reshape2)
tbl_freq_probe_temp <- dcast(tbl_freq_probe, SubID~State)
#Rename variables
library(plyr)
tbl_freq_probe_temp <- plyr::rename(tbl_freq_probe_temp, c("1"="ON_Freq", "2"="MW_Freq", "3"="MB_Freq", "4"="DR"))
# Drop don't remember probes
library(dplyr)
tbl_freq_probe_temp <- select(tbl_freq_probe_temp, -DR)
# Calculate percentage values for each of the mind-states (so the three mind-states of interest add up to 100%)
SUM_Probes = rowSums(tbl_freq_probe_temp[,c("ON_Freq","MW_Freq","MB_Freq")], na.rm = TRUE)
tbl_freq_probe_temp$ON_perc <- tbl_freq_probe_temp$ON_Freq / SUM_Probes
tbl_freq_probe_temp$MW_perc <- tbl_freq_probe_temp$MW_Freq / SUM_Probes
tbl_freq_probe_temp$MB_perc <- tbl_freq_probe_temp$MB_Freq / SUM_Probes



#Join survey and outcome variables
tbl_perf_ind <- merge(tbl_Correctness_ID_State, tbl_RT_ID_State, by="SubID")
tbl_perf_ind <- merge(tbl_perf_ind, tbl_freq_probe_temp, by="SubID")
tbl_perf_ind_merged <- merge(tbl_perf_ind,Survey.df, by.x="SubID", by.y="Participant_No")
tbl_perf_ind_merged <- merge(tbl_perf_ind_merged,NAB.df, by.x="SubID", by.y="Participant_No")
tbl_perf_ind_merged <- merge(tbl_perf_ind_merged,Guilford_AUT, by.x="SubID", by.y="ID")
tbl_perf_ind_merged <- merge(tbl_perf_ind_merged,behav_resp_DV, by.x="SubID", by.y="Sub")
tbl_perf_ind_merged <- merge(tbl_perf_ind_merged,OSPAN.new, by.x="SubID", by.y="script.subjectid")


attach(OSPAN.dat)
myvars <- c("script.subjectid", "values.ospan")
OSPAN.new <- OSPAN.dat[myvars]
detach(OSPAN.dat)

tbl_perf_ind_merged <- merge(tbl_perf_ind_merged,OSPAN.new, by.x="SubID", by.y="script.subjectid")

# Turn attention state frequencies into categories
# ON - Find quartiles
summary(tbl_perf_ind_merged$ON_perc)
# Recode
tbl_perf_ind_merged$ON_CAT[is.na(tbl_perf_ind_merged$ON_perc)] <- 'None'
tbl_perf_ind_merged$ON_CAT[tbl_perf_ind_merged$ON_perc <=0.33190] <- 'Low'
tbl_perf_ind_merged$ON_CAT[tbl_perf_ind_merged$ON_perc >=0.33190 & tbl_perf_ind_merged$ON_perc <=.43333] <- 'Avg-Low'
tbl_perf_ind_merged$ON_CAT[tbl_perf_ind_merged$ON_perc >=.43333 & tbl_perf_ind_merged$ON_perc <=.65417] <- 'Avg-High'
tbl_perf_ind_merged$ON_CAT[tbl_perf_ind_merged$ON_perc >=0.65417 & tbl_perf_ind_merged$ON_perc <=1] <- 'High'
# Check resulting data
tbl_perf_ind_merged$ON_CAT

# MW - Find quartiles
summary(tbl_perf_ind_merged$MW_perc)
# Recode
tbl_perf_ind_merged$MW_CAT[is.na(tbl_perf_ind_merged$MW_perc)] <- 'None'
tbl_perf_ind_merged$MW_CAT[tbl_perf_ind_merged$MW_perc <=0.2499] <- 'Low'
tbl_perf_ind_merged$MW_CAT[tbl_perf_ind_merged$MW_perc >=0.2499 & tbl_perf_ind_merged$MW_perc <=.3333] <- 'Avg-Low'
tbl_perf_ind_merged$MW_CAT[tbl_perf_ind_merged$MW_perc >=.3333 & tbl_perf_ind_merged$MW_perc <=.5359] <- 'Avg-High'
tbl_perf_ind_merged$MW_CAT[tbl_perf_ind_merged$MW_perc >=0.5359 & tbl_perf_ind_merged$MW_perc <=1] <- 'High'
# Check resulting data
tbl_perf_ind_merged$MW_CAT

# MB - Find quartiles
summary(tbl_perf_ind_merged$MB_perc)
# Recode
tbl_perf_ind_merged$MB_CAT[is.na(tbl_perf_ind_merged$MB_perc)] <- 'None'
tbl_perf_ind_merged$MB_CAT[tbl_perf_ind_merged$MB_perc <=0.08333] <- 'Low'
tbl_perf_ind_merged$MB_CAT[tbl_perf_ind_merged$MB_perc >=0.08333 & tbl_perf_ind_merged$MB_perc <=.13333] <- 'Avg-Low'
tbl_perf_ind_merged$MB_CAT[tbl_perf_ind_merged$MB_perc >=.13333 & tbl_perf_ind_merged$MB_perc <=.27624] <- 'Avg-High'
tbl_perf_ind_merged$MB_CAT[tbl_perf_ind_merged$MB_perc >=0.27624 & tbl_perf_ind_merged$MB_perc <=.40678] <- 'High'
# Check resulting data
tbl_perf_ind_merged$MB_CAT


# Specify factors
tbl_perf_ind_merged$ON_CAT <- as.factor(tbl_perf_ind_merged$ON_CAT)
tbl_perf_ind_merged$MW_CAT <- as.factor(tbl_perf_ind_merged$MW_CAT)
tbl_perf_ind_merged$MB_CAT <- as.factor(tbl_perf_ind_merged$MB_CAT)
tbl_perf_ind_merged$EOY17 <- as.factor(tbl_perf_ind_merged$EOY2017)


#### ANALYSIS OF SURVEY VARIABLES WITH ATTENTION STATES AND EMPLOYEE SCORES AS OUTCOMES ####

#DASS - Numeric
krusk_DepON <- kruskal.test(DASS_Depression ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_DepMW <- kruskal.test(DASS_Depression ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_DepMB <- kruskal.test(DASS_Depression ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_DepEOY17 <- kruskal.test(DASS_Depression ~ EOY2017, data = tbl_perf_ind_merged)

krusk_AnxON <- kruskal.test(DASS_Anxiety ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_AnxMW <- kruskal.test(DASS_Anxiety ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_AnxMB <- kruskal.test(DASS_Anxiety ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_EOY17 <- kruskal.test(DASS_Anxiety ~ EOY2017, data = tbl_perf_ind_merged)

krusk_StrON <- kruskal.test(DASS_Stress ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_StrMW <- kruskal.test(DASS_Stress ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_StrMB <- kruskal.test(DASS_Stress ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_StrEOY17 <- kruskal.test(DASS_Stress ~ EOY2017, data = tbl_perf_ind_merged)

#ESS - Numeric
krusk_ESSON <- kruskal.test(ESS ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_ESSMW <- kruskal.test(ESS ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_ESSMB <- kruskal.test(ESS ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_ESSEOY17 <- kruskal.test(ESS ~ EOY2017, data = tbl_perf_ind_merged)

#MWQ - Numeric
krusk_MWQON <- kruskal.test(MWQ_mean ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_MWQMW <- kruskal.test(MWQ_mean ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_MWQMB <- kruskal.test(MWQ_mean ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_MWQEOY17 <- kruskal.test(MWQ_mean ~ EOY2017, data = tbl_perf_ind_merged)

# Mini IPIP - Numeric
krusk_IPNeuON <- kruskal.test(MiniIPIP_Neuroticism ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_IPNeuMW <- kruskal.test(MiniIPIP_Neuroticism ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_IPNeuMB <- kruskal.test(MiniIPIP_Neuroticism ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_IPNeuEOY17 <- kruskal.test(MiniIPIP_Neuroticism ~ EOY2017, data = tbl_perf_ind_merged)

krusk_IPExtON <- kruskal.test(MiniIPIP_Extraversion ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_IPExtMW <- kruskal.test(MiniIPIP_Extraversion ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_IPExtMB <- kruskal.test(MiniIPIP_Extraversion ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_IPExtEOY17 <- kruskal.test(MiniIPIP_Extraversion ~ EOY2017, data = tbl_perf_ind_merged)

krusk_IPOpeON <- kruskal.test(MiniIPIP_Openness ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_IPOpeMW <- kruskal.test(MiniIPIP_Openness ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_IPOpeMB <- kruskal.test(MiniIPIP_Openness ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_IPOpeEOY17 <- kruskal.test(MiniIPIP_Openness ~ EOY2017, data = tbl_perf_ind_merged)

krusk_IPAgrON <- kruskal.test(MiniIPIP_Agreeableness ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_IPAgrMW <- kruskal.test(MiniIPIP_Agreeableness ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_IPAgrMB <- kruskal.test(MiniIPIP_Agreeableness ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_IPAgrEOY17 <- kruskal.test(MiniIPIP_Agreeableness ~ EOY2017, data = tbl_perf_ind_merged)

krusk_IPConON <- kruskal.test(MiniIPIP_Conscientiousness ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_IPConMW <- kruskal.test(MiniIPIP_Conscientiousness ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_IPConMB <- kruskal.test(MiniIPIP_Conscientiousness ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_IPConEOY17 <- kruskal.test(MiniIPIP_Conscientiousness ~ EOY2017, data = tbl_perf_ind_merged)

# PANAS - Numeric
krusk_PANPosON <- kruskal.test(PANAS_Positive_Affect ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_PANPosMW <- kruskal.test(PANAS_Positive_Affect ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_PANPosMB <- kruskal.test(PANAS_Positive_Affect ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_PANPosEOY17 <- kruskal.test(PANAS_Positive_Affect ~ EOY2017, data = tbl_perf_ind_merged)

krusk_PANNegON <- kruskal.test(PANAS_Negative_Affect ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_PANNegMW <- kruskal.test(PANAS_Negative_Affect ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_PANNegMB <- kruskal.test(PANAS_Negative_Affect ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_PANNegEOY17 <- kruskal.test(PANAS_Negative_Affect ~ EOY2017, data = tbl_perf_ind_merged)

#ASRS - Categorical
chisq_ASRSON <- chisq.test(tbl_perf_ind_merged$ASRS_CAT, tbl_perf_ind_merged$ON_CAT)
chisq_ASRSMW <- chisq.test(tbl_perf_ind_merged$ASRS_CAT, tbl_perf_ind_merged$MW_CAT)
chisq_ASRSMB <- chisq.test(tbl_perf_ind_merged$ASRS_CAT, tbl_perf_ind_merged$MB_CAT)
chisq_ASRSEOY17 <- chisq.test(tbl_perf_ind_merged$ASRS_CAT, tbl_perf_ind_merged$EOY2017)

#ASRS - Numeric
krusk_ASRSON <- kruskal.test(ASRS_ptA_num ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_ASRSMW <- kruskal.test(ASRS_ptA_num ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_ASRSMB <- kruskal.test(ASRS_ptA_num ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_ASRSEOY17 <- kruskal.test(ASRS_ptA_num ~ EOY2017, data = tbl_perf_ind_merged)

#BRS - Numeric
krusk_BRSON <- kruskal.test(BRS_mean ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_BRSMW <- kruskal.test(BRS_mean ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_BRSMB <- kruskal.test(BRS_mean ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_BRSEOY17 <- kruskal.test(BRS_mean ~ EOY2017, data = tbl_perf_ind_merged)

#WDQ - Problem Solving - Numeric
krusk_WDQProbsON <- kruskal.test(WDQ_Problem_Solving ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_WDQProbsMW <- kruskal.test(WDQ_Problem_Solving ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_WDQProbsMB <- kruskal.test(WDQ_Problem_Solving ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_WDQProbsEOY17 <- kruskal.test(WDQ_Problem_Solving ~ EOY2017, data = tbl_perf_ind_merged)

#Guilford's AUT - Numeric
krusk_AUTON <- kruskal.test(AUT_Score ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_AUTMW <- kruskal.test(AUT_Score ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_AUTMB <- kruskal.test(AUT_Score ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_AUTEOY17 <- kruskal.test(AUT_Score ~ EOY2017, data = tbl_perf_ind_merged)

#EES - Employee Engagement - Numeric
krusk_EESEngON <- kruskal.test(EES_Employee_Engagement ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_EESEngMW <- kruskal.test(EES_Employee_Engagement ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_EESEngMB <- kruskal.test(EES_Employee_Engagement ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_EESEngEOY17 <- kruskal.test(EES_Employee_Engagement ~ EOY2017, data = tbl_perf_ind_merged)

#EOY Performance Rating - Categorical
chisq_EOY17ON <- chisq.test(tbl_perf_ind_merged$EOY2017, tbl_perf_ind_merged$ON_CAT)
chisq_EOY17MW <- chisq.test(tbl_perf_ind_merged$EOY2017, tbl_perf_ind_merged$MW_CAT)
chisq_EOY17MB <- chisq.test(tbl_perf_ind_merged$EOY2017, tbl_perf_ind_merged$MB_CAT)

# Task performance % correct
krusk_ONPerfON <- kruskal.test(ON_Correctness ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_ONPerfMW <- kruskal.test(ON_Correctness ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_ONPerfMB <- kruskal.test(ON_Correctness ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_ONPerfEOY17 <- kruskal.test(ON_Correctness ~ EOY2017, data = tbl_perf_ind_merged)

krusk_MWPerfON <- kruskal.test(MW_Correctness ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_MWPerfMW <- kruskal.test(MW_Correctness ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_MWPerfMB <- kruskal.test(MW_Correctness ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_MWPerfEOY17 <- kruskal.test(MW_Correctness ~ EOY2017, data = tbl_perf_ind_merged)

krusk_MBPerfON <- kruskal.test(MB_Correctness ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_MBPerfMW <- kruskal.test(MB_Correctness ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_MBPerfMB <- kruskal.test(MB_Correctness ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_MBPerfEOY17 <- kruskal.test(MB_Correctness ~ EOY2017, data = tbl_perf_ind_merged)

# Task performance RT
krusk_ONRTON <- kruskal.test(ON_RT ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_ONRTMW <- kruskal.test(ON_RT ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_ONRTMB <- kruskal.test(ON_RT ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_ONRTEOY17 <- kruskal.test(ON_RT ~ EOY2017, data = tbl_perf_ind_merged)

krusk_MWRTON <- kruskal.test(MW_RT ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_MWRTMW <- kruskal.test(MW_RT ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_MWRTMB <- kruskal.test(MW_RT ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_MWRTEOY17 <- kruskal.test(MW_RT ~ EOY2017, data = tbl_perf_ind_merged)

krusk_MBRTON <- kruskal.test(MB_RT ~ ON_CAT, data = tbl_perf_ind_merged)
krusk_MBRTMW <- kruskal.test(MB_RT ~ MW_CAT, data = tbl_perf_ind_merged)
krusk_MBRTMB <- kruskal.test(MB_RT ~ MB_CAT, data = tbl_perf_ind_merged)
krusk_MBRTEOY17 <- kruskal.test(MB_RT ~ EOY2017, data = tbl_perf_ind_merged)



#krusk_ASRS_ONPerc <- kruskal.test(ON_perc ~ ASRS_CAT, data = tbl_perf_ind_merged)
#krusk_ASRS_MWPerc <- kruskal.test(MW_perc ~ ASRS_CAT, data = tbl_perf_ind_merged)
#krusk_ASRS_MBPerc <- kruskal.test(MB_perc ~ ASRS_CAT, data = tbl_perf_ind_merged)
#krusk_ASRS_EOY17Perc <- kruskal.test(EOY2017 ~ ASRS_CAT, data = tbl_perf_ind_merged)
#krusk_ASRS_ONPerc
#krusk_ASRS_MWPerc
#krusk_ASRS_MBPerc
#krusk_ASRS_EOY17Perc


# See results
#DASS - Numeric
krusk_DepON #Not sig
krusk_DepMW #Not sig
krusk_DepMB #Not sig
krusk_DepEOY17 #Not sig

krusk_AnxON #Sig at p=.05
krusk_AnxMW #Not sig
krusk_AnxMB #Not sig
krusk_EOY17 #Not sig

krusk_StrON # Sig at p=.1
krusk_StrMW #Not sig
krusk_StrMB # Sig at p=.1
krusk_StrEOY17 #Not sig

#ESS - Numeric
krusk_ESSON #Not sig
krusk_ESSMW #Not sig
krusk_ESSMB #Not sig
krusk_ESSEOY17 #Not sig

#MWQ - Numeric
krusk_MWQON #Not sig
krusk_MWQMW #Not sig
krusk_MWQMB #Not sig
krusk_MWQEOY17 #Not sig

# Mini IPIP - Numeric
krusk_IPNeuON #Not sig
krusk_IPNeuMW #Not sig
krusk_IPNeuMB #Not sig
krusk_IPNeuEOY17 # Sig at p=.1

krusk_IPExtON #Not sig
krusk_IPExtMW #Not sig
krusk_IPExtMB #Not sig
krusk_IPExtEOY17 #Not sig

krusk_IPOpeON #Not sig
krusk_IPOpeMW #Not sig
krusk_IPOpeMB #Not sig
krusk_IPOpeEOY17 #Not sig

krusk_IPAgrON #Not sig
krusk_IPAgrMW #Not sig
krusk_IPAgrMB #Not sig
krusk_IPAgrEOY17 #Not sig

krusk_IPConON #Not sig
krusk_IPConMW #Not sig
krusk_IPConMB #Not sig
krusk_IPConEOY17 #Not sig

# PANAS - Numeric
krusk_PANPosON #Not sig
krusk_PANPosMW #Not sig
krusk_PANPosMB # Sig at p=.1
krusk_PANPosEOY17 #Not sig

krusk_PANNegON #Not sig
krusk_PANNegMW #Not sig
krusk_PANNegMB #Not sig
krusk_PANNegEOY17 #Not sig

#ASRS - Categorical
chisq_ASRSON #Not sig
chisq_ASRSMW #Not sig
chisq_ASRSMB #Not sig
chisq_ASRSEOY17 #Not sig

#ASRS - Numeric
krusk_ASRSON #Not sig
krusk_ASRSMW #Not sig
krusk_ASRSMB #Not sig
krusk_ASRSEOY17 #Not sig

#BRS - Numeric
krusk_BRSON #Not sig
krusk_BRSMW #Not sig
krusk_BRSMB #Not sig
krusk_BRSEOY17 #Not sig

#WDQ - Problem Solving - Numeric
krusk_WDQProbsON # Sig at p=.1
krusk_WDQProbsMW # Sig at p=.1
krusk_WDQProbsMB #Not sig
krusk_WDQProbsEOY17 #Not sig

#Guilford's AUT - Numeric
krusk_AUTON #Sig at p=.1
krusk_AUTMW #Not sig
krusk_AUTMB #Not sig
krusk_AUTEOY17 #Not sig

#EES - Employee Engagement - Numeric
krusk_EESEngON #Not sig
krusk_EESEngMW #Not sig
krusk_EESEngMB #Not sig
krusk_EESEngEOY17 #Not sig

#EOY Performance Rating - Categorical
chisq_EOY17ON #Not sig
chisq_EOY17MW #Not sig
chisq_EOY17MB #Not sig

# Task performance % correct
krusk_ONPerfON # Not sig
krusk_ONPerfMW # Not sig
krusk_ONPerfMB # Not sig
krusk_ONPerfEOY17 # Not sig

krusk_MWPerfON # Not sig
krusk_MWPerfMW # Not sig
krusk_MWPerfMB # Not sig
krusk_MWPerfEOY17 # Not sig

krusk_MBPerfON # Not sig
krusk_MBPerfMW #Sig at p=.01
krusk_MBPerfMB # Not sig
krusk_MBPerfEOY17 # Not sig

# Task performance RT
krusk_ONRTON #Sig at p=.01
krusk_ONRTMW # Not sig
krusk_ONRTMB # Not sig
krusk_ONRTEOY17 # Not sig

krusk_MWRTON #Sig at p=.01
krusk_MWRTMW #Sig at p=.01
krusk_MWRTMB # Not sig
krusk_MWRTEOY17 # Not sig

krusk_MBRTON # Not sig
krusk_MBRTMW # Not sig
krusk_MBRTMB # Not sig
krusk_MBRTEOY17 # Not sig



library(FSA)
dunnTest(DASS_Anxiety ~ ON_CAT, data = tbl_perf_ind_merged)
dunnTest(DASS_Stress ~ ON_CAT, data = tbl_perf_ind_merged)
dunnTest(DASS_Stress ~ MB_CAT, data = tbl_perf_ind_merged)
dunnTest(MiniIPIP_Neuroticism ~ EOY2017, data = tbl_perf_ind_merged)
dunnTest(PANAS_Positive_Affect ~ MB_CAT, data = tbl_perf_ind_merged)
dunnTest(WDQ_Problem_Solving ~ ON_CAT, data = tbl_perf_ind_merged)
dunnTest(WDQ_Problem_Solving ~ MW_CAT, data = tbl_perf_ind_merged)
dunnTest(AUT_Score ~ ON_CAT, data = tbl_perf_ind_merged)

# Significant models
#krusk_AnxON <- kruskal.test(DASS_Anxiety ~ ON_CAT, data = tbl_perf_ind_merged)
#krusk_StrON <- kruskal.test(DASS_Stress ~ ON_CAT, data = tbl_perf_ind_merged)
#krusk_StrMB <- kruskal.test(DASS_Stress ~ MB_CAT, data = tbl_perf_ind_merged)
#krusk_IPNeuEOY17 <- kruskal.test(MiniIPIP_Neuroticism ~ EOY2017, data = tbl_perf_ind_merged)
#krusk_PANPosMB <- kruskal.test(PANAS_Positive_Affect ~ MB_CAT, data = tbl_perf_ind_merged)
#krusk_WDQProbsON <- kruskal.test(WDQ_Problem_Solving ~ ON_CAT, data = tbl_perf_ind_merged)
#krusk_WDQProbsMW <- kruskal.test(WDQ_Problem_Solving ~ MW_CAT, data = tbl_perf_ind_merged)
#krusk_AUTON <- kruskal.test(AUT_Score ~ ON_CAT, data = tbl_perf_ind_merged)

# Descriptive stats
attach(tbl_perf_ind_merged)
describeBy(DASS_Depression, group=ON_CAT)
describeBy(DASS_Anxiety, group=ON_CAT)
describeBy(DASS_Stress, group=ON_CAT)
describeBy(ESS, group=ON_CAT)
describeBy(MWQ_mean, group=ON_CAT)
describeBy(MiniIPIP_Neuroticism, group=ON_CAT)
describeBy(MiniIPIP_Extraversion, group=ON_CAT)
describeBy(MiniIPIP_Openness, group=ON_CAT)
describeBy(MiniIPIP_Agreeableness, group=ON_CAT)
describeBy(MiniIPIP_Conscientiousness, group=ON_CAT)
describeBy(PANAS_Positive_Affect, group=ON_CAT)
describeBy(PANAS_Negative_Affect, group=ON_CAT)
describeBy(ASRS_ptA_num, group=ON_CAT)
describeBy(BRS_mean, group=ON_CAT)
describeBy(WDQ_Problem_Solving, group=ON_CAT)
describeBy(EES_Employee_Engagement, group=ON_CAT)
describeBy(AUT_Score, group=ON_CAT)
describeBy(ON_Correctness, group=ON_CAT)
describeBy(MW_Correctness, group=ON_CAT)
describeBy(MB_Correctness, group=ON_CAT)
describeBy(ON_RT, group=ON_CAT)
describeBy(MW_RT, group=ON_CAT)
describeBy(MB_RT, group=ON_CAT)
describeBy(values.ospan, group=ON_CAT)


describeBy(DASS_Depression, group=MW_CAT)
describeBy(DASS_Anxiety, group=MW_CAT)
describeBy(DASS_Stress, group=MW_CAT)
describeBy(ESS, group=MW_CAT)
describeBy(MWQ_mean, group=MW_CAT)
describeBy(MiniIPIP_Neuroticism, group=MW_CAT)
describeBy(MiniIPIP_Extraversion, group=MW_CAT)
describeBy(MiniIPIP_Openness, group=MW_CAT)
describeBy(MiniIPIP_Agreeableness, group=MW_CAT)
describeBy(MiniIPIP_Conscientiousness, group=MW_CAT)
describeBy(PANAS_Positive_Affect, group=MW_CAT)
describeBy(PANAS_Negative_Affect, group=MW_CAT)
describeBy(ASRS_ptA_num, group=MW_CAT)
describeBy(BRS_mean, group=MW_CAT)
describeBy(WDQ_Problem_Solving, group=MW_CAT)
describeBy(EES_Employee_Engagement, group=MW_CAT)
describeBy(AUT_Score, group=MW_CAT)
describeBy(ON_Correctness, group=MW_CAT)
describeBy(MW_Correctness, group=MW_CAT)
describeBy(MB_Correctness, group=MW_CAT)
describeBy(ON_RT, group=MW_CAT)
describeBy(MW_RT, group=MW_CAT)
describeBy(MB_RT, group=MW_CAT)
describeBy(values.ospan, group=MW_CAT)



describeBy(DASS_Depression, group=MB_CAT)
describeBy(DASS_Anxiety, group=MB_CAT)
describeBy(DASS_Stress, group=MB_CAT)
describeBy(ESS, group=MB_CAT)
describeBy(MWQ_mean, group=MB_CAT)
describeBy(MiniIPIP_Neuroticism, group=MB_CAT)
describeBy(MiniIPIP_Extraversion, group=MB_CAT)
describeBy(MiniIPIP_Openness, group=MB_CAT)
describeBy(MiniIPIP_Agreeableness, group=MB_CAT)
describeBy(MiniIPIP_Conscientiousness, group=MB_CAT)
describeBy(PANAS_Positive_Affect, group=MB_CAT)
describeBy(PANAS_Negative_Affect, group=MB_CAT)
describeBy(ASRS_ptA_num, group=MB_CAT)
describeBy(BRS_mean, group=MB_CAT)
describeBy(WDQ_Problem_Solving, group=MB_CAT)
describeBy(EES_Employee_Engagement, group=MB_CAT)
describeBy(AUT_Score, group=MB_CAT)
describeBy(ON_Correctness, group=MB_CAT)
describeBy(MW_Correctness, group=MB_CAT)
describeBy(MB_Correctness, group=MB_CAT)
describeBy(ON_RT, group=MB_CAT)
describeBy(MW_RT, group=MB_CAT)
describeBy(MB_RT, group=MB_CAT)
describeBy(values.ospan, group=MB_CAT)


describe(DASS_Depression)
describe(DASS_Anxiety)
describe(DASS_Stress)
describe(ESS)
describe(MWQ_mean)
describe(MiniIPIP_Neuroticism)
describe(MiniIPIP_Extraversion)
describe(MiniIPIP_Openness)
describe(MiniIPIP_Agreeableness)
describe(MiniIPIP_Conscientiousness)
describe(PANAS_Positive_Affect)
describe(PANAS_Negative_Affect)
describe(ASRS_ptA_num)
describe(BRS_mean)
describe(WDQ_Problem_Solving)
describe(EES_Employee_Engagement)
describe(AUT_Score)
describe(ON_Correctness)
describe(MW_Correctness)
describe(MB_Correctness)
describe(ON_RT)
describe(MW_RT)
describe(MB_RT)






plot(DASS_Depression ~ ON_CAT, data=tbl_perf_ind_merged)
plot(DASS_Depression ~ MW_CAT, data=tbl_perf_ind_merged)
plot(DASS_Depression ~ MB_CAT, data=tbl_perf_ind_merged)


summary(krus.aov1)






#### end ####w





Hmisc::summarize()
aggregate(WanderIM_ProbeResults$Corr, by=WanderIM_ProbeResults$nBlock, FUN="mean") 

plot(WanderIM_ProbeResults$Corr, WanderIM_ProbeResults$nBlock, main="Correctness by block", 
     xlab="Correct", ylab="Block number", pch=19)



#### CORRELATIONS AND KRUSKAL WALLIS ON SCALE DVS ####
attach(tbl_perf_ind_merged)

### CORR GO

# perform tests
corrgo_test1 <- cor.test(CorrGo,DASS_Depression, method = "spearman")
corrgo_test2 <- cor.test(CorrGo,DASS_Anxiety, method = "spearman")
corrgo_test3 <- cor.test(CorrGo,DASS_Stress, method = "spearman")
corrgo_test4 <- cor.test(CorrGo,ESS, method = "spearman")
corrgo_test5 <- cor.test(CorrGo,MWQ_mean, method = "spearman")
corrgo_test6 <- cor.test(CorrGo,MiniIPIP_Openness, method = "spearman")
corrgo_test7 <- cor.test(CorrGo,MiniIPIP_Conscientiousness, method = "spearman")
corrgo_test8 <- cor.test(CorrGo,MiniIPIP_Extraversion, method = "spearman")
corrgo_test9 <- cor.test(CorrGo,MiniIPIP_Agreeableness, method = "spearman")
corrgo_test10 <- cor.test(CorrGo,MiniIPIP_Neuroticism, method = "spearman")
corrgo_test11 <- cor.test(CorrGo,PANAS_Positive_Affect, method = "spearman")
corrgo_test12 <- cor.test(CorrGo,PANAS_Negative_Affect, method = "spearman")
corrgo_test13 <- cor.test(CorrGo,ASRS_ptA_num, method = "spearman")
corrgo_test14 <- cor.test(CorrGo,BRS_mean, method = "spearman")
corrgo_test15 <- cor.test(CorrGo,WDQ_Autonomy, method = "spearman")
corrgo_test16 <- cor.test(CorrGo,WDQ_Information_Processing, method = "spearman")
corrgo_test17 <- cor.test(CorrGo,WDQ_Job_Complexity, method = "spearman")
corrgo_test18 <- cor.test(CorrGo,WDQ_Problem_Solving, method = "spearman")
corrgo_test19 <- cor.test(CorrGo,WDQ_Skill_Variety, method = "spearman")
corrgo_test20 <- cor.test(CorrGo,WDQ_Social_Support, method = "spearman")
corrgo_test21 <- cor.test(CorrGo,WDQ_Specialization, method = "spearman")
corrgo_test22 <- cor.test(CorrGo,WDQ_Task_Significance, method = "spearman")
corrgo_test23 <- cor.test(CorrGo,WDQ_Task_Variety, method = "spearman")
corrgo_test24 <- cor.test(CorrGo,EES_Employee_Engagement, method = "spearman")
corrgo_test25 <- cor.test(CorrGo,AUT_Score.y, method = "spearman")
corrgo_test26 <- cor.test(CorrGo,values.ospan, method = "spearman")
corrgo_test27 <- cor.test(CorrGo,EOY2017Rating_num.x, method = "spearman")

corrgo_test28 <- kruskal.test(CorrGo~D1...Gender)
corrgo_test29 <- cor.test(CorrGo,D2...Age, method = "spearman")
corrgo_test30 <- kruskal.test(CorrGo~D3...mindfulness)

# export stats and apply bonferroni adjustment
corrgo_test1
p.adjust(corrgo_test1$p.value, method = "fdr", n = length(corrgo_test1$p.value))
corrgo_test2
p.adjust(corrgo_test2$p.value, method = "fdr", n = length(corrgo_test2$p.value))
corrgo_test3
p.adjust(corrgo_test3$p.value, method = "fdr", n = length(corrgo_test3$p.value))
corrgo_test4
p.adjust(corrgo_test4$p.value, method = "fdr", n = length(corrgo_test4$p.value))
corrgo_test5
p.adjust(corrgo_test5$p.value, method = "fdr", n = length(corrgo_test5$p.value))
corrgo_test6
p.adjust(corrgo_test6$p.value, method = "fdr", n = length(corrgo_test6$p.value))
corrgo_test7
p.adjust(corrgo_test7$p.value, method = "fdr", n = length(corrgo_test7$p.value))
corrgo_test8
p.adjust(corrgo_test8$p.value, method = "fdr", n = length(corrgo_test8$p.value))
corrgo_test9
p.adjust(corrgo_test9$p.value, method = "fdr", n = length(corrgo_test9$p.value))
corrgo_test10
p.adjust(corrgo_test10$p.value, method = "fdr", n = length(corrgo_test10$p.value))
corrgo_test11
p.adjust(corrgo_test11$p.value, method = "fdr", n = length(corrgo_test11$p.value))
corrgo_test12
p.adjust(corrgo_test12$p.value, method = "fdr", n = length(corrgo_test12$p.value))
corrgo_test13
p.adjust(corrgo_test13$p.value, method = "fdr", n = length(corrgo_test13$p.value))
corrgo_test14
p.adjust(corrgo_test14$p.value, method = "fdr", n = length(corrgo_test14$p.value))
corrgo_test15
p.adjust(corrgo_test15$p.value, method = "fdr", n = length(corrgo_test15$p.value))
corrgo_test16
p.adjust(corrgo_test16$p.value, method = "fdr", n = length(corrgo_test16$p.value))
corrgo_test17
p.adjust(corrgo_test17$p.value, method = "fdr", n = length(corrgo_test17$p.value))
corrgo_test18
p.adjust(corrgo_test18$p.value, method = "fdr", n = length(corrgo_test18$p.value))
corrgo_test19
p.adjust(corrgo_test19$p.value, method = "fdr", n = length(corrgo_test19$p.value))
corrgo_test20
p.adjust(corrgo_test20$p.value, method = "fdr", n = length(corrgo_test20$p.value))
corrgo_test21
p.adjust(corrgo_test21$p.value, method = "fdr", n = length(corrgo_test21$p.value))
corrgo_test22
p.adjust(corrgo_test22$p.value, method = "fdr", n = length(corrgo_test22$p.value))
corrgo_test23
p.adjust(corrgo_test23$p.value, method = "fdr", n = length(corrgo_test23$p.value))
corrgo_test24
p.adjust(corrgo_test24$p.value, method = "fdr", n = length(corrgo_test24$p.value))
corrgo_test25
p.adjust(corrgo_test25$p.value, method = "fdr", n = length(corrgo_test25$p.value))
corrgo_test26
p.adjust(corrgo_test26$p.value, method = "fdr", n = length(corrgo_test26$p.value))
corrgo_test27
p.adjust(corrgo_test27$p.value, method = "fdr", n = length(corrgo_test27$p.value))
corrgo_test28
p.adjust(corrgo_test28$p.value, method = "fdr", n = length(corrgo_test28$p.value))
corrgo_test29
p.adjust(corrgo_test29$p.value, method = "fdr", n = length(corrgo_test29$p.value))
corrgo_test30
p.adjust(corrgo_test30$p.value, method = "fdr", n = length(corrgo_test30$p.value))


### CORRNOGO

# perform tests
CorrNoGo_test1 <- cor.test(CorrNoGo,DASS_Depression, method = "spearman")
CorrNoGo_test2 <- cor.test(CorrNoGo,DASS_Anxiety, method = "spearman")
CorrNoGo_test3 <- cor.test(CorrNoGo,DASS_Stress, method = "spearman")
CorrNoGo_test4 <- cor.test(CorrNoGo,ESS, method = "spearman")
CorrNoGo_test5 <- cor.test(CorrNoGo,MWQ_mean, method = "spearman")
CorrNoGo_test6 <- cor.test(CorrNoGo,MiniIPIP_Openness, method = "spearman")
CorrNoGo_test7 <- cor.test(CorrNoGo,MiniIPIP_Conscientiousness, method = "spearman")
CorrNoGo_test8 <- cor.test(CorrNoGo,MiniIPIP_Extraversion, method = "spearman")
CorrNoGo_test9 <- cor.test(CorrNoGo,MiniIPIP_Agreeableness, method = "spearman")
CorrNoGo_test10 <- cor.test(CorrNoGo,MiniIPIP_Neuroticism, method = "spearman")
CorrNoGo_test11 <- cor.test(CorrNoGo,PANAS_Positive_Affect, method = "spearman")
CorrNoGo_test12 <- cor.test(CorrNoGo,PANAS_Negative_Affect, method = "spearman")
CorrNoGo_test13 <- cor.test(CorrNoGo,ASRS_ptA_num, method = "spearman")
CorrNoGo_test14 <- cor.test(CorrNoGo,BRS_mean, method = "spearman")
CorrNoGo_test15 <- cor.test(CorrNoGo,WDQ_Autonomy, method = "spearman")
CorrNoGo_test16 <- cor.test(CorrNoGo,WDQ_Information_Processing, method = "spearman")
CorrNoGo_test17 <- cor.test(CorrNoGo,WDQ_Job_Complexity, method = "spearman")
CorrNoGo_test18 <- cor.test(CorrNoGo,WDQ_Problem_Solving, method = "spearman")
CorrNoGo_test19 <- cor.test(CorrNoGo,WDQ_Skill_Variety, method = "spearman")
CorrNoGo_test20 <- cor.test(CorrNoGo,WDQ_Social_Support, method = "spearman")
CorrNoGo_test21 <- cor.test(CorrNoGo,WDQ_Specialization, method = "spearman")
CorrNoGo_test22 <- cor.test(CorrNoGo,WDQ_Task_Significance, method = "spearman")
CorrNoGo_test23 <- cor.test(CorrNoGo,WDQ_Task_Variety, method = "spearman")
CorrNoGo_test24 <- cor.test(CorrNoGo,EES_Employee_Engagement, method = "spearman")
CorrNoGo_test25 <- cor.test(CorrNoGo,AUT_Score.y, method = "spearman")
CorrNoGo_test26 <- cor.test(CorrNoGo,values.ospan, method = "spearman")
CorrNoGo_test27 <- cor.test(CorrNoGo,EOY2017Rating_num.x, method = "spearman")

CorrNoGo_test28 <- kruskal.test(CorrNoGo~D1...Gender)
CorrNoGo_test29 <- cor.test(CorrNoGo,D2...Age, method = "spearman")
CorrNoGo_test30 <- kruskal.test(CorrNoGo~D3...mindfulness)

# export stats and apply bonferroni adjustment
CorrNoGo_test1
p.adjust(CorrNoGo_test1$p.value, method = "fdr", n = length(CorrNoGo_test1$p.value))
CorrNoGo_test2
p.adjust(CorrNoGo_test2$p.value, method = "fdr", n = length(CorrNoGo_test2$p.value))
CorrNoGo_test3
p.adjust(CorrNoGo_test3$p.value, method = "fdr", n = length(CorrNoGo_test3$p.value))
CorrNoGo_test4
p.adjust(CorrNoGo_test4$p.value, method = "fdr", n = length(CorrNoGo_test4$p.value))
CorrNoGo_test5
p.adjust(CorrNoGo_test5$p.value, method = "fdr", n = length(CorrNoGo_test5$p.value))
CorrNoGo_test6
p.adjust(CorrNoGo_test6$p.value, method = "fdr", n = length(CorrNoGo_test6$p.value))
CorrNoGo_test7
p.adjust(CorrNoGo_test7$p.value, method = "fdr", n = length(CorrNoGo_test7$p.value))
CorrNoGo_test8
p.adjust(CorrNoGo_test8$p.value, method = "fdr", n = length(CorrNoGo_test8$p.value))
CorrNoGo_test9
p.adjust(CorrNoGo_test9$p.value, method = "fdr", n = length(CorrNoGo_test9$p.value))
CorrNoGo_test10
p.adjust(CorrNoGo_test10$p.value, method = "fdr", n = length(CorrNoGo_test10$p.value))
CorrNoGo_test11
p.adjust(CorrNoGo_test11$p.value, method = "fdr", n = length(CorrNoGo_test11$p.value))
CorrNoGo_test12
p.adjust(CorrNoGo_test12$p.value, method = "fdr", n = length(CorrNoGo_test12$p.value))
CorrNoGo_test13
p.adjust(CorrNoGo_test13$p.value, method = "fdr", n = length(CorrNoGo_test13$p.value))
CorrNoGo_test14
p.adjust(CorrNoGo_test14$p.value, method = "fdr", n = length(CorrNoGo_test14$p.value))
CorrNoGo_test15
p.adjust(CorrNoGo_test15$p.value, method = "fdr", n = length(CorrNoGo_test15$p.value))
CorrNoGo_test16
p.adjust(CorrNoGo_test16$p.value, method = "fdr", n = length(CorrNoGo_test16$p.value))
CorrNoGo_test17
p.adjust(CorrNoGo_test17$p.value, method = "fdr", n = length(CorrNoGo_test17$p.value))
CorrNoGo_test18
p.adjust(CorrNoGo_test18$p.value, method = "fdr", n = length(CorrNoGo_test18$p.value))
CorrNoGo_test19
p.adjust(CorrNoGo_test19$p.value, method = "fdr", n = length(CorrNoGo_test19$p.value))
CorrNoGo_test20
p.adjust(CorrNoGo_test20$p.value, method = "fdr", n = length(CorrNoGo_test20$p.value))
CorrNoGo_test21
p.adjust(CorrNoGo_test21$p.value, method = "fdr", n = length(CorrNoGo_test21$p.value))
CorrNoGo_test22
p.adjust(CorrNoGo_test22$p.value, method = "fdr", n = length(CorrNoGo_test22$p.value))
CorrNoGo_test23
p.adjust(CorrNoGo_test23$p.value, method = "fdr", n = length(CorrNoGo_test23$p.value))
CorrNoGo_test24
p.adjust(CorrNoGo_test24$p.value, method = "fdr", n = length(CorrNoGo_test24$p.value))
CorrNoGo_test25
p.adjust(CorrNoGo_test25$p.value, method = "fdr", n = length(CorrNoGo_test25$p.value))
CorrNoGo_test26
p.adjust(CorrNoGo_test26$p.value, method = "fdr", n = length(CorrNoGo_test26$p.value))
CorrNoGo_test27
p.adjust(CorrNoGo_test27$p.value, method = "fdr", n = length(CorrNoGo_test27$p.value))
CorrNoGo_test28
p.adjust(CorrNoGo_test28$p.value, method = "fdr", n = length(CorrNoGo_test28$p.value))
CorrNoGo_test29
p.adjust(CorrNoGo_test29$p.value, method = "fdr", n = length(CorrNoGo_test29$p.value))
CorrNoGo_test30
p.adjust(CorrNoGo_test30$p.value, method = "fdr", n = length(CorrNoGo_test30$p.value))


### RTGO


# perform tests
RTGo_test1 <- cor.test(RTGo,DASS_Depression, method = "spearman")
RTGo_test2 <- cor.test(RTGo,DASS_Anxiety, method = "spearman")
RTGo_test3 <- cor.test(RTGo,DASS_Stress, method = "spearman")
RTGo_test4 <- cor.test(RTGo,ESS, method = "spearman")
RTGo_test5 <- cor.test(RTGo,MWQ_mean, method = "spearman")
RTGo_test6 <- cor.test(RTGo,MiniIPIP_Openness, method = "spearman")
RTGo_test7 <- cor.test(RTGo,MiniIPIP_Conscientiousness, method = "spearman")
RTGo_test8 <- cor.test(RTGo,MiniIPIP_Extraversion, method = "spearman")
RTGo_test9 <- cor.test(RTGo,MiniIPIP_Agreeableness, method = "spearman")
RTGo_test10 <- cor.test(RTGo,MiniIPIP_Neuroticism, method = "spearman")
RTGo_test11 <- cor.test(RTGo,PANAS_Positive_Affect, method = "spearman")
RTGo_test12 <- cor.test(RTGo,PANAS_Negative_Affect, method = "spearman")
RTGo_test13 <- cor.test(RTGo,ASRS_ptA_num, method = "spearman")
RTGo_test14 <- cor.test(RTGo,BRS_mean, method = "spearman")
RTGo_test15 <- cor.test(RTGo,WDQ_Autonomy, method = "spearman")
RTGo_test16 <- cor.test(RTGo,WDQ_Information_Processing, method = "spearman")
RTGo_test17 <- cor.test(RTGo,WDQ_Job_Complexity, method = "spearman")
RTGo_test18 <- cor.test(RTGo,WDQ_Problem_Solving, method = "spearman")
RTGo_test19 <- cor.test(RTGo,WDQ_Skill_Variety, method = "spearman")
RTGo_test20 <- cor.test(RTGo,WDQ_Social_Support, method = "spearman")
RTGo_test21 <- cor.test(RTGo,WDQ_Specialization, method = "spearman")
RTGo_test22 <- cor.test(RTGo,WDQ_Task_Significance, method = "spearman")
RTGo_test23 <- cor.test(RTGo,WDQ_Task_Variety, method = "spearman")
RTGo_test24 <- cor.test(RTGo,EES_Employee_Engagement, method = "spearman")
RTGo_test25 <- cor.test(RTGo,AUT_Score.y, method = "spearman")
RTGo_test26 <- cor.test(RTGo,values.ospan, method = "spearman")
RTGo_test27 <- cor.test(RTGo,EOY2017Rating_num.x, method = "spearman")

RTGo_test28 <- kruskal.test(RTGo~D1...Gender)
RTGo_test29 <- cor.test(RTGo,D2...Age, method = "spearman")
RTGo_test30 <- kruskal.test(RTGo~D3...mindfulness)

# export stats and apply bonferroni adjustment
RTGo_test1
p.adjust(RTGo_test1$p.value, method = "fdr", n = length(RTGo_test1$p.value))
RTGo_test2
p.adjust(RTGo_test2$p.value, method = "fdr", n = length(RTGo_test2$p.value))
RTGo_test3
p.adjust(RTGo_test3$p.value, method = "fdr", n = length(RTGo_test3$p.value))
RTGo_test4
p.adjust(RTGo_test4$p.value, method = "fdr", n = length(RTGo_test4$p.value))
RTGo_test5
p.adjust(RTGo_test5$p.value, method = "fdr", n = length(RTGo_test5$p.value))
RTGo_test6
p.adjust(RTGo_test6$p.value, method = "fdr", n = length(RTGo_test6$p.value))
RTGo_test7
p.adjust(RTGo_test7$p.value, method = "fdr", n = length(RTGo_test7$p.value))
RTGo_test8
p.adjust(RTGo_test8$p.value, method = "fdr", n = length(RTGo_test8$p.value))
RTGo_test9
p.adjust(RTGo_test9$p.value, method = "fdr", n = length(RTGo_test9$p.value))
RTGo_test10
p.adjust(RTGo_test10$p.value, method = "fdr", n = length(RTGo_test10$p.value))
RTGo_test11
p.adjust(RTGo_test11$p.value, method = "fdr", n = length(RTGo_test11$p.value))
RTGo_test12
p.adjust(RTGo_test12$p.value, method = "fdr", n = length(RTGo_test12$p.value))
RTGo_test13
p.adjust(RTGo_test13$p.value, method = "fdr", n = length(RTGo_test13$p.value))
RTGo_test14
p.adjust(RTGo_test14$p.value, method = "fdr", n = length(RTGo_test14$p.value))
RTGo_test15
p.adjust(RTGo_test15$p.value, method = "fdr", n = length(RTGo_test15$p.value))
RTGo_test16
p.adjust(RTGo_test16$p.value, method = "fdr", n = length(RTGo_test16$p.value))
RTGo_test17
p.adjust(RTGo_test17$p.value, method = "fdr", n = length(RTGo_test17$p.value))
RTGo_test18
p.adjust(RTGo_test18$p.value, method = "fdr", n = length(RTGo_test18$p.value))
RTGo_test19
p.adjust(RTGo_test19$p.value, method = "fdr", n = length(RTGo_test19$p.value))
RTGo_test20
p.adjust(RTGo_test20$p.value, method = "fdr", n = length(RTGo_test20$p.value))
RTGo_test21
p.adjust(RTGo_test21$p.value, method = "fdr", n = length(RTGo_test21$p.value))
RTGo_test22
p.adjust(RTGo_test22$p.value, method = "fdr", n = length(RTGo_test22$p.value))
RTGo_test23
p.adjust(RTGo_test23$p.value, method = "fdr", n = length(RTGo_test23$p.value))
RTGo_test24
p.adjust(RTGo_test24$p.value, method = "fdr", n = length(RTGo_test24$p.value))
RTGo_test25
p.adjust(RTGo_test25$p.value, method = "fdr", n = length(RTGo_test25$p.value))
RTGo_test26
p.adjust(RTGo_test26$p.value, method = "fdr", n = length(RTGo_test26$p.value))
RTGo_test27
p.adjust(RTGo_test27$p.value, method = "fdr", n = length(RTGo_test27$p.value))
RTGo_test28
p.adjust(RTGo_test28$p.value, method = "fdr", n = length(RTGo_test28$p.value))
RTGo_test29
p.adjust(RTGo_test29$p.value, method = "fdr", n = length(RTGo_test29$p.value))
RTGo_test30
p.adjust(RTGo_test30$p.value, method = "fdr", n = length(RTGo_test30$p.value))



### D-PRIME 


# perform tests
dp_test1 <- cor.test(dp,DASS_Depression, method = "spearman")
dp_test2 <- cor.test(dp,DASS_Anxiety, method = "spearman")
dp_test3 <- cor.test(dp,DASS_Stress, method = "spearman")
dp_test4 <- cor.test(dp,ESS, method = "spearman")
dp_test5 <- cor.test(dp,MWQ_mean, method = "spearman")
dp_test6 <- cor.test(dp,MiniIPIP_Openness, method = "spearman")
dp_test7 <- cor.test(dp,MiniIPIP_Conscientiousness, method = "spearman")
dp_test8 <- cor.test(dp,MiniIPIP_Extraversion, method = "spearman")
dp_test9 <- cor.test(dp,MiniIPIP_Agreeableness, method = "spearman")
dp_test10 <- cor.test(dp,MiniIPIP_Neuroticism, method = "spearman")
dp_test11 <- cor.test(dp,PANAS_Positive_Affect, method = "spearman")
dp_test12 <- cor.test(dp,PANAS_Negative_Affect, method = "spearman")
dp_test13 <- cor.test(dp,ASRS_ptA_num, method = "spearman")
dp_test14 <- cor.test(dp,BRS_mean, method = "spearman")
dp_test15 <- cor.test(dp,WDQ_Autonomy, method = "spearman")
dp_test16 <- cor.test(dp,WDQ_Information_Processing, method = "spearman")
dp_test17 <- cor.test(dp,WDQ_Job_Complexity, method = "spearman")
dp_test18 <- cor.test(dp,WDQ_Problem_Solving, method = "spearman")
dp_test19 <- cor.test(dp,WDQ_Skill_Variety, method = "spearman")
dp_test20 <- cor.test(dp,WDQ_Social_Support, method = "spearman")
dp_test21 <- cor.test(dp,WDQ_Specialization, method = "spearman")
dp_test22 <- cor.test(dp,WDQ_Task_Significance, method = "spearman")
dp_test23 <- cor.test(dp,WDQ_Task_Variety, method = "spearman")
dp_test24 <- cor.test(dp,EES_Employee_Engagement, method = "spearman")
dp_test25 <- cor.test(dp,AUT_Score.y, method = "spearman")
dp_test26 <- cor.test(dp,values.ospan, method = "spearman")
dp_test27 <- cor.test(dp,EOY2017Rating_num.x, method = "spearman")

dp_test28 <- kruskal.test(dp~D1...Gender)
dp_test29 <- cor.test(dp,D2...Age, method = "spearman")
dp_test30 <- kruskal.test(dp~D3...mindfulness)

# export stats and apply bonferroni adjustment
dp_test1
p.adjust(dp_test1$p.value, method = "fdr", n = length(dp_test1$p.value))
dp_test2
p.adjust(dp_test2$p.value, method = "fdr", n = length(dp_test2$p.value))
dp_test3
p.adjust(dp_test3$p.value, method = "fdr", n = length(dp_test3$p.value))
dp_test4
p.adjust(dp_test4$p.value, method = "fdr", n = length(dp_test4$p.value))
dp_test5
p.adjust(dp_test5$p.value, method = "fdr", n = length(dp_test5$p.value))
dp_test6
p.adjust(dp_test6$p.value, method = "fdr", n = length(dp_test6$p.value))
dp_test7
p.adjust(dp_test7$p.value, method = "fdr", n = length(dp_test7$p.value))
dp_test8
p.adjust(dp_test8$p.value, method = "fdr", n = length(dp_test8$p.value))
dp_test9
p.adjust(dp_test9$p.value, method = "fdr", n = length(dp_test9$p.value))
dp_test10
p.adjust(dp_test10$p.value, method = "fdr", n = length(dp_test10$p.value))
dp_test11
p.adjust(dp_test11$p.value, method = "fdr", n = length(dp_test11$p.value))
dp_test12
p.adjust(dp_test12$p.value, method = "fdr", n = length(dp_test12$p.value))
dp_test13
p.adjust(dp_test13$p.value, method = "fdr", n = length(dp_test13$p.value))
dp_test14
p.adjust(dp_test14$p.value, method = "fdr", n = length(dp_test14$p.value))
dp_test15
p.adjust(dp_test15$p.value, method = "fdr", n = length(dp_test15$p.value))
dp_test16
p.adjust(dp_test16$p.value, method = "fdr", n = length(dp_test16$p.value))
dp_test17
p.adjust(dp_test17$p.value, method = "fdr", n = length(dp_test17$p.value))
dp_test18
p.adjust(dp_test18$p.value, method = "fdr", n = length(dp_test18$p.value))
dp_test19
p.adjust(dp_test19$p.value, method = "fdr", n = length(dp_test19$p.value))
dp_test20
p.adjust(dp_test20$p.value, method = "fdr", n = length(dp_test20$p.value))
dp_test21
p.adjust(dp_test21$p.value, method = "fdr", n = length(dp_test21$p.value))
dp_test22
p.adjust(dp_test22$p.value, method = "fdr", n = length(dp_test22$p.value))
dp_test23
p.adjust(dp_test23$p.value, method = "fdr", n = length(dp_test23$p.value))
dp_test24
p.adjust(dp_test24$p.value, method = "fdr", n = length(dp_test24$p.value))
dp_test25
p.adjust(dp_test25$p.value, method = "fdr", n = length(dp_test25$p.value))
dp_test26
p.adjust(dp_test26$p.value, method = "fdr", n = length(dp_test26$p.value))
dp_test27
p.adjust(dp_test27$p.value, method = "fdr", n = length(dp_test27$p.value))
dp_test28
p.adjust(dp_test28$p.value, method = "fdr", n = length(dp_test28$p.value))
dp_test29
p.adjust(dp_test29$p.value, method = "fdr", n = length(dp_test29$p.value))
dp_test30
p.adjust(dp_test30$p.value, method = "fdr", n = length(dp_test30$p.value))



### CRITERION


# perform tests
crit_test1 <- cor.test(crit,DASS_Depression, method = "spearman")
crit_test2 <- cor.test(crit,DASS_Anxiety, method = "spearman")
crit_test3 <- cor.test(crit,DASS_Stress, method = "spearman")
crit_test4 <- cor.test(crit,ESS, method = "spearman")
crit_test5 <- cor.test(crit,MWQ_mean, method = "spearman")
crit_test6 <- cor.test(crit,MiniIPIP_Openness, method = "spearman")
crit_test7 <- cor.test(crit,MiniIPIP_Conscientiousness, method = "spearman")
crit_test8 <- cor.test(crit,MiniIPIP_Extraversion, method = "spearman")
crit_test9 <- cor.test(crit,MiniIPIP_Agreeableness, method = "spearman")
crit_test10 <- cor.test(crit,MiniIPIP_Neuroticism, method = "spearman")
crit_test11 <- cor.test(crit,PANAS_Positive_Affect, method = "spearman")
crit_test12 <- cor.test(crit,PANAS_Negative_Affect, method = "spearman")
crit_test13 <- cor.test(crit,ASRS_ptA_num, method = "spearman")
crit_test14 <- cor.test(crit,BRS_mean, method = "spearman")
crit_test15 <- cor.test(crit,WDQ_Autonomy, method = "spearman")
crit_test16 <- cor.test(crit,WDQ_Information_Processing, method = "spearman")
crit_test17 <- cor.test(crit,WDQ_Job_Complexity, method = "spearman")
crit_test18 <- cor.test(crit,WDQ_Problem_Solving, method = "spearman")
crit_test19 <- cor.test(crit,WDQ_Skill_Variety, method = "spearman")
crit_test20 <- cor.test(crit,WDQ_Social_Support, method = "spearman")
crit_test21 <- cor.test(crit,WDQ_Specialization, method = "spearman")
crit_test22 <- cor.test(crit,WDQ_Task_Significance, method = "spearman")
crit_test23 <- cor.test(crit,WDQ_Task_Variety, method = "spearman")
crit_test24 <- cor.test(crit,EES_Employee_Engagement, method = "spearman")
crit_test25 <- cor.test(crit,AUT_Score.y, method = "spearman")
crit_test26 <- cor.test(crit,values.ospan, method = "spearman")
crit_test27 <- cor.test(crit,EOY2017Rating_num.x, method = "spearman")

crit_test28 <- kruskal.test(crit~D1...Gender)
crit_test29 <- cor.test(crit,D2...Age, method = "spearman")
crit_test30 <- kruskal.test(crit~D3...mindfulness)

# export stats and apply bonferroni adjustment
crit_test1
p.adjust(crit_test1$p.value, method = "fdr", n = length(crit_test1$p.value))
crit_test2
p.adjust(crit_test2$p.value, method = "fdr", n = length(crit_test2$p.value))
crit_test3
p.adjust(crit_test3$p.value, method = "fdr", n = length(crit_test3$p.value))
crit_test4
p.adjust(crit_test4$p.value, method = "fdr", n = length(crit_test4$p.value))
crit_test5
p.adjust(crit_test5$p.value, method = "fdr", n = length(crit_test5$p.value))
crit_test6
p.adjust(crit_test6$p.value, method = "fdr", n = length(crit_test6$p.value))
crit_test7
p.adjust(crit_test7$p.value, method = "fdr", n = length(crit_test7$p.value))
crit_test8
p.adjust(crit_test8$p.value, method = "fdr", n = length(crit_test8$p.value))
crit_test9
p.adjust(crit_test9$p.value, method = "fdr", n = length(crit_test9$p.value))
crit_test10
p.adjust(crit_test10$p.value, method = "fdr", n = length(crit_test10$p.value))
crit_test11
p.adjust(crit_test11$p.value, method = "fdr", n = length(crit_test11$p.value))
crit_test12
p.adjust(crit_test12$p.value, method = "fdr", n = length(crit_test12$p.value))
crit_test13
p.adjust(crit_test13$p.value, method = "fdr", n = length(crit_test13$p.value))
crit_test14
p.adjust(crit_test14$p.value, method = "fdr", n = length(crit_test14$p.value))
crit_test15
p.adjust(crit_test15$p.value, method = "fdr", n = length(crit_test15$p.value))
crit_test16
p.adjust(crit_test16$p.value, method = "fdr", n = length(crit_test16$p.value))
crit_test17
p.adjust(crit_test17$p.value, method = "fdr", n = length(crit_test17$p.value))
crit_test18
p.adjust(crit_test18$p.value, method = "fdr", n = length(crit_test18$p.value))
crit_test19
p.adjust(crit_test19$p.value, method = "fdr", n = length(crit_test19$p.value))
crit_test20
p.adjust(crit_test20$p.value, method = "fdr", n = length(crit_test20$p.value))
crit_test21
p.adjust(crit_test21$p.value, method = "fdr", n = length(crit_test21$p.value))
crit_test22
p.adjust(crit_test22$p.value, method = "fdr", n = length(crit_test22$p.value))
crit_test23
p.adjust(crit_test23$p.value, method = "fdr", n = length(crit_test23$p.value))
crit_test24
p.adjust(crit_test24$p.value, method = "fdr", n = length(crit_test24$p.value))
crit_test25
p.adjust(crit_test25$p.value, method = "fdr", n = length(crit_test25$p.value))
crit_test26
p.adjust(crit_test26$p.value, method = "fdr", n = length(crit_test26$p.value))
crit_test27
p.adjust(crit_test27$p.value, method = "fdr", n = length(crit_test27$p.value))
crit_test28
p.adjust(crit_test28$p.value, method = "fdr", n = length(crit_test28$p.value))
crit_test29
p.adjust(crit_test29$p.value, method = "fdr", n = length(crit_test29$p.value))
crit_test30
p.adjust(crit_test30$p.value, method = "fdr", n = length(crit_test30$p.value))



### ON PERCENTAGE


# perform tests
ON_test1 <- cor.test(ON,DASS_Depression, method = "spearman")
ON_test2 <- cor.test(ON,DASS_Anxiety, method = "spearman")
ON_test3 <- cor.test(ON,DASS_Stress, method = "spearman")
ON_test4 <- cor.test(ON,ESS, method = "spearman")
ON_test5 <- cor.test(ON,MWQ_mean, method = "spearman")
ON_test6 <- cor.test(ON,MiniIPIP_Openness, method = "spearman")
ON_test7 <- cor.test(ON,MiniIPIP_Conscientiousness, method = "spearman")
ON_test8 <- cor.test(ON,MiniIPIP_Extraversion, method = "spearman")
ON_test9 <- cor.test(ON,MiniIPIP_Agreeableness, method = "spearman")
ON_test10 <- cor.test(ON,MiniIPIP_Neuroticism, method = "spearman")
ON_test11 <- cor.test(ON,PANAS_Positive_Affect, method = "spearman")
ON_test12 <- cor.test(ON,PANAS_Negative_Affect, method = "spearman")
ON_test13 <- cor.test(ON,ASRS_ptA_num, method = "spearman")
ON_test14 <- cor.test(ON,BRS_mean, method = "spearman")
ON_test15 <- cor.test(ON,WDQ_Autonomy, method = "spearman")
ON_test16 <- cor.test(ON,WDQ_Information_Processing, method = "spearman")
ON_test17 <- cor.test(ON,WDQ_Job_Complexity, method = "spearman")
ON_test18 <- cor.test(ON,WDQ_Problem_Solving, method = "spearman")
ON_test19 <- cor.test(ON,WDQ_Skill_Variety, method = "spearman")
ON_test20 <- cor.test(ON,WDQ_Social_Support, method = "spearman")
ON_test21 <- cor.test(ON,WDQ_Specialization, method = "spearman")
ON_test22 <- cor.test(ON,WDQ_Task_Significance, method = "spearman")
ON_test23 <- cor.test(ON,WDQ_Task_Variety, method = "spearman")
ON_test24 <- cor.test(ON,EES_Employee_Engagement, method = "spearman")
ON_test25 <- cor.test(ON,AUT_Score.y, method = "spearman")
ON_test26 <- cor.test(ON,values.ospan, method = "spearman")
ON_test27 <- cor.test(ON,EOY2017Rating_num.x, method = "spearman")

ON_test28 <- kruskal.test(ON~D1...Gender)
ON_test29 <- cor.test(ON,D2...Age, method = "spearman")
ON_test30 <- kruskal.test(ON~D3...mindfulness)

# export stats and apply bonferroni adjustment
ON_test1
p.adjust(ON_test1$p.value, method = "fdr", n = length(ON_test1$p.value))
ON_test2
p.adjust(ON_test2$p.value, method = "fdr", n = length(ON_test2$p.value))
ON_test3
p.adjust(ON_test3$p.value, method = "fdr", n = length(ON_test3$p.value))
ON_test4
p.adjust(ON_test4$p.value, method = "fdr", n = length(ON_test4$p.value))
ON_test5
p.adjust(ON_test5$p.value, method = "fdr", n = length(ON_test5$p.value))
ON_test6
p.adjust(ON_test6$p.value, method = "fdr", n = length(ON_test6$p.value))
ON_test7
p.adjust(ON_test7$p.value, method = "fdr", n = length(ON_test7$p.value))
ON_test8
p.adjust(ON_test8$p.value, method = "fdr", n = length(ON_test8$p.value))
ON_test9
p.adjust(ON_test9$p.value, method = "fdr", n = length(ON_test9$p.value))
ON_test10
p.adjust(ON_test10$p.value, method = "fdr", n = length(ON_test10$p.value))
ON_test11
p.adjust(ON_test11$p.value, method = "fdr", n = length(ON_test11$p.value))
ON_test12
p.adjust(ON_test12$p.value, method = "fdr", n = length(ON_test12$p.value))
ON_test13
p.adjust(ON_test13$p.value, method = "fdr", n = length(ON_test13$p.value))
ON_test14
p.adjust(ON_test14$p.value, method = "fdr", n = length(ON_test14$p.value))
ON_test15
p.adjust(ON_test15$p.value, method = "fdr", n = length(ON_test15$p.value))
ON_test16
p.adjust(ON_test16$p.value, method = "fdr", n = length(ON_test16$p.value))
ON_test17
p.adjust(ON_test17$p.value, method = "fdr", n = length(ON_test17$p.value))
ON_test18
p.adjust(ON_test18$p.value, method = "fdr", n = length(ON_test18$p.value))
ON_test19
p.adjust(ON_test19$p.value, method = "fdr", n = length(ON_test19$p.value))
ON_test20
p.adjust(ON_test20$p.value, method = "fdr", n = length(ON_test20$p.value))
ON_test21
p.adjust(ON_test21$p.value, method = "fdr", n = length(ON_test21$p.value))
ON_test22
p.adjust(ON_test22$p.value, method = "fdr", n = length(ON_test22$p.value))
ON_test23
p.adjust(ON_test23$p.value, method = "fdr", n = length(ON_test23$p.value))
ON_test24
p.adjust(ON_test24$p.value, method = "fdr", n = length(ON_test24$p.value))
ON_test25
p.adjust(ON_test25$p.value, method = "fdr", n = length(ON_test25$p.value))
ON_test26
p.adjust(ON_test26$p.value, method = "fdr", n = length(ON_test26$p.value))
ON_test27
p.adjust(ON_test27$p.value, method = "fdr", n = length(ON_test27$p.value))
ON_test28
p.adjust(ON_test28$p.value, method = "fdr", n = length(ON_test28$p.value))
ON_test29
p.adjust(ON_test29$p.value, method = "fdr", n = length(ON_test29$p.value))
ON_test30
p.adjust(ON_test30$p.value, method = "fdr", n = length(ON_test30$p.value))



### MW PERCENTAGE


# perform tests
MW_test1 <- cor.test(MW,DASS_Depression, method = "spearman")
MW_test2 <- cor.test(MW,DASS_Anxiety, method = "spearman")
MW_test3 <- cor.test(MW,DASS_Stress, method = "spearman")
MW_test4 <- cor.test(MW,ESS, method = "spearman")
MW_test5 <- cor.test(MW,MWQ_mean, method = "spearman")
MW_test6 <- cor.test(MW,MiniIPIP_Openness, method = "spearman")
MW_test7 <- cor.test(MW,MiniIPIP_Conscientiousness, method = "spearman")
MW_test8 <- cor.test(MW,MiniIPIP_Extraversion, method = "spearman")
MW_test9 <- cor.test(MW,MiniIPIP_Agreeableness, method = "spearman")
MW_test10 <- cor.test(MW,MiniIPIP_Neuroticism, method = "spearman")
MW_test11 <- cor.test(MW,PANAS_Positive_Affect, method = "spearman")
MW_test12 <- cor.test(MW,PANAS_Negative_Affect, method = "spearman")
MW_test13 <- cor.test(MW,ASRS_ptA_num, method = "spearman")
MW_test14 <- cor.test(MW,BRS_mean, method = "spearman")
MW_test15 <- cor.test(MW,WDQ_Autonomy, method = "spearman")
MW_test16 <- cor.test(MW,WDQ_Information_Processing, method = "spearman")
MW_test17 <- cor.test(MW,WDQ_Job_Complexity, method = "spearman")
MW_test18 <- cor.test(MW,WDQ_Problem_Solving, method = "spearman")
MW_test19 <- cor.test(MW,WDQ_Skill_Variety, method = "spearman")
MW_test20 <- cor.test(MW,WDQ_Social_Support, method = "spearman")
MW_test21 <- cor.test(MW,WDQ_Specialization, method = "spearman")
MW_test22 <- cor.test(MW,WDQ_Task_Significance, method = "spearman")
MW_test23 <- cor.test(MW,WDQ_Task_Variety, method = "spearman")
MW_test24 <- cor.test(MW,EES_Employee_Engagement, method = "spearman")
MW_test25 <- cor.test(MW,AUT_Score.y, method = "spearman")
MW_test26 <- cor.test(MW,values.ospan, method = "spearman")
MW_test27 <- cor.test(MW,EOY2017Rating_num.x, method = "spearman")

MW_test28 <- kruskal.test(MW~D1...Gender)
MW_test29 <- cor.test(MW,D2...Age, method = "spearman")
MW_test30 <- kruskal.test(MW~D3...mindfulness)

# export stats and apply bonferroni adjustment
MW_test1
p.adjust(MW_test1$p.value, method = "fdr", n = length(MW_test1$p.value))
MW_test2
p.adjust(MW_test2$p.value, method = "fdr", n = length(MW_test2$p.value))
MW_test3
p.adjust(MW_test3$p.value, method = "fdr", n = length(MW_test3$p.value))
MW_test4
p.adjust(MW_test4$p.value, method = "fdr", n = length(MW_test4$p.value))
MW_test5
p.adjust(MW_test5$p.value, method = "fdr", n = length(MW_test5$p.value))
MW_test6
p.adjust(MW_test6$p.value, method = "fdr", n = length(MW_test6$p.value))
MW_test7
p.adjust(MW_test7$p.value, method = "fdr", n = length(MW_test7$p.value))
MW_test8
p.adjust(MW_test8$p.value, method = "fdr", n = length(MW_test8$p.value))
MW_test9
p.adjust(MW_test9$p.value, method = "fdr", n = length(MW_test9$p.value))
MW_test10
p.adjust(MW_test10$p.value, method = "fdr", n = length(MW_test10$p.value))
MW_test11
p.adjust(MW_test11$p.value, method = "fdr", n = length(MW_test11$p.value))
MW_test12
p.adjust(MW_test12$p.value, method = "fdr", n = length(MW_test12$p.value))
MW_test13
p.adjust(MW_test13$p.value, method = "fdr", n = length(MW_test13$p.value))
MW_test14
p.adjust(MW_test14$p.value, method = "fdr", n = length(MW_test14$p.value))
MW_test15
p.adjust(MW_test15$p.value, method = "fdr", n = length(MW_test15$p.value))
MW_test16
p.adjust(MW_test16$p.value, method = "fdr", n = length(MW_test16$p.value))
MW_test17
p.adjust(MW_test17$p.value, method = "fdr", n = length(MW_test17$p.value))
MW_test18
p.adjust(MW_test18$p.value, method = "fdr", n = length(MW_test18$p.value))
MW_test19
p.adjust(MW_test19$p.value, method = "fdr", n = length(MW_test19$p.value))
MW_test20
p.adjust(MW_test20$p.value, method = "fdr", n = length(MW_test20$p.value))
MW_test21
p.adjust(MW_test21$p.value, method = "fdr", n = length(MW_test21$p.value))
MW_test22
p.adjust(MW_test22$p.value, method = "fdr", n = length(MW_test22$p.value))
MW_test23
p.adjust(MW_test23$p.value, method = "fdr", n = length(MW_test23$p.value))
MW_test24
p.adjust(MW_test24$p.value, method = "fdr", n = length(MW_test24$p.value))
MW_test25
p.adjust(MW_test25$p.value, method = "fdr", n = length(MW_test25$p.value))
MW_test26
p.adjust(MW_test26$p.value, method = "fdr", n = length(MW_test26$p.value))
MW_test27
p.adjust(MW_test27$p.value, method = "fdr", n = length(MW_test27$p.value))
MW_test28
p.adjust(MW_test28$p.value, method = "fdr", n = length(MW_test28$p.value))
MW_test29
p.adjust(MW_test29$p.value, method = "fdr", n = length(MW_test29$p.value))
MW_test30
p.adjust(MW_test30$p.value, method = "fdr", n = length(MW_test30$p.value))




### MB PERCENTAGE


# perform tests
MB_test1 <- cor.test(MB,DASS_Depression, method = "spearman")
MB_test2 <- cor.test(MB,DASS_Anxiety, method = "spearman")
MB_test3 <- cor.test(MB,DASS_Stress, method = "spearman")
MB_test4 <- cor.test(MB,ESS, method = "spearman")
MB_test5 <- cor.test(MB,MWQ_mean, method = "spearman")
MB_test6 <- cor.test(MB,MiniIPIP_Openness, method = "spearman")
MB_test7 <- cor.test(MB,MiniIPIP_Conscientiousness, method = "spearman")
MB_test8 <- cor.test(MB,MiniIPIP_Extraversion, method = "spearman")
MB_test9 <- cor.test(MB,MiniIPIP_Agreeableness, method = "spearman")
MB_test10 <- cor.test(MB,MiniIPIP_Neuroticism, method = "spearman")
MB_test11 <- cor.test(MB,PANAS_Positive_Affect, method = "spearman")
MB_test12 <- cor.test(MB,PANAS_Negative_Affect, method = "spearman")
MB_test13 <- cor.test(MB,ASRS_ptA_num, method = "spearman")
MB_test14 <- cor.test(MB,BRS_mean, method = "spearman")
MB_test15 <- cor.test(MB,WDQ_Autonomy, method = "spearman")
MB_test16 <- cor.test(MB,WDQ_Information_Processing, method = "spearman")
MB_test17 <- cor.test(MB,WDQ_Job_Complexity, method = "spearman")
MB_test18 <- cor.test(MB,WDQ_Problem_Solving, method = "spearman")
MB_test19 <- cor.test(MB,WDQ_Skill_Variety, method = "spearman")
MB_test20 <- cor.test(MB,WDQ_Social_Support, method = "spearman")
MB_test21 <- cor.test(MB,WDQ_Specialization, method = "spearman")
MB_test22 <- cor.test(MB,WDQ_Task_Significance, method = "spearman")
MB_test23 <- cor.test(MB,WDQ_Task_Variety, method = "spearman")
MB_test24 <- cor.test(MB,EES_Employee_Engagement, method = "spearman")
MB_test25 <- cor.test(MB,AUT_Score.y, method = "spearman")
MB_test26 <- cor.test(MB,values.ospan, method = "spearman")
MB_test27 <- cor.test(MB,EOY2017Rating_num.x, method = "spearman")

MB_test28 <- kruskal.test(MB~D1...Gender)
MB_test29 <- cor.test(MB,D2...Age, method = "spearman")
MB_test30 <- kruskal.test(MB~D3...mindfulness)

# export stats and apply bonferroni adjustment
MB_test1
p.adjust(MB_test1$p.value, method = "fdr", n = length(MB_test1$p.value))
MB_test2
p.adjust(MB_test2$p.value, method = "fdr", n = length(MB_test2$p.value))
MB_test3
p.adjust(MB_test3$p.value, method = "fdr", n = length(MB_test3$p.value))
MB_test4
p.adjust(MB_test4$p.value, method = "fdr", n = length(MB_test4$p.value))
MB_test5
p.adjust(MB_test5$p.value, method = "fdr", n = length(MB_test5$p.value))
MB_test6
p.adjust(MB_test6$p.value, method = "fdr", n = length(MB_test6$p.value))
MB_test7
p.adjust(MB_test7$p.value, method = "fdr", n = length(MB_test7$p.value))
MB_test8
p.adjust(MB_test8$p.value, method = "fdr", n = length(MB_test8$p.value))
MB_test9
p.adjust(MB_test9$p.value, method = "fdr", n = length(MB_test9$p.value))
MB_test10
p.adjust(MB_test10$p.value, method = "fdr", n = length(MB_test10$p.value))
MB_test11
p.adjust(MB_test11$p.value, method = "fdr", n = length(MB_test11$p.value))
MB_test12
p.adjust(MB_test12$p.value, method = "fdr", n = length(MB_test12$p.value))
MB_test13
p.adjust(MB_test13$p.value, method = "fdr", n = length(MB_test13$p.value))
MB_test14
p.adjust(MB_test14$p.value, method = "fdr", n = length(MB_test14$p.value))
MB_test15
p.adjust(MB_test15$p.value, method = "fdr", n = length(MB_test15$p.value))
MB_test16
p.adjust(MB_test16$p.value, method = "fdr", n = length(MB_test16$p.value))
MB_test17
p.adjust(MB_test17$p.value, method = "fdr", n = length(MB_test17$p.value))
MB_test18
p.adjust(MB_test18$p.value, method = "fdr", n = length(MB_test18$p.value))
MB_test19
p.adjust(MB_test19$p.value, method = "fdr", n = length(MB_test19$p.value))
MB_test20
p.adjust(MB_test20$p.value, method = "fdr", n = length(MB_test20$p.value))
MB_test21
p.adjust(MB_test21$p.value, method = "fdr", n = length(MB_test21$p.value))
MB_test22
p.adjust(MB_test22$p.value, method = "fdr", n = length(MB_test22$p.value))
MB_test23
p.adjust(MB_test23$p.value, method = "fdr", n = length(MB_test23$p.value))
MB_test24
p.adjust(MB_test24$p.value, method = "fdr", n = length(MB_test24$p.value))
MB_test25
p.adjust(MB_test25$p.value, method = "fdr", n = length(MB_test25$p.value))
MB_test26
p.adjust(MB_test26$p.value, method = "fdr", n = length(MB_test26$p.value))
MB_test27
p.adjust(MB_test27$p.value, method = "fdr", n = length(MB_test27$p.value))
MB_test28
p.adjust(MB_test28$p.value, method = "fdr", n = length(MB_test28$p.value))
MB_test29
p.adjust(MB_test29$p.value, method = "fdr", n = length(MB_test29$p.value))
MB_test30
p.adjust(MB_test30$p.value, method = "fdr", n = length(MB_test30$p.value))


plot(MB,PANAS_Positive_Affect)
plot(MB,WDQ_Information_Processing)
plot(MB,WDQ_Job_Complexity)

barplot


plot(MB_test11)

#### THOMAS RESULTS FEEDBACK ####

## Frequency of attention states
# Aggregate to Frequency of MB/MW/ON
## Retain task type variable
tbl_freq_probe1 <- subset(WanderIM_ProbeResults, DistProbe >-2)
attach(tbl_freq_probe1)
tbl_freq_probe1 <- aggregate(nTrial~SubID+State+Task, FUN = "length")
#Pivot so attention each attention state is own variable
library(reshape2)
tbl_freq_probe1 <- dcast(tbl_freq_probe1, SubID+Task~State)
#Rename variables
library(plyr)
tbl_freq_probe1 <- plyr::rename(tbl_freq_probe1, c("1"="ON_Freq", "2"="MW_Freq", "3"="MB_Freq", "4"="DR"))
# Drop don't remember probes
library(dplyr)
tbl_freq_probe1 <- select(tbl_freq_probe1, -DR)
detach(tbl_freq_probe1)

tbl_freq_probe1$ON_Freq[is.na(tbl_freq_probe1$ON_Freq)] <- 0
tbl_freq_probe1$MW_Freq[is.na(tbl_freq_probe1$MW_Freq)] <- 0
tbl_freq_probe1$MB_Freq[is.na(tbl_freq_probe1$MB_Freq)] <- 0

# t-tests to show no effect of task type on proportion of mind-state frequency
as.factor(tbl_freq_probe1$Task)
is.numeric(tbl_freq_probe1$ON_Freq)
on.perc.ttest <- t.test(tbl_freq_probe1$ON_Freq~tbl_freq_probe1$Task, paired = TRUE)
on.perc.ttest
aggregate(tbl_freq_probe1$ON_Freq~tbl_freq_probe1$Task, FUN=sd)
aggregate(tbl_freq_probe1$ON_Freq~tbl_freq_probe1$Task, FUN=mean)

mw.perc.ttest <- t.test(tbl_freq_probe1$MW_Freq~tbl_freq_probe1$Task, paired = TRUE)
mw.perc.ttest
aggregate(tbl_freq_probe1$MW_Freq~tbl_freq_probe1$Task, FUN=sd)
aggregate(tbl_freq_probe1$MW_Freq~tbl_freq_probe1$Task, FUN=mean)

mb.perc.ttest <- t.test(tbl_freq_probe1$MB_Freq~tbl_freq_probe1$Task, paired = TRUE, na.action = na.omit)
mb.perc.ttest
aggregate(tbl_freq_probe1$MB_Freq~tbl_freq_probe1$Task, FUN=sd)
aggregate(tbl_freq_probe1$MB_Freq~tbl_freq_probe1$Task, FUN=mean)


# Aggregate to origin of MW
## Retain task type variable
tbl_freq_probe2 <- subset(WanderIM_ProbeResults, DistProbe >-2)
attach(tbl_freq_probe2)
tbl_freq_probe2 <- aggregate(nTrial~SubID+Orig+Task, FUN = "length")
#Pivot so attention each attention state is own variable
library(reshape2)
tbl_freq_probe2 <- dcast(tbl_freq_probe2, SubID+Task~Orig)
#Rename variables
library(plyr)
tbl_freq_probe2 <- plyr::rename(tbl_freq_probe2, c("1"="Room", "2"="Pers", "3"="TaskOrig"))
detach(tbl_freq_probe2)

# t-tests to show no effect of task type on proportion of mind-state frequency
as.factor(tbl_freq_probe2$Task)
is.numeric(tbl_freq_probe2$Room)

tbl_freq_probe2$Room[is.na(tbl_freq_probe2$Room)] <- 0
tbl_freq_probe2$Pers[is.na(tbl_freq_probe2$Pers)] <- 0
tbl_freq_probe2$TaskOrig[is.na(tbl_freq_probe2$TaskOrig)] <- 0

room.ttest <- t.test(tbl_freq_probe2$Room~tbl_freq_probe2$Task, paired = TRUE)
room.ttest
aggregate(tbl_freq_probe2$Room~tbl_freq_probe2$Task, FUN=sd)

pers.ttest <- t.test(tbl_freq_probe2$Pers~tbl_freq_probe2$Task, paired = TRUE)
pers.ttest
aggregate(tbl_freq_probe2$Pers~tbl_freq_probe2$Task, FUN=sd)

task.ttest <- t.test(tbl_freq_probe2$TaskOrig~tbl_freq_probe2$Task, paired = TRUE)
task.ttest
aggregate(tbl_freq_probe2$TaskOrig~tbl_freq_probe2$Task, FUN=sd)

