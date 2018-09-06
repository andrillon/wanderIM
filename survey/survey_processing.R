#### IMPORT DATA ####
# -- NOTE
# Data is extracted from qualtrics in legacy format as a csv file with numeric values

# -- Import data file exactly as exported from qualtrics and remove first line of metadata
#Import csv file as exported from qualtrics
library(readr)
Attention_Scales <- read_csv("./data/input/Attention_Scales.csv")
#Remove first line of metadata
Attention_Scales <- Attention_Scales[-1,]
#Write data to new text file for import
write.table(Attention_Scales, file = "./data/output/Attention_Scales.txt",row.names=FALSE, na="",col.names=TRUE, sep="\t")
#Remove cleaned csv file
remove(Attention_Scales)

# -- Import modified text file just created
rawdf <- read.delim("./data/output/Attention_Scales.txt")
View(rawdf)

# Import Manually scored Guilford's AUT data
library(readr)
Guilford_AUT <- read_csv("data/input/AUT Input R.csv")

# Import OSPAN
OSPAN.dat <- read_csv("data/input/automatedospan_summary_18_07_09.csv")
attach(OSPAN.dat)
myvars <- c("script.subjectid", "values.ospan")
OSPAN.dat <- OSPAN.dat[myvars]
detach(OSPAN.dat)

#Import NAB data
NAB.df <- read.delim("data/input/Shortened Participant Output.txt")
View(NAB.df)
as.factor(NAB.df$Participant_No)

attach(NAB.df)
myvars <- c("Participant_No", "EOY2017Rating", "Org.Tenure", "Position.Tenure", "EmployeeLevel")
NAB.df <- NAB.df[myvars]
detach(NAB.df)

#Recode performance variable from nab.df
library(stringr)
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 

NAB.df$EOY2017Rating_num <- numextract(NAB.df$EOY2017Rating)

# -- Libraries
library("car")
library("dplyr")
library("plyr")

#### DATA CLEANING ####

# --- Data cleaning
#Demographics
rawdf$Participant_No <- rawdf$V4
rawdf$D1...Gender <- mapvalues(rawdf$D1...Gender, from = c(1,2,3,4), to = c("Male","Female","Other","Prefer not to answer"))
rawdf$D3...mindfulness <- mapvalues(rawdf$D3...mindfulness, from = c(1,2,3,4,5), to = c("Every day","A few times per week","A few times per month","A few times per year or less", "I have never practiced mindfulness"))

#DASS recodes
rawdf$DASS_1 <- mapvalues(rawdf$DASS_1, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_2 <- mapvalues(rawdf$DASS_2, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_3 <- mapvalues(rawdf$DASS_3, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_4 <- mapvalues(rawdf$DASS_4, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_5 <- mapvalues(rawdf$DASS_5, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_6 <- mapvalues(rawdf$DASS_6, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_7 <- mapvalues(rawdf$DASS_7, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_8 <- mapvalues(rawdf$DASS_8, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_9 <- mapvalues(rawdf$DASS_9, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_10 <- mapvalues(rawdf$DASS_10, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_11 <- mapvalues(rawdf$DASS_11, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_12 <- mapvalues(rawdf$DASS_12, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_13 <- mapvalues(rawdf$DASS_13, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_14 <- mapvalues(rawdf$DASS_14, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_15 <- mapvalues(rawdf$DASS_15, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_16 <- mapvalues(rawdf$DASS_16, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_17 <- mapvalues(rawdf$DASS_17, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_18 <- mapvalues(rawdf$DASS_18, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_19 <- mapvalues(rawdf$DASS_19, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_20 <- mapvalues(rawdf$DASS_20, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$DASS_21 <- mapvalues(rawdf$DASS_21, from = c(1,2,3,4), to = c(0,1,2,3))

#ESS value recodes
rawdf$ESS_1 <- mapvalues(rawdf$ESS_1, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$ESS_2 <- mapvalues(rawdf$ESS_2, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$ESS_3 <- mapvalues(rawdf$ESS_3, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$ESS_4 <- mapvalues(rawdf$ESS_4, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$ESS_5 <- mapvalues(rawdf$ESS_5, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$ESS_6 <- mapvalues(rawdf$ESS_6, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$ESS_7 <- mapvalues(rawdf$ESS_7, from = c(1,2,3,4), to = c(0,1,2,3))
rawdf$ESS_8 <- mapvalues(rawdf$ESS_8, from = c(1,2,3,4), to = c(0,1,2,3))

# MWQ questions do not require recoding

# Mini IPIP - items  6, 7, 8, 9, 10, 15, 16, 17, 18, 19, and 20 need to first be reverse scored before scoring the five domains
rawdf$Mini.IPIP_6 <- mapvalues(rawdf$Mini.IPIP_6, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$Mini.IPIP_7 <- mapvalues(rawdf$Mini.IPIP_7, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$Mini.IPIP_8 <- mapvalues(rawdf$Mini.IPIP_8, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$Mini.IPIP_9 <- mapvalues(rawdf$Mini.IPIP_9, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$Mini.IPIP_10 <- mapvalues(rawdf$Mini.IPIP_10, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$Mini.IPIP_15<- mapvalues(rawdf$Mini.IPIP_15, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$Mini.IPIP_16 <- mapvalues(rawdf$Mini.IPIP_16, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$Mini.IPIP_17 <- mapvalues(rawdf$Mini.IPIP_17, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$Mini.IPIP_18 <- mapvalues(rawdf$Mini.IPIP_18, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$Mini.IPIP_19 <- mapvalues(rawdf$Mini.IPIP_19, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$Mini.IPIP_20 <- mapvalues(rawdf$Mini.IPIP_20, from = c(1,2,3,4,5), to = c(5,4,3,2,1))

# PANAS questions do not require recoding

# ASRS questions - recode Part A into binary for 'dark shaded boxes' - analysis is done in a binary manner: >4 means the patient has symptoms highly consistent with ADHD in adults and further investigation is warranted.
rawdf$ASRS_A1 <- mapvalues(rawdf$ASRS_1, from = c(1,2,3,4,5), to = c(0,0,1,1,1))
rawdf$ASRS_A2 <- mapvalues(rawdf$ASRS_2, from = c(1,2,3,4,5), to = c(0,0,1,1,1))
rawdf$ASRS_A3 <- mapvalues(rawdf$ASRS_3, from = c(1,2,3,4,5), to = c(0,0,1,1,1))
rawdf$ASRS_A4 <- mapvalues(rawdf$ASRS_4, from = c(1,2,3,4,5), to = c(0,0,0,1,1))
rawdf$ASRS_A5 <- mapvalues(rawdf$ASRS_5, from = c(1,2,3,4,5), to = c(0,0,0,1,1))
rawdf$ASRS_A6 <- mapvalues(rawdf$ASRS_6, from = c(1,2,3,4,5), to = c(0,0,0,1,1))

# BRS - The BRS is scored by reverse coding items 2, 4, and 6 and finding the mean of the six items
rawdf$BRS_2 <- mapvalues(rawdf$BRS_2, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$BRS_4 <- mapvalues(rawdf$BRS_4, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$BRS_6 <- mapvalues(rawdf$BRS_6, from = c(1,2,3,4,5), to = c(5,4,3,2,1))

# WDQ - Recode: Scoring export from qualtrics has added +90 to the task items and +5 to the social items
rawdf$WDQ._.Task_1 <- mapvalues(rawdf$WDQ._.Task_1, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_2 <- mapvalues(rawdf$WDQ._.Task_2, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_3 <- mapvalues(rawdf$WDQ._.Task_3, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_4 <- mapvalues(rawdf$WDQ._.Task_4, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_5 <- mapvalues(rawdf$WDQ._.Task_5, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_6 <- mapvalues(rawdf$WDQ._.Task_6, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_7 <- mapvalues(rawdf$WDQ._.Task_7, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_8 <- mapvalues(rawdf$WDQ._.Task_8, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
#rename WDQ._.Task_41
rawdf$WDQ._.Task_9 <- mapvalues(rawdf$WDQ._.Task_41, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_10 <- mapvalues(rawdf$WDQ._.Task_10, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_11 <- mapvalues(rawdf$WDQ._.Task_11, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_12 <- mapvalues(rawdf$WDQ._.Task_12, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_13 <- mapvalues(rawdf$WDQ._.Task_13, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_14 <- mapvalues(rawdf$WDQ._.Task_14, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_15 <- mapvalues(rawdf$WDQ._.Task_15, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_16 <- mapvalues(rawdf$WDQ._.Task_16, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_17 <- mapvalues(rawdf$WDQ._.Task_17, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_18 <- mapvalues(rawdf$WDQ._.Task_18, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_19 <- mapvalues(rawdf$WDQ._.Task_19, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_20 <- mapvalues(rawdf$WDQ._.Task_20, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_21 <- mapvalues(rawdf$WDQ._.Task_21, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_22 <- mapvalues(rawdf$WDQ._.Task_22, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
#rename WDQ._.Task_42
rawdf$WDQ._.Task_23 <- rawdf$WDQ._.Task_42
rawdf$WDQ._.Task_23 <- mapvalues(rawdf$WDQ._.Task_23, from = c(91,92,93,94,95), to = c(1,2,3,4,5))
rawdf$WDQ._.Task_24 <- mapvalues(rawdf$WDQ._.Task_24, from = c(91,92,93,94,95), to = c(1,2,3,4,5))

#note to remove knowledge23
rawdf$WDQ._.Knowledge_41 <- rawdf$WDQ._.Knowledge_23
#note to remove knowledge9
rawdf$WDQ._.Knowledge_42 <- rawdf$WDQ._.Knowledge_9
#note to remove knowledge44
rawdf$WDQ._.Knowledge_43 <- rawdf$WDQ._.Knowledge_44

rawdf$WDQ...Social_1 <- mapvalues(rawdf$WDQ...Social_1, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_2 <- mapvalues(rawdf$WDQ...Social_2, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_3 <- mapvalues(rawdf$WDQ...Social_3, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_4 <- mapvalues(rawdf$WDQ...Social_4, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_5 <- mapvalues(rawdf$WDQ...Social_5, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_6 <- mapvalues(rawdf$WDQ...Social_6, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_7 <- mapvalues(rawdf$WDQ...Social_7, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_8 <- mapvalues(rawdf$WDQ...Social_8, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_9 <- mapvalues(rawdf$WDQ...Social_9, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_10 <- mapvalues(rawdf$WDQ...Social_10, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_11 <- mapvalues(rawdf$WDQ...Social_11, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_12 <- mapvalues(rawdf$WDQ...Social_12, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_13 <- mapvalues(rawdf$WDQ...Social_13, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_14 <- mapvalues(rawdf$WDQ...Social_14, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_15 <- mapvalues(rawdf$WDQ...Social_15, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
# remove social21
rawdf$WDQ...Social_16 <- mapvalues(rawdf$WDQ...Social_21, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
rawdf$WDQ...Social_17 <- mapvalues(rawdf$WDQ...Social_17, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
#remove social22
rawdf$WDQ...Social_18 <- mapvalues(rawdf$WDQ...Social_22, from = c(6,7,8,9,10), to = c(1,2,3,4,5))
#remove social23
rawdf$WDQ...Social_19 <- mapvalues(rawdf$WDQ...Social_23, from = c(6,7,8,9,10), to = c(1,2,3,4,5))

#Remove WDQ variables incorrectly named by qualtrics
rawdf <- subset(rawdf, select = -c(WDQ._.Task_41, WDQ._.Task_42, WDQ._.Knowledge_23, WDQ._.Knowledge_9, WDQ._.Knowledge_44, WDQ...Social_21, WDQ...Social_22, WDQ...Social_23))

# WDQ - reverse scoring questions
rawdf$WDQ._.Knowledge_25 <- mapvalues(rawdf$WDQ._.Knowledge_25, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$WDQ._.Knowledge_26 <- mapvalues(rawdf$WDQ._.Knowledge_26, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$WDQ._.Knowledge_27 <- mapvalues(rawdf$WDQ._.Knowledge_27, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
rawdf$WDQ._.Knowledge_28 <- mapvalues(rawdf$WDQ._.Knowledge_28, from = c(1,2,3,4,5), to = c(5,4,3,2,1))

#### CALCULATIONS ####

# -- Create calculated variables for the processing of scales
#DASS numeric variables
rawdf$DASS_Depression <- rawdf$DASS_3 + rawdf$DASS_5 + rawdf$DASS_10 + rawdf$DASS_13 + rawdf$DASS_16 + rawdf$DASS_17 + rawdf$DASS_21
rawdf$DASS_Anxiety <- rawdf$DASS_2 + rawdf$DASS_4 + rawdf$DASS_7 + rawdf$DASS_9 + rawdf$DASS_15 + rawdf$DASS_19 + rawdf$DASS_20
rawdf$DASS_Stress <- rawdf$DASS_1 + rawdf$DASS_6 + rawdf$DASS_8 + rawdf$DASS_11 + rawdf$DASS_12 + rawdf$DASS_14 + rawdf$DASS_18

#DASS categorical variables - calculate and recode to categorical conventional severity labels
DASS_Depression_CAT_calc <- rawdf$DASS_Depression*2
attach(rawdf)
rawdf$DASS_Depression_CAT[DASS_Depression_CAT_calc <=9] <- 'Normal'
rawdf$DASS_Depression_CAT[DASS_Depression_CAT_calc >=10 & DASS_Depression_CAT_calc <=13] <- 'Mild'
rawdf$DASS_Depression_CAT[DASS_Depression_CAT_calc >=14 & DASS_Depression_CAT_calc <=20] <- 'Moderate'
rawdf$DASS_Depression_CAT[DASS_Depression_CAT_calc >=21 & DASS_Depression_CAT_calc <=27] <- 'Severe'
rawdf$DASS_Depression_CAT[DASS_Depression_CAT_calc >=28] <- 'Extremely Severe'
detach(rawdf)
remove(DASS_Depression_CAT_calc)

DASS_Anxiety_CAT_calc <- rawdf$DASS_Anxiety*2
attach(rawdf)
rawdf$DASS_Anxiety_CAT[DASS_Anxiety_CAT_calc <=7] <- 'Normal'
rawdf$DASS_Anxiety_CAT[DASS_Anxiety_CAT_calc >=8 & DASS_Anxiety_CAT_calc <=9] <- 'Mild'
rawdf$DASS_Anxiety_CAT[DASS_Anxiety_CAT_calc >=10 & DASS_Anxiety_CAT_calc <=14] <- 'Moderate'
rawdf$DASS_Anxiety_CAT[DASS_Anxiety_CAT_calc >=15 & DASS_Anxiety_CAT_calc <=19] <- 'Severe'
rawdf$DASS_Anxiety_CAT[DASS_Anxiety_CAT_calc >=20] <- 'Extremely Severe'
detach(rawdf)
remove(DASS_Anxiety_CAT_calc)

DASS_Stress_CAT_calc <- rawdf$DASS_Stress*2
attach(rawdf)
rawdf$DASS_Stress_CAT[DASS_Stress_CAT_calc <=14] <- 'Normal'
rawdf$DASS_Stress_CAT[DASS_Stress_CAT_calc >=15 & DASS_Stress_CAT_calc <=18] <- 'Mild'
rawdf$DASS_Stress_CAT[DASS_Stress_CAT_calc >=19 & DASS_Stress_CAT_calc <=25] <- 'Moderate'
rawdf$DASS_Stress_CAT[DASS_Stress_CAT_calc >=26 & DASS_Stress_CAT_calc <=33] <- 'Severe'
rawdf$DASS_Stress_CAT[DASS_Stress_CAT_calc >=34] <- 'Extremely Severe'
detach(rawdf)
remove(DASS_Stress_CAT_calc)

#ESS numeric variable
rawdf$ESS <- rawdf$ESS_1 + rawdf$ESS_2 + rawdf$ESS_3 + rawdf$ESS_4 + rawdf$ESS_5 + rawdf$ESS_6 + rawdf$ESS_7 + rawdf$ESS_8

#ESS Categorical Interpretation Calculation
attach(rawdf)
rawdf$ESS_CAT[ESS <=5] <- 'Lower Normal Daytime Sleepiness'
rawdf$ESS_CAT[ESS >=6 & rawdf$ESS <=10] <- 'Higher Normal Daytime Sleepiness'
rawdf$ESS_CAT[ESS >=11 & rawdf$ESS <=12] <- 'Mild Excessive Daytime Sleepiness'
rawdf$ESS_CAT[ESS >=13 & rawdf$ESS <=15] <- 'Moderate Excessive Daytime Sleepiness'
rawdf$ESS_CAT[ESS >=16] <- 'Severe Excessive Daytime Sleepiness'
detach(rawdf)

#MWQ numeric variable - documentation on the MWQ only refers to "high scores" - was unable to find clear scoring instructions
rawdf$MWQ_mean <- (rawdf$MWQ_1 + rawdf$MWQ_2 + rawdf$MWQ_3 + rawdf$MWQ_4 + rawdf$MWQ_5)/5

# Mini IPIP - reverse scored items have already been computed above
rawdf$MiniIPIP_Neuroticism <- (rawdf$Mini.IPIP_4 + rawdf$Mini.IPIP_9 + rawdf$Mini.IPIP_14 + rawdf$Mini.IPIP_19)/4
rawdf$MiniIPIP_Extraversion <- (rawdf$Mini.IPIP_1 + rawdf$Mini.IPIP_6 + rawdf$Mini.IPIP_11 + rawdf$Mini.IPIP_16)/4
rawdf$MiniIPIP_Openness <- (rawdf$Mini.IPIP_5 + rawdf$Mini.IPIP_10 + rawdf$Mini.IPIP_15 + rawdf$Mini.IPIP_20)/4
rawdf$MiniIPIP_Agreeableness <- (rawdf$Mini.IPIP_2 + rawdf$Mini.IPIP_7 + rawdf$Mini.IPIP_12 + rawdf$Mini.IPIP_17)/4
rawdf$MiniIPIP_Conscientiousness <- (rawdf$Mini.IPIP_3 + rawdf$Mini.IPIP_8 + rawdf$Mini.IPIP_13 + rawdf$Mini.IPIP_18)/4

# PANAS Scoring
#Positive Affect - Scores can range from 10 – 50. Higher scores represent higher levels of positive affect. Mean Scores: Weekly   33.3 (SD   7.2)
rawdf$PANAS_Positive_Affect <- rawdf$PANAS_1 + rawdf$PANAS_3 + rawdf$PANAS_5 + rawdf$PANAS_9 + rawdf$PANAS_10 + rawdf$PANAS_12 + 
  rawdf$PANAS_14 + rawdf$PANAS_16 + rawdf$PANAS_17 + rawdf$PANAS_19
#Negative Affect - Scores can range from 10 – 50. Lower scores represent lower levels of negative affect. Mean Score: Weekly   17.4 (SD   6.2)
rawdf$PANAS_Negative_Affect <- rawdf$PANAS_2 + rawdf$PANAS_4 + rawdf$PANAS_6 + rawdf$PANAS_7 + rawdf$PANAS_8 + rawdf$PANAS_11 + 
  rawdf$PANAS_13 + rawdf$PANAS_15 + rawdf$PANAS_18 + rawdf$PANAS_20

#ASRS 
#binary Scoring calc
rawdf$ASRS_ptA_num <- rawdf$ASRS_A1 + rawdf$ASRS_A2 + rawdf$ASRS_A3 + rawdf$ASRS_A4 + rawdf$ASRS_A5 + rawdf$ASRS_A6
#Recode categorical variable for section A shaded boxes >4 :>4 means the patient has symptoms highly consistent with ADHD in adults and further investigation is warranted
attach(rawdf)
rawdf$ASRS_CAT[ASRS_ptA_num >=4] <- 'Symptoms highly consistent with ADHD'
rawdf$ASRS_CAT[ASRS_ptA_num <4] <- 'Symptoms not consistent with ADHD'
detach(rawdf)
#The remaining ASRS questions (questions 7-18:section B) provide additional cues and can serve as further probes into the patient’s symptoms

# BRS (reverse scoring has already been performed above) The BRS is scored by reverse coding items 2, 4, and 6 and finding the mean of the six items
rawdf$BRS_mean <- (rawdf$BRS_1 + rawdf$BRS_2 + rawdf$BRS_3 + rawdf$BRS_4 + rawdf$BRS_5 + rawdf$BRS_6)/6

#WDQ - Mean scores calculated from relevant dimensions
rawdf$WDQ_Autonomy <- (rawdf$WDQ._.Task_1 + rawdf$WDQ._.Task_2 + rawdf$WDQ._.Task_3 + rawdf$WDQ._.Task_4 + rawdf$WDQ._.Task_5 +
                              rawdf$WDQ._.Task_6 + rawdf$WDQ._.Task_7 + rawdf$WDQ._.Task_8 + rawdf$WDQ._.Task_9)/9
rawdf$WDQ_Task_Variety <- (rawdf$WDQ._.Task_10 + rawdf$WDQ._.Task_11 + rawdf$WDQ._.Task_12 + rawdf$WDQ._.Task_13)/4
rawdf$WDQ_Task_Significance <-  (rawdf$WDQ._.Task_14 + rawdf$WDQ._.Task_15 + rawdf$WDQ._.Task_16 + rawdf$WDQ._.Task_17)/4
rawdf$WDQ_Job_Complexity <- (rawdf$WDQ._.Knowledge_25 + rawdf$WDQ._.Knowledge_26 + rawdf$WDQ._.Knowledge_27 + rawdf$WDQ._.Knowledge_28)/4
rawdf$WDQ_Information_Processing <- (rawdf$WDQ._.Knowledge_29 + rawdf$WDQ._.Knowledge_30 + rawdf$WDQ._.Knowledge_31 + rawdf$WDQ._.Knowledge_32)/4
rawdf$WDQ_Problem_Solving <- (rawdf$WDQ._.Knowledge_33 + rawdf$WDQ._.Knowledge_34 + rawdf$WDQ._.Knowledge_35 + rawdf$WDQ._.Knowledge_36)/4
rawdf$WDQ_Skill_Variety <- (rawdf$WDQ._.Knowledge_37 + rawdf$WDQ._.Knowledge_38 + rawdf$WDQ._.Knowledge_39 + rawdf$WDQ._.Knowledge_40)/4
# Specialization is calculated from only 3 of the 4 questions included in the original WDQ due to missing data
rawdf$WDQ_Specialization <- (rawdf$WDQ._.Knowledge_40 + rawdf$WDQ._.Knowledge_42 + rawdf$WDQ._.Knowledge_43)/3
rawdf$WDQ_Social_Support <- (rawdf$WDQ...Social_1 + rawdf$WDQ...Social_2 + rawdf$WDQ...Social_3 + rawdf$WDQ...Social_4)/4

# EES - scored by summing the three subdimensions. The authors (Schuck) also note that there is a higher order construct of employee engagement measured by the three subconstructs. Therefore, employee engagement will also be computed.
# Scores on the three subscales range from 4-20 while scores on the higher order employee engagement construct range from 16-60
rawdf$EES_Emotional_Engagement <- rawdf$EES_1 + rawdf$EES_2 + rawdf$EES_3 + rawdf$EES_4
rawdf$EES_Behavioural_Engagement <-  rawdf$EES_5 + rawdf$EES_6 + rawdf$EES_7 + rawdf$EES_8
rawdf$EES_Cognitive_Engagement <-  rawdf$EES_9 + rawdf$EES_10 + rawdf$EES_11 + rawdf$EES_12
rawdf$EES_Employee_Engagement <- rawdf$EES_Emotional_Engagement + rawdf$EES_Behavioural_Engagement + rawdf$EES_Cognitive_Engagement


# -- Create a shortened output containing only computed and other relevant variables
Computed_Scales <- subset(rawdf, select = c(
  Participant_No, D1...Gender, D2...Age, D3...mindfulness,
  DASS_Depression, DASS_Anxiety, DASS_Stress, DASS_Depression_CAT, DASS_Anxiety_CAT, DASS_Stress_CAT,
  ESS, ESS_CAT,
  MWQ_mean,
  MiniIPIP_Neuroticism, MiniIPIP_Extraversion, MiniIPIP_Openness, MiniIPIP_Agreeableness, MiniIPIP_Conscientiousness,
  PANAS_Positive_Affect, PANAS_Negative_Affect,
  ASRS_ptA_num, ASRS_CAT,
  BRS_mean,
  WDQ_Autonomy, WDQ_Task_Variety, WDQ_Task_Significance, WDQ_Job_Complexity, WDQ_Information_Processing, WDQ_Problem_Solving, WDQ_Skill_Variety, WDQ_Specialization, WDQ_Social_Support,
  EES_Emotional_Engagement, EES_Behavioural_Engagement, EES_Cognitive_Engagement, EES_Employee_Engagement
))


#### JOIN DATA FILES ####

Computed_Scales <- merge(Computed_Scales, Guilford_AUT, by.x="Participant_No", by.y="ID")
Computed_Scales <- merge(Computed_Scales, OSPAN.dat, by.x="Participant_No", by.y="script.subjectid")
Computed_Scales <- merge(Computed_Scales,NAB.df, by.x="Participant_No", by.y="Participant_No")


#### EXPORT ####

View(Computed_Scales)

write.table(rawdf, file = "./data/output/MWI_Cleaned_Raw_Srv_Data.csv",row.names=FALSE, na="",col.names=TRUE, sep=",")
write.table(Computed_Scales, file = "./data/output/MWI_Computed_Scales.csv",row.names=FALSE, na="",col.names=TRUE, sep=",")


