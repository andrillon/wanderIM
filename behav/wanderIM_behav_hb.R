t<-read.table("~/Data/WanderIM/behav/WanderIM_ProbeResults_HB.txt", header = TRUE,
              sep = ",", quote = "\"'",
              dec = ".",
              na.strings = "NaN",
              colClasses=c("factor","numeric","factor","numeric","factor","factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))
