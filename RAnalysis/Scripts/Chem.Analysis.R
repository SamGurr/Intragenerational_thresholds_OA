# Intragenerational_thresholds (github repo)
#Data published in ___
#Title: 
#Contact: Sam Gurr @ samuel_gurr@uri.edu
#Supported by: FFAR
#See Readme file for details on data files and metadata

rm(list=ls()) # removes all prior objects

#R Version: R version 3.3.1
#RStudio Version: 1.0.44
######Read in required libraries#####
# Install packages if not already in your library

# Load packages and pacage version/date/import/depends info
library(dplyr)          # Version 0.7.6, Packaged: 2018-06-27, Depends: R (>= 3.1.2)Imports: assertthat (>= 0.2.0), bindrcpp (>= 0.2.0.9000), glue (>=1.1.1), magrittr (>= 1.5), methods, pkgconfig (>= 2.0.1), R6(>= 2.2.2), Rcpp (>= 0.12.15), rlang (>= 0.2.0), tibble (>=1.3.1), tidyselect (>= 0.2.3), utils
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')

#####Required Data files#####

#CRM_TA_Data.csv
#Daily_Temp_pH_Sal.csv
#~/MyProjects/BioMin_HIS/RAnalysis/Data/pH_Calibration_Files

#############################################################
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/") #set working directory
#setwd("~/MyProjects/Geoduck_Conditioning/RAnalysis/Data") #set working directory
mainDir<-'C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/' #set main directory
# mainDir<-'~/MyProjects/Geoduck_Conditioning/RAnalysis/' #set main directory
#############################################################

#call cumulative spreadsheet of discrete seawater chemistry
chem<-read.csv("Output/Seawater_chemistry_table_Output_All.csv", header=T, sep=",", na.string="NA", as.is=T) 
chem # view table 

short.table <- chem %>%  # call full table of carbonate chemistry
  distinct(Date, Treatment, pH, pCO2, TA, Aragonite.Sat) # call all unique values for target carbon chemistry
short.table <- as.data.frame(short.table) # call as data frame

date.table <- data.frame(unique(short.table$Date)) # call data frame of unique date value to loop
date.table # view table of dates to call in for loops



# ---------------------ANOVA LOOP for treatment efffects on carbonate chemistry---------------------------------


aov.results.table <- data.frame(matrix(nrow = 1,ncol = 17)) # template table to output anova results
colnames(aov.results.table) <-c('Date', 'pval.pH','pval.pCO2', 'pval.TA', 'pval.Aragonite.Sat',
                                'df.pH','df.pCO2', 'df.TA', 'df.Aragonite.Sat',
                                'df.Res.pH','df.Res.pCO2', 'df.Res.TA', 'df.Res.Aragonite.Sat',
                                'F.pH','F.pCO2', 'F.TA', 'F.Aragonite.Sat') # names for columns in the for loop
aov.results.FINAL <- data.frame()

for(i in 1:nrow(date.table)) {
  each.day.table <- short.table %>% 
    filter(Date == date.table[i,1])
  #ph aov 
  pH.anova.table <- summary(aov(lm(pH ~ Treatment, data = each.day.table)))
  x.pH <- as.data.frame(pH.anova.table[[1]])
  #pCO2 aov 
  pCO2.anova.table <- summary(aov(lm(pCO2 ~ Treatment, data = each.day.table)))
  x.pCO2 <- as.data.frame(pCO2.anova.table[[1]])
  #TA aov 
  TA.anova.table <- summary(aov(lm(TA ~ Treatment, data = each.day.table)))
  x.TA <- as.data.frame(TA.anova.table[[1]])
  #aragonite sat aov 
  Aragonite.Sat.anova.table <- summary(aov(lm(Aragonite.Sat ~ Treatment, data = each.day.table)))
  x.Aragonite.Sat <- as.data.frame(Aragonite.Sat.anova.table[[1]])
  # date
  aov.results.table$Date <- date.table[i,1] # all files have date in the form of yyyymmdd at the start of each csv name
  # p.values call = table[1,5]
  aov.results.table$pval.pH <- x.pH[1,5] 
  aov.results.table$pval.pCO2 <- x.pCO2[1,5]
  aov.results.table$pval.TA <- x.TA[1,5]
  aov.results.table$pval.Aragonite.Sat <- x.Aragonite.Sat[1,5]
  # df treatment  call = table[1,1]
  aov.results.table$df.pH <- x.pH[1,1]
  aov.results.table$df.pCO2 <- x.pCO2[1,1]
  aov.results.table$df.TA <- x.TA[1,1]
  aov.results.table$df.Aragonite.Sat <- x.Aragonite.Sat[1,1]
  # F call = table[1,4]
  aov.results.table$df.Res.pH <- x.pH[1,4]
  aov.results.table$df.Res.pCO2 <- x.pCO2[1,4]
  aov.results.table$df.Res.TA <- x.TA[1,4]
  aov.results.table$df.Res.Aragonite.Sat <- x.Aragonite.Sat[1,4]
  # df residuals  call = table[2,1]
  aov.results.table$F.pH <- x.pH[2,1]
  aov.results.table$F.pCO2 <- x.pCO2[2,1]
  aov.results.table$F.TA <- x.TA[2,1]
  aov.results.table$F.Aragonite.Sat <- x.Aragonite.Sat[2,1]
  
  aov.df <- data.frame(aov.results.table) # name dataframe for this single row
  aov.results.FINAL <- rbind(aov.results.FINAL,aov.df) #bind to a cumulative list dataframe
  print(aov.results.FINAL) # print final data table
}

aov.results.FINAL # view table
# write a table for anova results
write.table(aov.results.FINAL,"C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/Output/Carb.chem.anova.table.csv",sep=",", row.names=FALSE)  # write out to the path names outputNAME



# ---------------- Tables of mean ± st. error by date and summary -------------------------------------



for(i in 1:nrow(date.table)) {
  chem.by.date <- chem %>% 
    filter(Date == date.table[i,1])

chem.long <- melt(chem.by.date, id.vars=c("Date", "Tank", "Treatment")) # uses tidyr to make a long table from wide

Sum.Treatments <- ddply(chem.long, c("Date", "Treatment", "variable"), summarise, #Calculate descriptive stats by Treatment and exposure 
                        N = length(na.omit(value)), #count the sample size removing NA
                        mean = mean(value), #calculate average 
                        sem = sd(value)/sqrt(N)) # calcualte the standard error of the mean

#Exposure1.chem <-subset(Sum.Treatments, Date == i) #separate out exposure 1
Sum.long <- reshape(Sum.Treatments, idvar="Treatment", direction="wide", timevar = "variable", drop = c("Date","N")) #reshape data format for table layout
write.csv(Sum.long, paste0("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/Output/chem.tables/", date.table[i,1],".csv"), row.names = FALSE)
} # end of inside for loop


