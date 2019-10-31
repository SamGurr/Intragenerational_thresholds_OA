#Author: Sam Gurr 
#Edited by: Sam Gurr
#Date Last Modified: 20191030
#Purpose: anaysis of TAC data (standard curve and calculations)
#See Readme file for details

rm(list=ls())

# Install packages if not already in your library
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
# Load packages and pacage version/date/import/depends info
library(dplyr) # Version 0.7.6, Packaged: 2018-06-27, Depends: R (>= 3.1.2)Imports: assertthat (>= 0.2.0), bindrcpp (>= 0.2.0.9000), glue (>=1.1.1), magrittr (>= 1.5), methods, pkgconfig (>= 2.0.1), R6(>= 2.2.2), Rcpp (>= 0.12.15), rlang (>= 0.2.0), tibble (>=1.3.1), tidyselect (>= 0.2.3), utils
library(lsmeans)        # Version: 2.27-62, Date/Publication: 2018-05-11, Depends: methods, R (>= 3.2)
#set working directory--------------------------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis") #set working

### CONICAL Seawater chemistry Data - Analysis, Graphs, Tables (APEX DATA) ####
Reference<-read.csv("Data/Phys.assays/TAC_data/Assay.Reference.csv", header=T, sep=",", na.string="NA", as.is=T) 
Master.Sample.Tracker <- read.csv("Data/Phys.assays/Master_table.csv", header=T, sep=",", na.string="NA", as.is=T) # upload the master tracker with the treatment and ID
Master.Sample.Tracker <- Master.Sample.Tracker %>% dplyr::select(Date.fixed, Exprmt.Day, Tank.ID,
                                                                 Treatment, NEW.Tube.ID, µl.Total.Homog.Volume, µl.TAC) # slect only necessary columns
colnames(Master.Sample.Tracker)[5] <- "ID" # cahnge to "ID" to ease a merge later in the script
path.p<-"Data/Phys.assays/TAC_data/csv.files" #the location of all TAC spec files (.csv format)
file.names.full<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders
file.names.full # look at the names of the csv files you will call in the following for loop
# file.names <- file.names.full[c(25:28)] # call the files you want to analyze and rbind to the current cumunaltive file
# file.names # look at the names of the csv files you will call in the following for loop



# FOR LOOP to for a cumulative dataset of all TAC runs ---------------------------------------------------------------------------------- #



# bind all data with the same Reagent.Date; name dataset DATE_TAC.RUNS
DF.TAC <- data.frame(matrix(nrow = 96, ncol = 5))  # first create a master data table to bind in the for loop
colnames(DF.TAC)<-c('Date', 'Run', 'Well', 'Absorbance.RAW', 'Initial.Post') # column names of matches all data
DF.TAC.TOTAL <- data.frame() # start dataframe 

for(i in 1:length(file.names.full)) { # for every file name in list 
  TAC.Data <-read.table(file.path(path.p,file.names.full[i]), 
                        header=T, sep=",", na.string="NA", 
                        fill = TRUE, as.is=TRUE) # reads in the data files of that name one by one
  DF.TAC$Date <- substr(file.names.full[i], 1,8) # date from file name
  DF.TAC$Run <- substr(file.names.full[i], 14,14) # Run number from file name
  DF.TAC$Well <- TAC.Data[,1] # row containing the well ID
  DF.TAC$Initial.Post <- substr(file.names.full[i], 16,16) # first letter of 'initial' or 'post' in file name
  DF.TAC$Absorbance.RAW <- TAC.Data[,2] # raw absorbance value
  DF.loop <- data.frame(DF.TAC) # assign the table to a different name
  DF.TAC.TOTAL <- rbind(DF.TAC.TOTAL,DF.loop) # bind the newly named table to a cumulative dataframe
  print(DF.TAC.TOTAL) # print to monitor progress
  # script returns to next file in order and creates a new DF.TAC to rename as DF.loop and bind to DF.TAC.TOTAL wihtin the for loop
}



# Analysis Steps
# 1. Merge output with ID references  --------------------------------------------------------------------------------------------------------- #
# 2. Seperate datasets into 'initial' and 'post' copper reaction ------------------------------------------------------------------------------ #
# 3. Obtain median (where applies) of all 'initial' reads prior to Copper rxn   --------------------------------------------------------------- #
# 4. Merge the two datasets to simplify the correction for initial reads (post - initial = diff)  --------------------------------------------- #
# 5. Seperate Standards - assign theoretical concentration values (from Oxiselect protocol) --------------------------------------------------- #




# merge with Reference and ommit NAs
REF.MERGE.TAC.DF <- merge(DF.TAC.TOTAL, Reference, by = c('Date','Run', 'Well')) # merge with the reference data
REF.MERGE.TAC.DF <- na.omit(REF.MERGE.TAC.DF) # ommit NAs
# INITIAL READS calculate median for duplicated values
INITIAL.Calc.DF <- REF.MERGE.TAC.DF %>%  dplyr::filter(Initial.Post %in% 'i') %>%  # call only initial data
  dplyr::group_by(Date, Run, ID, Initial.Post) %>% # call column to summarize 
  dplyr::summarise(initial_mean_ABS = mean(Absorbance.RAW), # call mean raw absorbance
                   initial_count_ABS =n()) %>%  # call for the count of reads
  dplyr::arrange(desc(initial_count_ABS)) # makes table in descending order 
INITIAL.Calc.DF # view the table of mean values (actually the median of duplicates - prints the single value for readings not duplicated)
# POST READS (after the copper reaction)
POST.DF <- REF.MERGE.TAC.DF %>%  dplyr::filter(Initial.Post %in% 'p') # call only post data
# Merge the Initial with the Post data - initial medians will be duplicated where multiple reads were completed on the same individuals ID
MERGE.DF <- merge(INITIAL.Calc.DF, POST.DF, by = c("Date", "ID", "Run"))
MERGE.DF$Absorbance.diff <- MERGE.DF$Absorbance.RAW - MERGE.DF$initial_mean_ABS
# Master Data Frame
MASTER.TAC.DF <- MERGE.DF %>% select(Date, Reagent.Date, ID, Type, Absorbance.diff)

# Standards
Stand.reference <- data.frame(matrix(nrow = 10, ncol = 2)) # create dataframe with 10 rwos (standards) and 2 columnds (ID and mM Uric Acid)
colnames(Stand.reference)<-c('ID', 'mM.Uric.Acid') # assign column names
Stand.reference$ID <- c("S.1","S.2","S.3","S.4","S.5","S.6","S.7","S.8","S.9","S.10") # assign IDs (must be the same format as IDs in reference csv file)
Stand.reference$mM.Uric.Acid <- c(1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.0156, 0.0078, 0.0039, 0.0) # assign mM uric acid (Oxiselect protocol)
Standards <- MASTER.TAC.DF %>% dplyr::filter(Type %in% "Standard") # call all rows that are standards
Standards.Median <- Standards %>% 
  dplyr::group_by(Date, ID) %>% # call column to summarize 
  dplyr::summarise(Stand.ABS.median = mean(Absorbance.diff)) # call mean (median)
Standards.Median # view the table of mean values (actually the median of duplicates - all corrected for the median of the 'initial'
Standards.MASTER <- merge(Standards.Median, Stand.reference, by = "ID")

Standards.MASTER <- Standards.MASTER %>% dplyr::filter(mM.Uric.Acid < 1.0)
Standards.model <-lm(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.MASTER) #runs a linear regression of mV as a function of temperature
b <- summary(Standards.model)$coefficients[1]
m <- summary(Standards.model)$coefficients[2]
R2<-summary(Standards.model)$r.squared
plot(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.MASTER)
summary(lm(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.MASTER))
abline(lm(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.MASTER))
legend('topleft', legend = bquote(R^2 == .(format(R2, digits = 3))), bty='n')


MASTER.TAC.DF.SAMPLES <- MASTER.TAC.DF %>% filter(!Type == "Standard") # remove standards
MASTER.TAC.DF.SAMPLES$mM.Uric.Acid.Equivalents <- (MASTER.TAC.DF.SAMPLES$Absorbance.diff - b)/m # Calculate URIC ACID EQUIVALENTS (UAE)
plot(Absorbance.diff ~ mM.Uric.Acid.Equivalents, data=MASTER.TAC.DF.SAMPLES) # plot the data
MASTER.TAC.DF.SAMPLES$µM.Copper.Reducing.Equivalents <- MASTER.TAC.DF.SAMPLES$mM.Uric.Acid.Equivalents*2189 # Calculate µM Copper Reducing Equivlents (CRE)


MASTER.TAC.MEANS <- MASTER.TAC.DF.SAMPLES %>%  
  dplyr::group_by(Date, ID) %>% # call column to summarize 
  dplyr::summarise(mean.UAE = mean(mM.Uric.Acid.Equivalents), 
                   mean.CRE = mean(µM.Copper.Reducing.Equivalents)) 
MASTER.TAC.MEANS # view the table of mean values (actually the median of duplicates - prints the single value for readings not duplicated)
# POST READS (after the copper reaction)
MASTER.TAC.FINAL <- merge(MASTER.TAC.MEANS, Master.Sample.Tracker, by = "ID")
MASTER.TAC.FINAL$µl.Total.Homog.Volume <- as.numeric(MASTER.TAC.FINAL$µl.Total.Homog.Volume)
MASTER.TAC.FINAL$µl.TAC <- as.numeric(MASTER.TAC.FINAL$µl.TAC)
MASTER.TAC.FINAL$mean.UAE <- as.numeric(MASTER.TAC.FINAL$mean.UAE)
MASTER.TAC.FINAL$mean.CRE <- as.numeric(MASTER.TAC.FINAL$mean.CRE)
MASTER.TAC.FINAL$mean.UAE.volcorrect <- MASTER.TAC.FINAL$mean.UAE * ((MASTER.TAC.FINAL$µl.Total.Homog.Volume)/(MASTER.TAC.FINAL$µl.Total.Homog.Volume-MASTER.TAC.FINAL$µl.TAC))
MASTER.TAC.FINAL$mean.CRE.volcorrect <- MASTER.TAC.FINAL$mean.CRE * ((MASTER.TAC.FINAL$µl.Total.Homog.Volume)/(MASTER.TAC.FINAL$µl.Total.Homog.Volume-MASTER.TAC.FINAL$µl.TAC))

# Run exploratory models
Data.D.1.4.7 <- MASTER.TAC.FINAL %>% filter(Date.fixed < 20190801)
Model_1_TAC <- aov(mean.CRE.volcorrect~Treatment, data = Data.D.1.4.7)
summary(Model_1_TAC)
# post-hoc
D.1.4.7.TAC <- lsmeans(Model_1_TAC, pairwise ~  Treatment)# pariwise Tukey Post-hoc test between repeated treatments
D.1.4.7.TAC # view post hoc summary
# HIGHER ANTIOXIDANTS EXPRESSED IN THE ELEVATED HISTORY ANIMALS
Data.D.15.18.21 <- MASTER.TAC.FINAL %>% filter(Date.fixed > 20190807)
Model_2_TAC <- aov(mean.CRE.volcorrect~Treatment, data = Data.D.15.18.21)
summary(Model_2_TAC)
# post-hoc
D.1.4.7.TAC <- lsmeans(Model_1_TAC, pairwise ~  Treatment)# pariwise Tukey Post-hoc test between repeated treatments
D.1.4.7.TAC # view post hoc summary