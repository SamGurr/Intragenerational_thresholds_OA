#Author: Sam Gurr 
#Edited by: Sam Gurr
#Date Last Modified: 20191030
#Purpose: anaysis of TAC data (standard curve and calculations)
#See Readme file for details

rm(list=ls())

# Install packages if not already in your library-----------------------------------------------------------------------------------------------
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("lsmeans" %in% rownames(installed.packages()) == 'FALSE') install.packages('lsmeans') 
if ("ggpubr" %in% rownames(installed.packages()) == 'FALSE') install.packages('ggpubr') 
# Load packages and pacage version/date/import/depends info
library(dplyr) # Version 0.7.6, Packaged: 2018-06-27, Depends: R (>= 3.1.2)Imports: assertthat (>= 0.2.0), bindrcpp (>= 0.2.0.9000), glue (>=1.1.1), magrittr (>= 1.5), methods, pkgconfig (>= 2.0.1), R6(>= 2.2.2), Rcpp (>= 0.12.15), rlang (>= 0.2.0), tibble (>=1.3.1), tidyselect (>= 0.2.3), utils
library(lsmeans) # Version: 2.27-62, Date/Publication: 2018-05-11, Depends: methods, R (>= 3.2)
library(ggpubr)
#set working directory--------------------------------------------------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis") #set working

### CONICAL Seawater chemistry Data - Analysis, Graphs, Tables (APEX DATA) ####
TAC.Reference<-read.csv("Data/Phys.assays/TAC_data/Assay.Reference.csv", header=T, sep=",", na.string="NA", as.is=T) # reference to total antioxidant capacity assay - 96-well plate 
TP.Reference<-read.csv("Data/Phys.assays/Total.protein_data/Assay.Reference.csv", header=T, sep=",", na.string="NA", as.is=T) # reference to total protein assay - 96-well plate 
Master.Sample.Tracker <- read.csv("Data/Phys.assays/Master_table.csv", header=T, sep=",", na.string="NA", as.is=T) # upload the master tracker with the treatment and ID
colnames(Master.Sample.Tracker)[7] <- "ID" # change to "ID" to ease a merge later in the script

################################################################################################################################## #
################################################################################################################################## #
###############################################################  TOTAL PROTEIN ###################################################
##################################################### Pierce Rapid Gold BCA Protein Assay Kit #################################### #
################################################################################################################################## #

# FOR LOOP to for a cumulative dataset of all TAC runs ---------------------------------------------------------------------------------- #

# PRE for loop prep
path.TP.files<-"Data/Phys.assays/Total.protein_data/csv.files" #the location of all TP spec files (.csv format)
TP.file.names.full<-basename(list.files(path = path.TP.files, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders
TP.file.names.full # look at the names of the csv files you will call in the following for loop
# INSIDE for loop prep 
DF.TP <- data.frame(matrix(nrow = 96, ncol = 4))  # first create a master data table to bind in the for loop
colnames(DF.TP)<-c('Date', 'Run', 'Well', 'Absorbance.RAW') # column names, matches all data
DF.TP.TOTAL <- data.frame() # start dataframe 
# RUN FOR LOOP
for(i in 1:length(TP.file.names.full)) { # for every file name in list 
  TP.Data <-read.table(file.path(path.TP.files, TP.file.names.full[i]), 
                        header=T, sep=",", na.string="NA", 
                        fill = TRUE, as.is=TRUE) # reads in the data files of that name one by one
  DF.TP$Date <- substr(TP.file.names.full[i], 1,8) # date from file name
  DF.TP$Run <- substr(TP.file.names.full[i], 22,22) # Run number from file name
  DF.TP$Well <- TP.Data[,1] # row containing the well ID
  DF.TP$Absorbance.RAW <- TP.Data[,2] # raw absorbance value
  DF.loop <- data.frame(DF.TP) # assign the table to a different name
  DF.TP.TOTAL <- rbind(DF.TP.TOTAL,DF.loop) # bind the newly named table to a cumulative dataframe
  print(DF.TP.TOTAL) # print to monitor progress
  # script returns to next file in order and creates a new DF.TP to rename as DF.loop and bind to DF.TP.TOTAL wihtin the for loop
}
DF.TP.TOTAL # view cumulative data table

# Analysis Steps
# 1. Merge output with ID references  --------------------------------------------------------------------------------------------------------- #
# 2. Seperate datasets into 'initial' and 'post' copper reaction ------------------------------------------------------------------------------ #
# 3. Obtain median (where applies) of all 'initial' reads prior to Copper rxn   --------------------------------------------------------------- #
# 4. Merge the two datasets to simplify the correction for initial reads (post - initial = diff)  --------------------------------------------- #
# 5. Seperate Standards - assign theoretical concentration values (from Oxiselect protocol) --------------------------------------------------- #
# 6. Calculate UAE and CRE (UAE using y=mx+b from standard curve AND CRE is UAE*2189)
# 7. Merge with the master sheet and correct for homogenate volume, afdw, and total protein
# 8. Run models for treatment effect(s)


# 1. ------------------------------------------------------------------------------------------------------------------------------------------ #


TP.MERGE.Reference <- merge(DF.TP.TOTAL, TP.Reference, by = c('Date', 'Run', 'Well'))
TP.MERGE.Reference <- na.omit(TP.MERGE.Reference) # ommit NAs
TP.MERGE.Reference$Absorbance.RAW <- as.numeric(TP.MERGE.Reference$Absorbance.RAW) # make the absorbance values numeric (raw data from spec as characters)
tail(TP.MERGE.Reference)
# 2. ------------------------------------------------------------------------------------------------------------------------------------------ #

TP.Blanks.means <- TP.MERGE.Reference %>% dplyr::filter(ID == 'I') %>% 
  dplyr::group_by(Date, Run) %>% 
  dplyr::summarise(mean.blank = mean(Absorbance.RAW),
                   count.blank =n()) # call mean raw absorbance
TP.Blanks.means # view blanks


# 3. ------------------------------------------------------------------------------------------------------------------------------------------ #

TP.Standards.means <- TP.MERGE.Reference %>% dplyr::filter(Type == 'Standard') %>% 
     dplyr::group_by(Date, Run, ID) %>% 
     dplyr::summarise(mean.absorbance = mean(Absorbance.RAW),
                   count.standard =n()) # call mean raw absorbance
TP.Standards.means # view file


# 4. ------------------------------------------------------------------------------------------------------------------------------------------ #
TP.Standards.concentration <- data.frame(matrix(nrow = 9, ncol = 2)) # create dataframe with 9 rows (standards) and 2 columns (based on # of standards)
colnames(TP.Standards.concentration)<-c('ID', 'µg.mL_BCA') # assign column names
TP.Standards.concentration$ID <- c("A","B","C","D","E","F","G","H","I") # assign IDs (must be the same format as IDs in reference csv file)
TP.Standards.concentration$mg.mL_BCA <- c(2, 1.5, 1.0, 0.75, 0.5, 0.25, 0.125, 0.025, 0) # assign mg.mL_BCA (Pierce BCA protocol Rapid Gold Assay)

TP.Standards.MERGE.1 <- merge(TP.Standards.concentration, TP.Standards.means, by=c('ID')) # merge with the means by ID to add theoretical BCA values
TP.Standards.MERGE.FINAL <- merge(TP.Standards.MERGE.1, TP.Blanks.means, by=c('Date', 'Run')) # merge with the mean blanks to correct BELOW
TP.Standards.MERGE.FINAL$mean.absorbance.CORRECTED <- TP.Standards.MERGE.FINAL$mean.absorbance - TP.Standards.MERGE.FINAL$mean.blank # correct to mean blank absorbance as "mean.absorbance.CORRECTED"
# NOTE: max standard is 2.6290
# SAMPLES ANNNALYZED ON 11/04
TP.Standards.Nov <- TP.Standards.MERGE.FINAL %>% dplyr::filter(Date == '20191104')
TP.Standards.model.Nov <- lm(mean.absorbance.CORRECTED ~ mg.mL_BCA, data=TP.Standards.Nov) #runs a linear regression of mV as a function of temperature
TP.b.Nov <- summary(TP.Standards.model.Nov)$coefficients[1] # 0.04205469
TP.m.Nov <- summary(TP.Standards.model.Nov)$coefficients[2] # 1.334554
TP.R2.Nov<-summary(TP.Standards.model.Nov)$r.squared # 0.9967063
plot(mean.absorbance.CORRECTED ~ mg.mL_BCA, data=TP.Standards.Nov)
summary(lm(mean.absorbance.CORRECTED ~ mg.mL_BCA, data=TP.Standards.Nov)) # Adjusted R-squared:  0.9962 
abline(lm(mean.absorbance.CORRECTED ~ mg.mL_BCA, data=TP.Standards.Nov))
legend('topleft', legend = c((bquote(R^2 == .(format(TP.R2.Nov), digits = 5))), 
                             (bquote(slope == .(format(TP.m.Nov), digits = 5))),
                             (bquote(intercept == .(format(TP.b.Nov), digits = 5)))),  bty='n')
# SAMPLES ANNNALYZED ON 12/16
TP.Standards.Dec <- TP.Standards.MERGE.FINAL %>% dplyr::filter(Date == '20191216')
TP.Standards.model.Dec <- lm(mean.absorbance.CORRECTED ~ mg.mL_BCA, data=TP.Standards.Dec) #runs a linear regression of mV as a function of temperature
TP.b.Dec <- summary(TP.Standards.model.Dec)$coefficients[1] # 0.04445541
TP.m.Dec <- summary(TP.Standards.model.Dec)$coefficients[2] # 1.299334
TP.R2.Dec <-summary(TP.Standards.model.Dec)$r.squared #  0.9953494
plot(mean.absorbance.CORRECTED ~ mg.mL_BCA, data=TP.Standards.Dec)
summary(lm(mean.absorbance.CORRECTED ~ mg.mL_BCA, data=TP.Standards.Dec)) # Adjusted R-squared:  0.9947 
abline(lm(mean.absorbance.CORRECTED ~ mg.mL_BCA, data=TP.Standards.Dec))
legend('topleft', legend = c((bquote(R^2 == .(format(TP.R2.Dec), digits = 5))), 
                             (bquote(slope == .(format(TP.m.Dec), digits = 5))),
                             (bquote(intercept == .(format(TP.b.Dec), digits = 5)))),  bty='n')

# 6. ------------------------------------------------------------------------------------------------------------------------------------------ #


# remove standards from data, Calculate Total Protein
TP.Samples.DF <- TP.MERGE.Reference %>% dplyr::filter(!Type == "Standard") # remove standards
TP.MASTER.DF <- merge(TP.Samples.DF, TP.Blanks.means, by=c('Date')) # align blank values by date and run to correct absorbance values
tail(TP.MASTER.DF) # view newer data
TP.MASTER.DF$Absorbance.CORRECTED <- (TP.MASTER.DF$Absorbance.RAW -  TP.MASTER.DF$mean.blank) # subtract from the blank on the corresponding date (aligned from merge)


TP.MASTER.DF.Nov <-  TP.MASTER.DF %>% dplyr::filter(Date == '20191104') 
TP.MASTER.DF.Nov$mg.mL_BCA  <- ((TP.MASTER.DF.Nov$Absorbance.CORRECTED - TP.b.Nov)/TP.m.Nov) # Calculate µg.mL_BCA 
plot(Absorbance.CORRECTED ~ mg.mL_BCA, data=TP.MASTER.DF.Nov) # plot the data

TP.MASTER.DF.Dec <-  TP.MASTER.DF %>% dplyr::filter(Date == '20191216') 
TP.MASTER.DF.Dec$mg.mL_BCA  <- ((TP.MASTER.DF.Dec$Absorbance.CORRECTED - TP.b.Dec)/TP.m.Dec) # Calculate µg.mL_BCA 
plot(Absorbance.CORRECTED ~ mg.mL_BCA, data=TP.MASTER.DF.Dec) # plot the data


TP.MASTER.DF.2 <- rbind(TP.MASTER.DF.Nov,TP.MASTER.DF.Dec) # view total protein dataframe
TP.MASTER.DF.2 <- na.omit(TP.MASTER.DF.2) # removes both duplicate measurements of ID 1161 and ONE measurement of ID 205 (over the spectrum of spec)

# 7. ------------------------------------------------------------------------------------------------------------------------------------------ #


# Summarize by ID to obtain the median (assay completed in duplicates)
MASTER.TP.MEANS <- TP.MASTER.DF.2 %>%  
 # dplyr::filter(Absorbance.CORRECTED < 2.6290) %>% 
  dplyr::group_by(Date, ID) %>% # call column to summarize 
  dplyr::summarise(mean.mg.mL_BCA = mean(mg.mL_BCA),
                   count.mean.µg.mL_BCA = n())
MASTER.TP.MEANS # view the table of mean values (actually the median of duplicates - prints the single value for readings not duplicated, e.g. 205)

################################################################################################################################## #
################################################################################################################################## #
###############################################################  TOTAL ANTIOXIDANT CAPACITY  #####################################
################################################################### OXISELECT ASSSAY     ######################################### #
################################################################################################################################## #
# NOTE! this assay has both initial and post reaction absorbance readings  ####################################################### #

# FOR LOOP to for a cumulative dataset of all TAC runs ---------------------------------------------------------------------------------- #

# PRE for loop prep
path.TAC.files<-"Data/Phys.assays/TAC_data/csv.files" #the location of all TAC spec files (.csv format)
TAC.file.names.full<-basename(list.files(path = path.TAC.files, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders
TAC.file.names.full # look at the names of the csv files you will call in the following for loop
# INSIDE for loop prep 
DF.TAC <- data.frame(matrix(nrow = 96, ncol = 5))  # first create a master data table to bind in the for loop
colnames(DF.TAC)<-c('Date', 'Run', 'Well', 'Absorbance.RAW', 'Initial.Post') # column names, matches all data
DF.TAC.TOTAL <- data.frame() # start dataframe 
# RUN FOR LOOP
for(i in 1:length(TAC.file.names.full)) { # for every file name in list 
  TAC.Data <-read.table(file.path(path.TAC.files, TAC.file.names.full[i]), 
                        header=T, sep=",", na.string="NA", 
                        fill = TRUE, as.is=TRUE) # reads in the data files of that name one by one
  DF.TAC$Date <- substr(TAC.file.names.full[i], 1,8) # date from file name
  DF.TAC$Run <- substr(TAC.file.names.full[i], 14,14) # Run number from file name
  DF.TAC$Well <- TAC.Data[,1] # row containing the well ID
  DF.TAC$Initial.Post <- substr(TAC.file.names.full[i], 16,16) # first letter of 'initial' or 'post' in file name
  DF.TAC$Absorbance.RAW <- TAC.Data[,2] # raw absorbance value
  DF.loop <- data.frame(DF.TAC) # assign the table to a different name
  DF.TAC.TOTAL <- rbind(DF.TAC.TOTAL,DF.loop) # bind the newly named table to a cumulative dataframe
  print(DF.TAC.TOTAL) # print to monitor progress
  # script returns to next file in order and creates a new DF.TAC to rename as DF.loop and bind to DF.TAC.TOTAL wihtin the for loop
}
DF.TAC.TOTAL # view cumulative data table
tail(DF.TAC.TOTAL, n = 50)
# Analysis Steps
# 1. Merge output with ID references  --------------------------------------------------------------------------------------------------------- #
# 2. Seperate datasets into 'initial' and 'post' copper reaction ------------------------------------------------------------------------------ #
# 3. Obtain median (where applies) of all 'initial' reads prior to Copper rxn   --------------------------------------------------------------- #
# 4. Merge the two datasets to simplify the correction for initial reads (post - initial = diff)  --------------------------------------------- #
# 5. Seperate Standards - assign theoretical concentration values (from Oxiselect protocol) --------------------------------------------------- #
# 6. Calculate UAE and CRE (UAE using y=mx+b from standard curve AND CRE is UAE*2189)
# 7. Merge with the master sheet and correct for homogenate volume, afdw, and total protein
# 8. Run models for treatment effect(s)


# 1. ------------------------------------------------------------------------------------------------------------------------------------------ #


# merge with Reference and ommit NAs
TAC.MERGE.Reference <- merge(DF.TAC.TOTAL, TAC.Reference, by = c('Date','Run', 'Well')) # merge with the reference data
tail(TAC.MERGE.Reference, n = 50)
TAC.MERGE.Reference <- na.omit(TAC.MERGE.Reference) # ommit NAs
tail(TAC.MERGE.Reference, n = 50)

# 2 and 3. ------------------------------------------------------------------------------------------------------------------------------------ #


# INITIAL READS calculate median for duplicated values
INITIAL.Calc.DF <- TAC.MERGE.Reference %>%  dplyr::filter(Initial.Post %in% 'i')# %>%  # call only initial data
INITIAL.Calc.DF # view the table of mean values (actually the median of duplicates - prints the single value for readings not duplicated)
tail(INITIAL.Calc.DF, n = 50)

# POST READS (after the copper reaction)
POST.DF <- TAC.MERGE.Reference %>%  dplyr::filter(Initial.Post %in% 'p') # call only post data
tail(POST.DF, n = 50)

# 4. ------------------------------------------------------------------------------------------------------------------------------------------ #


# Merge the Initial with the Post data - initial medians will be duplicated where multiple reads were completed on the same individuals ID
MERGE.DF <- merge(INITIAL.Calc.DF, POST.DF, by = c("Date", "Well", "Run", "Reagent.Date", "Type"))
MERGE.DF$Absorbance.diff <- MERGE.DF$Absorbance.RAW.y - MERGE.DF$Absorbance.RAW.x 

MERGE.DF.1 <- MERGE.DF %>% dplyr::filter(Run == 1) 
run1less.than.stand <- MERGE.DF.1 %>%  dplyr::filter(Absorbance.diff < 0.03) # find the values that will be negative downstream (less absorbance than the standard at 0)
MERGE.DF.2 <- MERGE.DF %>% dplyr::filter(Run == 2) 
run2less.than.stand <- MERGE.DF.2 %>%  dplyr::filter(Absorbance.diff < 0.03) # find the values that will be negative downstream (less absorbance than the standard at 0)
MERGE.DF.3 <- MERGE.DF %>% dplyr::filter(Run == 3) 
run3less.than.stand <- MERGE.DF.3 %>%  dplyr::filter(Absorbance.diff < 0.03) # one negative value of -0.001 - make this zero in the following...

MERGE.DF$Absorbance.diff <- ifelse(MERGE.DF$Absorbance.diff < 0.0000, 0, MERGE.DF$Absorbance.diff) # make negative value ZERO (only one in the third run)

#MERGE.DF.d <-  MERGE.DF %>% dplyr::filter(ID.x == 1483)


# 5. ------------------------------------------------------------------------------------------------------------------------------------------ #
MASTER.TAC.DF <- MERGE.DF %>% dplyr::select(Date, Reagent.Date, ID.x, Type, Absorbance.diff) # select target columns
colnames(MASTER.TAC.DF)[3] <- "ID" # change column name to ID

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
# NOTE: highest absorbance of standard 0.4925 
# lowest standard absorbance is 0.0300
#Standards.MASTER.2 <- Standards.MASTER %>% filter(Stand.ABS.median > 0.0301) # remove the zero standard
#Standards.MASTER.2$Stand.ABS.median <- Standards.MASTER.2$Stand.ABS.median + 0.03 # add the absorbance of the zero standard to all other standards 
# Note: this will give the accurate calculattion from the regression to calculation UAE based on all the data added by 0.03 (the zero standard absorbance)


# Standard Curves
# Standard curve - OCTOBER run
Standards.Oct <- Standards.MASTER %>% dplyr::filter(Date %in% '20191030')
Standards.model.Oct <-lm(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.Oct) #runs a linear regression of mV as a function of temperature
plot(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.Oct) # high standard is off

Standards.Oct.2 <- Standards.Oct %>% dplyr::filter(mM.Uric.Acid < 1.0) # remove the high standard for better fit
Standards.model.Oct.2 <-lm(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.Oct.2) #runs a linear regression of mV as a function of temperature
plot(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.Oct.2) # better  fit
summary(lm(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.Oct.2)) # Adjusted R-squared:  0.9906 
abline(lm(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.Oct.2))
#b.Oct <- 0.03#summary(Standards.model)$coefficients[1]
b.Oct <- summary(Standards.model.Oct.2)$coefficients[1] # 0.03590498
m.Oct <- summary(Standards.model.Oct.2)$coefficients[2] # 0.5781389
R2.Oct<-summary(Standards.model.Oct.2)$r.squared # 0.9917731
legend('topleft', legend = bquote(R^2 == .(format(R2.Oct, digits = 3))), bty='n',)


# Standard curve - Dec run
Standards.Dec <- Standards.MASTER %>% dplyr::filter(Date %in% '20191210')
Standards.model.Dec <-lm(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.Dec) #runs a linear regression of mV as a function of temperature
plot(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.Dec) # high standard is off
abline(lm(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.Dec))
summary(Standards.model.Dec) # 0.9991
b.Dec <- summary(Standards.model.Dec)$coefficients[1] # 0.03982677
m.Dec <- summary(Standards.model.Dec)$coefficients[2] # 1.272379
R2.Dec<-summary(Standards.model.Dec)$r.squared # 0.9991627
legend('topleft', legend = bquote(R^2 == .(format(R2.Dec, digits = 3))), bty='n',)


# 6. ------------------------------------------------------------------------------------------------------------------------------------------ #

# Master Data Frame 
zero.standards <- Standards.MASTER %>%  filter(mM.Uric.Acid == 0) # 0.0300 and 0.0335
MASTER.TAC.DF.Oct <- MASTER.TAC.DF %>%  filter(Date == 20191030) # october samples 
MASTER.TAC.DF.Oct$zero.stand <- 0.0300
MASTER.TAC.DF.Dec <- MASTER.TAC.DF %>%  filter(Date %in% c('20191210', '20191212')) # december samples
MASTER.TAC.DF.Dec$zero.stand <- 0.0335
MASTER.TAC.DF.2 <- rbind(MASTER.TAC.DF.Oct, MASTER.TAC.DF.Dec)
MASTER.TAC.DF.2$Absorbance.diff.add.zero.stand <- MASTER.TAC.DF.2$Absorbance.diff + MASTER.TAC.DF.2$zero.stand # new column add the ZERO absorbance of the standard at zero
MASTER.TAC.DF.means <- MASTER.TAC.DF.2 %>% group_by(Date,ID,Type) %>% 
  dplyr::summarise(Absorbance.mean = mean(Absorbance.diff.add.zero.stand))# summarize data table of the mean values for the absorbance corrected for the zero standard

# remove standards from data, Calculate UAE and CRE
MASTER.TAC.DF.SAMPLES <- MASTER.TAC.DF.means %>% dplyr::filter(!Type == "Standard") # remove standards
SAMPLES.under.standards <- MASTER.TAC.DF.SAMPLES %>% dplyr::filter(Absorbance.mean < 0.0300)

MASTER.TAC.DF.SAMPLES.Oct <- MASTER.TAC.DF.SAMPLES %>%  filter(Date == 20191030) # october samples 
#MASTER.TAC.DF.SAMPLES.Oct$mM.Uric.Acid.Equivalents <- (MASTER.TAC.DF.SAMPLES.Oct$Absorbance.mean - b.Oct)/m.Oct # Calculate URIC ACID EQUIVALENTS (UAE) - weighed low mass of Cu standard - may have confoundded our calib curve
MASTER.TAC.DF.SAMPLES.Oct$mM.Uric.Acid.Equivalents <- (MASTER.TAC.DF.SAMPLES.Oct$Absorbance.mean - b.Dec)/m.Dec # based on calibration curve in december run Calculate URIC ACID EQUIVALENTS (UAE)
MASTER.TAC.DF.SAMPLES.Dec  <- MASTER.TAC.DF.SAMPLES %>%  filter(Date %in% c('20191210', '20191212')) # december samples
MASTER.TAC.DF.SAMPLES.Dec$mM.Uric.Acid.Equivalents <- (MASTER.TAC.DF.SAMPLES.Dec$Absorbance.mean - b.Dec)/m.Dec # Calculate URIC ACID EQUIVALENTS (UAE)

MASTER.TAC.DF.SAMPLES.2 <- rbind(MASTER.TAC.DF.SAMPLES.Oct, MASTER.TAC.DF.SAMPLES.Dec)

plot(Absorbance.mean ~ mM.Uric.Acid.Equivalents, data=MASTER.TAC.DF.SAMPLES.2) # plot the data
summary(lm(Absorbance.mean ~ mM.Uric.Acid.Equivalents, data=MASTER.TAC.DF.SAMPLES.2))

MASTER.TAC.DF.SAMPLES.2$µM.Copper.Reducing.Equivalents <- MASTER.TAC.DF.SAMPLES.2$mM.Uric.Acid.Equivalents*2189 # Calculate µM Copper Reducing Equivlents (CRE)

# Summarize by ID to obtain the median (assay completed in duplicates)
MASTER.TAC.MEANS <- MASTER.TAC.DF.SAMPLES.2 %>%  
  dplyr::group_by(Date, ID) %>% # call column to summarize 
  dplyr::summarise(mean.UAE = mean(mM.Uric.Acid.Equivalents), 
                   mean.CRE = mean(µM.Copper.Reducing.Equivalents)) 
MASTER.TAC.MEANS # view the table of mean values (actually the median of duplicates - prints the single value for readings not duplicated)
tail(MASTER.TAC.MEANS)


################################################################################################################################## #
################################################################################################################################## #
###############################################################  MERGE FOR FINAL DATASET  ##########################################
################################################################### OXISELECT ASSSAY     ######################################### #
################################################################################################################################## #
MASTER.TAC.MEANS # calculated total antioxidant capacity
sum(complete.cases(MASTER.TAC.MEANS)) # 144 - still have ~ 20 samples to run
MASTER.TP.MEANS # calculated total antioxidant capacity
sum(complete.cases(MASTER.TP.MEANS)) # 161 - three less than total master sheet due to bad spec measurements
Master.Sample.Tracker # master sheet of all data (total homogenized volume, AFDW, volumes for each assay, dolution information, etc.)
sum(complete.cases(Master.Sample.Tracker$ID))  # 164

# Merge with the master sheet (containing Treatment IDs and Correction factors)
TAC.TP.MERGE <- merge(MASTER.TP.MEANS, MASTER.TAC.MEANS, by = 'ID')
Master.Sample.Tracker$ID <- Master.Sample.Tracker$NEW.Tube.ID
MASTER.PHYS.ASSAY <- merge(TAC.TP.MERGE, Master.Sample.Tracker, by = 'ID')
# correct for homogenate volume as "TOTAL"
# TOTAL.total.protein mg/mL - (total homogenate vol / vol aliquot for protein assay) * (protein measured in aliquot*dilution factor)
sapply(MASTER.PHYS.ASSAY, class) # check the class - occasionally variables that are numeric are characters
MASTER.PHYS.ASSAY$µl.Total.Homog.Volume <- as.numeric(MASTER.PHYS.ASSAY$µl.Total.Homog.Volume) # convert to numeric
MASTER.PHYS.ASSAY$µl.Total.protein  <- as.numeric(MASTER.PHYS.ASSAY$µl.Total.protein)  # convert to numeric
MASTER.PHYS.ASSAY$mean.mg.mL_BCA # = mg / ml raw signal from the calibration curve; 
MASTER.PHYS.ASSAY$TOTAL_PROTEIN_per_well_corrected <- (MASTER.PHYS.ASSAY$mean.mg.mL_BCA*MASTER.PHYS.ASSAY$Dil.factor.Total.Protein) # mg  protein ONLY IN THE 20ul sample (account for dilution from NaOH and HCl; 20 ul used for the Rapid Gold assay
# NOTE: TOTAL_PROTEIN_per_SAMPLE is the mg of protein in the 20 µl sample for the rapid gold assay kit

# TOTAL.mean.CRE - µm cOPPER reDUCING eQUIVALENTS MG PROTEIN -1 (IN THE 20 UL SAMPLE)
#MASTER.PHYS.ASSAY$µl.TAC <- as.numeric(MASTER.PHYS.ASSAY$µl.TAC) # convert to numericvolume aliquot (~50 - 100 ul; in spreadsheet) from total homogenate volume (HV = PBS + whole geoduck)
MASTER.PHYS.ASSAY$mean.CRE <- as.numeric(MASTER.PHYS.ASSAY$mean.CRE) # uM concentration NOT DILUTED convert to numeric

MASTER.PHYS.ASSAY$Dil.factor.TAC <- as.numeric(MASTER.PHYS.ASSAY$Dil.factor.TAC) # convert to numeric dilution factor for a few samples that were additionally diluted prior to NaOH and HCl
MASTER.PHYS.ASSAY$µM.CRE.mg.protein<- (MASTER.PHYS.ASSAY$mean.CRE*MASTER.PHYS.ASSAY$Dil.factor.TAC) / (MASTER.PHYS.ASSAY$TOTAL_PROTEIN_per_well_corrected) # µM CRE mg protein = concentration of CRE in the 20 µl sample / protein in 20 ul sample

# Total protein WHOLE ANIMALS corrected for AFDW (AS mgProtein.mgAFDW)
MASTER.PHYS.ASSAY$TOTAL.PROTEIN_whole.animal.mg <- (MASTER.PHYS.ASSAY$µl.Total.Homog.Volume/20)*MASTER.PHYS.ASSAY$TOTAL_PROTEIN_per_well_corrected
MASTER.PHYS.ASSAY$mgProtein.mgAFDW <- MASTER.PHYS.ASSAY$TOTAL.PROTEIN_whole.animal.mg /  MASTER.PHYS.ASSAY$mgTOTAL_AFDW

PHYS.DATA <- MASTER.PHYS.ASSAY %>% dplyr::select(Date.fixed, ID, Treatment, Exprmt.Day, Tank.ID, 
                                                 mgTOTAL_AFDW,
                                                 mgProtein.mgAFDW,
                                                 µM.CRE.mg.protein)

#PHYS.DATA <- na.omit(PHYS.DATA)

count.days <- PHYS.DATA %>%  dplyr::group_by(Exprmt.Day) %>% 
  dplyr::summarise(count_replicated.by.day =n()) # day 1, 12; day 4, 12; day 7, 35; day 15,  16, day 18, 15; day 21, 49

count.treat.D7 <- PHYS.DATA %>%  dplyr::filter(Exprmt.Day %in% 'Day7') %>% 
  dplyr::group_by(Treatment) %>% 
  dplyr::summarise(count_replicated.by.day =n())
count.treat.D7

count.treat.D21 <- PHYS.DATA %>%  dplyr::filter(Exprmt.Day %in% 'DAY21') %>% 
  dplyr::group_by(Treatment) %>% 
  dplyr::summarise(count_replicated.by.day =n())
count.treat.D21


# write out tables
write.table(PHYS.DATA, file="C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/Output/Phys.Assay.Table.csv", sep=",", row.names = FALSE) #save data to output file

