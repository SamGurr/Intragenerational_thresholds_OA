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
TP.Standards.model <-lm(mean.absorbance.CORRECTED ~ mg.mL_BCA, data=TP.Standards.MERGE.FINAL) #runs a linear regression of mV as a function of temperature
TP.b <- summary(TP.Standards.model)$coefficients[1]
TP.m <- summary(TP.Standards.model)$coefficients[2]
TP.R2<-summary(TP.Standards.model)$r.squared
plot(mean.absorbance.CORRECTED ~ mg.mL_BCA, data=TP.Standards.MERGE.FINAL)
summary(lm(mean.absorbance.CORRECTED ~ mg.mL_BCA, data=TP.Standards.MERGE.FINAL))
abline(lm(mean.absorbance.CORRECTED ~ mg.mL_BCA, data=TP.Standards.MERGE.FINAL))
legend('topleft', legend = c((bquote(R^2 == .(format(TP.R2), digits = 5))), 
                             (bquote(slope == .(format(TP.m), digits = 5))),
                             (bquote(intercept == .(format(TP.b), digits = 5)))),  bty='n')


# 6. ------------------------------------------------------------------------------------------------------------------------------------------ #


# remove standards from data, Calculate Total Protein
TP.Samples.DF <- TP.MERGE.Reference %>% dplyr::filter(!Type == "Standard") # remove standards
TP.MASTER.DF <- merge(TP.Samples.DF, TP.Blanks.means, by=c('Date')) # align blank values by date and run to correct absorbance values
TP.MASTER.DF$Absorbance.CORRECTED <- (TP.MASTER.DF$Absorbance.RAW -  TP.MASTER.DF$mean.blank) # subtract from the blank on the corresponding date (aligned from merge)
TP.MASTER.DF$mg.mL_BCA  <- ((TP.MASTER.DF$Absorbance.CORRECTED - TP.b)/TP.m) # Calculate µg.mL_BCA 
plot(Absorbance.CORRECTED ~ mg.mL_BCA, data=TP.MASTER.DF) # plot the data

TP.MASTER.DF # view total protein dataframe
TP.MASTER.DF <- na.omit(TP.MASTER.DF) # removes both duplicate measurements of ID 1161 and ONE measurement of ID 205 (over the spectrum of spec)

# 7. ------------------------------------------------------------------------------------------------------------------------------------------ #


# Summarize by ID to obtain the median (assay completed in duplicates)
MASTER.TP.MEANS <- TP.MASTER.DF %>%  
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
TAC.MERGE.Reference <- na.omit(TAC.MERGE.Reference) # ommit NAs


# 2 and 3. ------------------------------------------------------------------------------------------------------------------------------------ #


# INITIAL READS calculate median for duplicated values
INITIAL.Calc.DF <- TAC.MERGE.Reference %>%  dplyr::filter(Initial.Post %in% 'i')# %>%  # call only initial data
INITIAL.Calc.DF # view the table of mean values (actually the median of duplicates - prints the single value for readings not duplicated)

# POST READS (after the copper reaction)
POST.DF <- TAC.MERGE.Reference %>%  dplyr::filter(Initial.Post %in% 'p') # call only post data
POST.DF

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

MERGE.DF.d <-  MERGE.DF %>% dplyr::filter(ID.x == 1483)

# Master Data Frame 
MASTER.TAC.DF <- MERGE.DF %>% dplyr::select(Date, Reagent.Date, ID.x, Type, Absorbance.diff) # select target columns
colnames(MASTER.TAC.DF)[3] <- "ID" # change column name to ID
MASTER.TAC.DF$Absorbance.diff.add.zero.stand <- MASTER.TAC.DF$Absorbance.diff + 0.03 # new column add the ZERO absorbance of the standard at zero
MASTER.TAC.DF.means <- MASTER.TAC.DF %>% group_by(Date,ID,Type) %>% 
  dplyr::summarise(Absorbance.mean = mean(Absorbance.diff.add.zero.stand))# summarize data table of the mean values for the absorbance corrected for the zero standard


# 5. ------------------------------------------------------------------------------------------------------------------------------------------ #


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

Standards.MASTER <- Standards.MASTER %>% dplyr::filter(mM.Uric.Acid < 1.0)
Standards.model <-lm(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.MASTER) #runs a linear regression of mV as a function of temperature
b <- 0.03#summary(Standards.model)$coefficients[1]
m <- summary(Standards.model)$coefficients[2]
R2<-summary(Standards.model)$r.squared
plot(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.MASTER)
summary(lm(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.MASTER))
abline(lm(Stand.ABS.median ~ mM.Uric.Acid, data=Standards.MASTER))
legend('topleft', legend = bquote(R^2 == .(format(R2, digits = 3))), bty='n')


# 6. ------------------------------------------------------------------------------------------------------------------------------------------ #


# remove standards from data, Calculate UAE and CRE
MASTER.TAC.DF.SAMPLES <- MASTER.TAC.DF.means %>% dplyr::filter(!Type == "Standard") # remove standards
SAMPLES.under.standards <- MASTER.TAC.DF.SAMPLES %>% dplyr::filter(Absorbance.mean < 0.0300)
MASTER.TAC.DF.SAMPLES$mM.Uric.Acid.Equivalents <- (MASTER.TAC.DF.SAMPLES$Absorbance.mean - b)/m # Calculate URIC ACID EQUIVALENTS (UAE)

plot(Absorbance.mean ~ mM.Uric.Acid.Equivalents, data=MASTER.TAC.DF.SAMPLES) # plot the data
summary(lm(Absorbance.mean ~ mM.Uric.Acid.Equivalents, data=MASTER.TAC.DF.SAMPLES))

MASTER.TAC.DF.SAMPLES$µM.Copper.Reducing.Equivalents <- MASTER.TAC.DF.SAMPLES$mM.Uric.Acid.Equivalents*2189 # Calculate µM Copper Reducing Equivlents (CRE)

# Summarize by ID to obtain the median (assay completed in duplicates)
MASTER.TAC.MEANS <- MASTER.TAC.DF.SAMPLES %>%  
  dplyr::group_by(Date, ID) %>% # call column to summarize 
  dplyr::summarise(mean.UAE = mean(mM.Uric.Acid.Equivalents), 
                   mean.CRE = mean(µM.Copper.Reducing.Equivalents)) 
MASTER.TAC.MEANS # view the table of mean values (actually the median of duplicates - prints the single value for readings not duplicated)



################################################################################################################################## #
################################################################################################################################## #
###############################################################  MERGE FOR FINAL DATASET  ##########################################
################################################################### OXISELECT ASSSAY     ######################################### #
################################################################################################################################## #
MASTER.TAC.MEANS # calculated total antioxidant capacity
MASTER.TP.MEANS # calculated total antioxidant capacity
Master.Sample.Tracker # master sheet of all data (total homogenized volume, AFDW, volumes for each assay, dolution information, etc.)

# Merge with the master sheet (containing Treatment IDs and Correction factors)
TAC.TP.MERGE <- merge(MASTER.TP.MEANS, MASTER.TAC.MEANS, by = 'ID')
MASTER.PHYS.ASSAY <- merge(TAC.TP.MERGE, Master.Sample.Tracker, by = 'ID')

# correct for homogenate volume as "TOTAL"
# TOTAL.total.protein mg/mL - (total homogenate vol / vol aliquot for protein assay) * (protein measured in aliquot*dilution factor)
sapply(MASTER.PHYS.ASSAY, class) # check the class - occasionally variables that are numeric are characters
MASTER.PHYS.ASSAY$µl.Total.Homog.Volume <- as.numeric(MASTER.PHYS.ASSAY$µl.Total.Homog.Volume)
MASTER.PHYS.ASSAY$µl.Total.protein  <- as.numeric(MASTER.PHYS.ASSAY$µl.Total.protein)  
MASTER.PHYS.ASSAY$TOTAL.total.protein.mg.mL <- (MASTER.PHYS.ASSAY$µl.Total.Homog.Volume/MASTER.PHYS.ASSAY$µl.Total.protein)*(MASTER.PHYS.ASSAY$mean.mg.mL_BCA *MASTER.PHYS.ASSAY$Dil.factor.Total.Protein)

# TOTAL.mean.CRE - (total homogenate vol / vol aliquot for protein assay) * (TAC measured in aliquot*dilution factor)
MASTER.PHYS.ASSAY$µl.TAC <- as.numeric(MASTER.PHYS.ASSAY$µl.TAC)
MASTER.PHYS.ASSAY$mean.CRE <- as.numeric(MASTER.PHYS.ASSAY$mean.CRE)
MASTER.PHYS.ASSAY$Dil.factor.TAC <- as.numeric(MASTER.PHYS.ASSAY$Dil.factor.TAC)
MASTER.PHYS.ASSAY$TOTAL.mean.CRE <-  (MASTER.PHYS.ASSAY$µl.Total.Homog.Volume/MASTER.PHYS.ASSAY$µl.TAC)*(MASTER.PHYS.ASSAY$mean.CRE*MASTER.PHYS.ASSAY$Dil.factor.TAC)

# Total protein corrected for AFDW
MASTER.PHYS.ASSAY$mgProtein.mgAFDW <- (MASTER.PHYS.ASSAY$TOTAL.total.protein.mg.mL / ((MASTER.PHYS.ASSAY$TOTAL_AFDW)*1000))
# Total antioxidant capacity corrected for Total protein per gAFDW
MASTER.PHYS.ASSAY$mMCuRedEqu.mgPro.mgAFDW <-  (MASTER.PHYS.ASSAY$TOTAL.mean.CRE / MASTER.PHYS.ASSAY$mgProtein.mgAFDW)/1000


PHYS.DATA <- MASTER.PHYS.ASSAY %>% dplyr::select(Date.fixed, ID, Treatment, Exprmt.Day, Tank.ID, TOTAL_AFDW, mgProtein.mgAFDW, mMCuRedEqu.mgPro.mgAFDW )
PHYS.DATA <- na.omit(PHYS.DATA)

count.days <- PHYS.DATA %>%  dplyr::group_by(Exprmt.Day) %>% 
  dplyr::summarise(count_replicated.by.day =n())

# write out tables
write.table(PHYS.DATA, file="C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/Output/Phys.Assay.Table.csv", sep=",", row.names = FALSE) #save data to output file

