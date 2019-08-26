#Title: Respiration Calculations
#Author: Sam Gurr 
#Edited by: Sam Gurr
#Date Last Modified: 20190705
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
library(devtools)
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') install_github('colin-olito/LoLinR') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("lmtest" %in% rownames(installed.packages()) == 'FALSE') install.packages('lmtest') 
if ("car" %in% rownames(installed.packages()) == 'FALSE') install.packages('car') 

#Read in required libraries
##### Include Versions of libraries
#install_github('colin-olito/LoLinR')
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')
library('car')
library('lmtest')
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')
library('car')
library('lmtest')
library(ggplot2)        # Version 2.2.1, Packaged: 2016-12-30, Depends: R (>= 3.1)Imports: digest, grid, gtable (>= 0.1.1), MASS, plyr (>= 1.7.1),reshape2, scales (>= 0.4.1), stats, tibble, lazyeval
library(ggpubr)         # Version: 0.1.8 Date: 2018-08-30, Depends: R (>= 3.1.0), ggplot2, magrittrImports: ggrepel, grid, ggsci, stats, utils, tidyr, purrr, dplyr(>=0.7.1), cowplot, ggsignif, scales, gridExtra, glue, polynom
library(Rmisc)          # Version: 1.5 Packaged: 2013-10-21, Depends: lattice, plyr
library(plotrix)        # Version: 3.7-4, Date/Publication: 2018-10-03

# Set Working Directory:
# setwd("~/MyProjects/Geoduck_Conditioning/RAnalysis/") #set working
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/")
#Load Size Info
Size.Info <- read.csv(file="Data/Shell_Length/Shell_length_data.csv", header=T) #read Size.info data
Tank.info <- read.csv(file="Data/Tank.ID.reference.csv", header=T) #read Size.info data
Tank.info.2 <- read.csv(file="Data/Tank.ID.reference.subsequent.csv", header=T) #read Size.info data
# call the cumulative resp table of Lolin raw outputs
cumulative_resp_table <- read.csv(file="Data/SDR_data/Cumulative_resp_alpha0.4.csv", header=T) #read Size.info data
cumulative_resp_table <- na.omit(cumulative_resp_table)
# call the Size info of shell length and number of individuals per trial ( Size info called above)
x <- merge(cumulative_resp_table, Size.Info, by=c("Date","SDR_position", "RUN"))

################################################################## #
# SECTION I - create new columns for shell size data  ############ #
################################################################## #

#------------- Add columns to new table (x.2) for MEAN shell length, SUM shell length, & COUNT ----------

# TOTAL SHELL LENGTH, MEAN SHELL LENGTH
SUM.MEAN.LENGTH <- x %>% 
  dplyr::select(Date, SDR_position, RUN, Tank.ID, Length, SDR.vial.volume.ml) %>% # select desired columns
  dplyr::filter(!is.na(Length)) %>% # ommit nas (blanks and empty resp vials)
  dplyr::group_by(SDR_position, Date, RUN) %>% # call column to summarize 
  dplyr::summarise(mean_length = mean(Length), # call mean shell length
            sum_length = sum(Length), # call sum of shell lengths
            count_individuals =n()) %>%  # call for the count of repeats (number of individuals)
  dplyr::arrange(desc(sum_length)) # makes table in descending order 
SUM.MEAN.LENGTH # view table
# merge the table of values to the original - creates new columns for mean, sum, and count (NOTE: blanks are gone)
x.2 <- merge(SUM.MEAN.LENGTH, x, by=c("Date","SDR_position", "RUN")) # NOTE: this repeats for every distinct length value
x.2 #view new table

# modify x.2 -> x.3 by removing unnecessary columns & ommitting repeated values
x.3 <- x.2 %>%  # call table
  dplyr::select(-c(Notes, Length)) # ommit notes and length 
x.3$Date_SDR_Run <- paste(x.3$Date, x.3$SDR_position, x.3$RUN, sep = "_") # merge a unique column
x.3 <- x.3[!duplicated(x.3$Date_SDR_Run), ] # ommit all duplicates based on unique column
x.3 # view table

Tank.info
x.4 <- merge(Tank.info, x.3, by=c("Tank.ID"))
x.4 <- merge(x.3, Tank.info.2, by=c("Tank.ID")) # Creates new column for "Treatment.ID"
# NOTE: Treatment.ID shows EH and AH for history (pedivveliger to jjuvenile rearing) followed the initial and subseqent exposures to A, M, S

# In summary...
# x = all raw respiration data with length, but no size analysis
# use table x to calculate the mean blank values
# x.2 = all data needed, use this data to calculate respiration rates
# after the mean blank value (from table x) is obtained
# x.3 = clean and ready for analysis

########################################### #
# SECTION II - analyze blanks  ############ #
########################################### #

# ------------ (1) ------------#
# loop to average oxygen consumption by treatment and date 

dates.runs <- x %>%  # call table
  filter(!is.na(Sw.Condition)) %>% # ommit empty vials
  distinct(Date, Sw.Condition) # call all unique values for date run and sw condition
dates.runs <- as.data.frame(dates.runs) # make a dataframe

# call dataframe and build table to rbind in for loop
blanks_total <- data.frame() # start dataframe 
blanks.table <- data.frame(matrix(nrow = 1,ncol = 5)) # make a table template
colnames(blanks.table)<-c('Date', 'Sw.Condition', 'BLANK.mean.Lpc', 'BLANK.mean.Leq' , 'BLANK.mean.Lz') # names for comuns in the for loop

# for loop. objective = obtian a mean value for all blanks specific to date, run #, seawater treatment
for(i in 1:nrow(dates.runs)) {
data <- x %>% 
  dplyr::select(Date, Type, Tank.ID, Sw.Condition, Lpc,  Leq, Lz, alpha) %>% 
  dplyr::filter(!is.na(Sw.Condition)) %>% # ommits empty resp vials
  dplyr::filter(Date == dates.runs[i,1])  %>%
  #dplyr::filter(RUN == dates.runs[i,2]) %>%
  dplyr::filter(Sw.Condition == dates.runs[i,2])

blanks <- data %>%
  dplyr::filter(Type == "Blank") %>% 
  dplyr::summarise(mean_Lpc = mean(abs(Lpc)),
                  mean_Leq = mean(abs(Leq)), 
                  mean_Lz = mean(abs(Lz)))

      blanks.table$Date <- dates.runs[i,1] # all files have date in the form of yyyymmdd at the start of each csv name
      #blanks.table$RUN <- dates.runs[i,2] # assign the run to the number in the title for the trials completed that day
      blanks.table$Sw.Condition <- data[1,4]
      blanks.table$BLANK.mean.Lpc <- blanks[1,1]
      blanks.table$BLANK.mean.Leq <- blanks[1,2]
      blanks.table$BLANK.mean.Lz <- blanks[1,3]
     # blanks.table$alpha <- data[1,9] # set at start of script - reresents the proportion of data for final estimate of slopes (Lpc, Leq, Lz)
      
      df <- data.frame(blanks.table) # name dataframe for this singl e row
      blanks_total <- rbind(blanks_total,df) #bind to a cumulative list dataframe
      print(blanks_total) # print to monitor progress
}

blanks_total # view blanks table

# ------------ (2) ------------#
# get the mean and st.dev for the volume of water in the blanks
# this will be used to calculate biovolume 

all.blank.data <- x %>% 
  dplyr::filter(Type == "Blank") # all only blanks

all.blank.data$row.num <- seq.int(nrow(all.blank.data)) # make new column for run numbers

average.blank.volume <- all.blank.data %>% 
  dplyr::filter(row.num > 18) %>%  # call only data after vival volume was measured, before 20190725 just assumed 4 ml
  dplyr::summarise(mean_vol = mean(SDR.vial.volume.ml), # call mean vol
                   sd_vol=sd(SDR.vial.volume.ml)) # call st dev vol
average.blank.volume # view table

ml.blank <- average.blank.volume[1,1] # ml.blank is the call to calculate biovolume

#################################################### #
# SECTION III calculate mean biovolume  ############ #
#################################################### #
as.data.frame(colnames(x.4)) # view column order to make for loop

cumulative.biovolume <- data.frame(matrix(nrow = 1,ncol = 1))
colnames(cumulative.biovolume)<-c('mean.biovolume.ml')

for(i in 1:nrow(x.4)) {
  x.4$mean.biovolume.ml[i] <- ((ml.blank - x.4[i,14])/x.4[i,7])
}
  
##################################################### #
# SECTION IV - calc respiration rate  ############### #
##################################################### #

#----------------------------Respiration rate calculation -------------------------------------------
# data summary
x # used to calc blanks total, not used here
x.4 # target table to calculate resp rates (has only "samples", no blanks or empty vials)
blanks_total # correct for microbial respiration 

CALC.table <- merge(x.4, blanks_total, by=c("Date", "Sw.Condition")) # NOTE: this repeats for every distinct length value

# call dataframe and build table to rbind in for loop
Final.resp.table <- data.frame() # start dataframe 
Cumulative.resp <- data.frame(matrix(nrow = 1,ncol = 20)) # make a table template
colnames(Cumulative.resp)<-c('Date', 'Tank.ID', 'Treatment.ID.initial', 'head.tank.INITIAL' , 'tray.ID.INTIAL',
                             'Treatment.ID.SUBSQ','head.tank.SUBSQ' , 'tray.ID.SUBSQ',
                             'RUN', 'SDR_position', 'Sw.Condition', 'mean_length', 'sum_length', 
                             'count_individuals', 'mean_biovolume', 
                             'resp.RAW.µg.L.hr', 'resp.MEAN.µg.L.hr.mm', 'resp.TOTAL.µg.L.hr.mm', 
                            'resp.COUNT.µg.L.hr.indiv', 'resp.MEAN.µg.L.hr.ml') # names for comuns in the for loop

# NOTE: the raw data is in umol L-1 4ml-1
# "resp.RAW.µgLhr" calculation = ((((((abs(Lpc)) - (BLANK.mean.Lpc))*(4/1000))*(60))*31.998))
# (1) corrects for blank (2) converts to Liters (3) converts to hours (4) converts moles to grams; final unit = µg O2 L-1 h-1
# "resp.MEAN.SHELL" calculation = resp.RAW.µgLhr / mean shell length per vial (individual mean)
# "resp.TOTAL.SHELL" calculation = resp.RAW.µgLhr / sum of shell lengths per vial
# "resp.COUNT" calculation = resp.RAW.µgLhr / number of individuals

as.data.frame(colnames(CALC.table)) # view column order to make for loop

for(i in 1:nrow(CALC.table)) {
  
  vial.vol <- CALC.table[i,14] # volume in the sdr vial measured after respiration run (variable due to size of geoduck(s))
  
  resp.RAW.µgLhr <-
    ((((((abs(CALC.table[i,9])) - (CALC.table[i,27]))*(vial.vol/1000))*(60))*31.998))
  resp.MEAN.SHELL <-
    resp.RAW.µgLhr/CALC.table[i,6]
  resp.TOTAL.SHELL <-
    resp.RAW.µgLhr/CALC.table[i,7]
  resp.COUNT <-
    resp.RAW.µgLhr/CALC.table[i,8]
  resp.MEAN.BIOVOLUME <-
    resp.RAW.µgLhr/CALC.table[i,26]
  
  Cumulative.resp$Date                      <- CALC.table[i,1]
  Cumulative.resp$RUN                       <- CALC.table[i,5]
  Cumulative.resp$SDR_position              <- CALC.table[i,4]
  Cumulative.resp$mean_length               <- CALC.table[i,6]
  Cumulative.resp$sum_length                <- CALC.table[i,7]
  Cumulative.resp$count_individuals         <- CALC.table[i,8]
  Cumulative.resp$mean_biovolume            <-CALC.table[i,26]
  Cumulative.resp$Tank.ID                   <- CALC.table[i,3]
  Cumulative.resp$Sw.Condition              <- CALC.table[i,2]
  Cumulative.resp$Treatment.ID.initial      <- CALC.table[i,20]
  Cumulative.resp$head.tank.INITIAL         <- CALC.table[i,17]
  Cumulative.resp$tray.ID.INTIAL            <- CALC.table[i,18]
  Cumulative.resp$Treatment.ID.SUBSQ        <- CALC.table[i,24]
  Cumulative.resp$head.tank.SUBSQ           <- CALC.table[i,21]
  Cumulative.resp$tray.ID.SUBSQ             <- CALC.table[i,22]
  #Cumulative.resp$Treat.initial             <- CALC.table[i,15]
  #Cumulative.resp$Treat.Secondary           <- CALC.table[i,16]
  Cumulative.resp$resp.RAW.µg.L.hr          <- resp.RAW.µgLhr
  Cumulative.resp$resp.MEAN.µg.L.hr.mm      <- resp.MEAN.SHELL
  Cumulative.resp$resp.TOTAL.µg.L.hr.mm     <- resp.TOTAL.SHELL
  Cumulative.resp$resp.COUNT.µg.L.hr.indiv  <- resp.COUNT
  Cumulative.resp$resp.MEAN.µg.L.hr.ml      <- resp.MEAN.BIOVOLUME
  
  df.resp <- data.frame(Cumulative.resp) # name dataframe for this singl e row
  Final.resp.table <- rbind(Final.resp.table,df.resp) #bind to a cumulative list dataframe
  #print(Final.resp.table) # print to monitor progress
}

Final.resp.table <- Final.resp.table[!is.na(Final.resp.table$resp.RAW.µg.L.hr),] # ommit NAs - there are two in which vial volume was spilled and not measured
Final.resp.table # view table 
#Final.resp.table <- read.csv(file="Data/SDR_data/Final.resp.rates.csv", header=T) #read Size.info data
Final.resp.table$row.num <- seq.int(nrow(Final.resp.table)) # make new column for run numbers

#----------------------OUTPUT - save plots and cumulative tables-----------------------------------------

write.table(blanks_total,"Data/SDR_data/All.blank.resp.rates.csv",sep=",", row.names=FALSE)  # write out to the path names outputNAME
write.table(Final.resp.table,"Data/SDR_data/Final.resp.rates.csv",sep=",", row.names=FALSE)  # write out to the path names outputNAME

