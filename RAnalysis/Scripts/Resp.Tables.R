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

# choose the data during the 21-day experiment to work with in the script
Final.resp.table <- read.csv(file="Data/SDR_data/Final.resp.rates.csv", header=T) #read table
Final.resp.table.EXP <- Final.resp.table %>% dplyr::filter(row.num > 22) # table for all data on and after 20190723 

Short_resp_table <- Final.resp.table.EXP[,c("Date","Tank.ID", "Treatment.ID.initial", "Treatment.ID.SUBSQ", "resp.COUNT.µg.L.hr.indiv", "resp.MEAN.µg.L.hr.mm")]
Short_resp_table$Treatment.history <- substr(Short_resp_table$Treatment.ID.initial, 1,2)
colnames(Short_resp_table) <- c("Date","Tank.ID", "Treatment.EXP_1", "Treatment.EXP_2", "resp.COUNT.µg.L.hr.indiv", "resp.MEAN.µg.L.hr.mm","Treatment.history")


Resp.pre <- Short_resp_table %>% dplyr::filter(Date %in% 20190723) 
Resp.pre$Treatment.EXP_1 <- "NA" # first exposure NA 
Resp.pre$Treatment.EXP_2 <- "NA" # second exposure NA
Resp.D.1.14 <- Short_resp_table %>% dplyr::filter(Date %in% 20190725:20190807)  # filter data with six treatments
Resp.D.1.14$Treatment.EXP_1 <- substr(Resp.D.1.14$Treatment.EXP_1, 3,3) # first exposure NA
Resp.D.1.14$Treatment.EXP_2 <- "NA" # second exposure NA
Resp.D.15.21 <- Short_resp_table %>% dplyr::filter(Date > 20190807)  # filter data with twelve treatments
Resp.D.15.21$Treatment.history <- substr(Resp.D.15.21$Treatment.EXP_2, 1,2) # assign treatment history as first two characters of the 4 digit treatment in exp_2
Resp.D.15.21$Treatment.EXP_1 <- substr(Resp.D.15.21$Treatment.EXP_2, 3,3) # assign treatment history as first three characters of the 4 digit treatment in exp_2
Resp.D.15.21$Treatment.EXP_2 <- substr(Resp.D.15.21$Treatment.EXP_2, 4,4) # assign treatment history as first three characters of the 4 digit treatment in exp_2

TableFINAL <- rbind(Resp.pre, Resp.D.1.14, Resp.D.15.21)

# write output table
write.table(TableFINAL,"Data/SDR_data/Final_table_for_resp_analysis.csv",sep=",", row.names=FALSE)  # write out to the path names outputNAME




# RELATIVE RESPIRATION RATE TO MEAN AMBIENT RATE FOR EACH TREATMENT --------------------------------------




# resp non controls contains all resp data ny animals NOT in ambient (in MOD or SEV pCO2 conditions)
Resp.non.controls.1 <- Final.resp.table.EXP %>% 
  dplyr::filter(Date < 20190808)  %>% 
  filter(Treatment.ID.initial  %in% c("EH", "AHS", "AHM", "EHS", "EHM"))
Resp.non.controls.2 <- Final.resp.table.EXP %>% 
  dplyr::filter(Date > 20190807)  %>% 
  filter(Treatment.ID.SUBSQ   %in% c("AHAM", "AHMM", "AHSM", "EHAM", "EHMM", "EHSM"))

Resp.non.controls.all <- rbind(Resp.non.controls.1, Resp.non.controls.2) # compile to one data table

Resp.non.controls.AMBIENT.HIST <- Resp.non.controls.all %>% 
  filter(Treatment.ID.SUBSQ %in% c("AHAM", "AHMM", "AHSM")) # only animals from ambient history
Resp.non.controls.ELEVATED.HIST <- Resp.non.controls.all %>% 
  filter(Treatment.ID.SUBSQ %in% c("EHAM", "EHMM", "EHSM")) # only animals from elevated history 

# data table of means for all resp rates by animals in AMBIENT COND's
ambient.means.I <- Final.resp.table.EXP %>% 
  dplyr::filter(Date < 20190808)  %>% 
  dplyr::filter(Treatment.ID.initial %in% c("AH", "AHA", "EHA")) %>% 
  dplyr::group_by(Treatment.ID.initial, Date) %>% # call column to summarize 
  dplyr::summarise(mean_resp.MEAN.µg.L.hr.mm = mean(resp.MEAN.µg.L.hr.mm), # call resp corrected for mean shell length
                   mean_resp.COUNT.µg.L.hr.indiv = mean(resp.COUNT.µg.L.hr.indiv)) # call resp corrected for num. indivs.

ambient.means.II <- Final.resp.table.EXP %>% 
  dplyr::filter(Date > 20190807)  %>% 
  dplyr::filter(Treatment.ID.SUBSQ %in% c("AHAA", "AHMA", "AHSA", "EHAA", "EHMA", "EHSA"))  %>% 
  dplyr::group_by(Treatment.ID.SUBSQ, Date) %>% # call column to summarize 
  dplyr::summarise(mean_resp.MEAN.µg.L.hr.mm = mean(resp.MEAN.µg.L.hr.mm), # call resp corrected for mean shell length
                   mean_resp.COUNT.µg.L.hr.indiv = mean(resp.COUNT.µg.L.hr.indiv)) # call resp corrected for num. indivs.

ambient.all <- rbind(ambient.means.I, ambient.means.II) # compile to one dataframe
# merge the initial and subsequent treatments to a single column in ambient.all
# this is not an elegant way to do it but it ggets the job done...
table <- data.frame(matrix(nrow = 31, ncol = 1)) # create new dataframe with equal rows as ambient all
colnames(table)<-c('Treatment') # name the single column of new dataframe as Treatment
table[1:13,1] <- as.character(ambient.all$Treatment.ID.initial[1:13]) # assign the initial treat
table[14:31,1] <- as.character(ambient.all$Treatment.ID.SUBSQ[14:31]) # assign the subseq treat - this is the target column in ambient all
ambient.all$Treatment <- table$Treatment # add this column to the table - review ambient.all 
ambient.all.2 <- ambient.all %>% 
  dplyr::ungroup() %>%
  dplyr::select(-c('Treatment.ID.initial','Treatment.ID.SUBSQ')) # select all columns for ambient.all.2 BUT the initial and subs treat (now compiled in Treatment)

ambient.means.DAYS.1.14 <- ambient.all.2 %>% # ambient menas for days 1 to 14
  dplyr::filter(Date < 20190808)
ambient.means.DAYS.15.21 <- ambient.all.2 %>% # ambient means for days 15 - 21
  dplyr::filter(Date > 20190807)

Resp.non.controls.all.2 <- Resp.non.controls.all %>% 
  select(Date, Treatment.ID.initial, Treatment.ID.SUBSQ, resp.MEAN.µg.L.hr.mm, resp.COUNT.µg.L.hr.indiv) # select only necessary columns

Resp.non.controls.DAYS.1.14 <- Resp.non.controls.all.2 %>% # filter the first 14 days under 6 treatments - this will simplify merging by treatment in next step
  dplyr::filter(Date < 20190808) 
Resp.non.controls.DAYS.1.14$Treatment <- Resp.non.controls.DAYS.1.14$Treatment.ID.initial # "Treatment" to merge with the ambient means from days 1-14

Resp.non.controls.DAYS.15.21 <- Resp.non.controls.all.2 %>%  # filter the last 7 days (15 - 21) under 12 treatments - this will simplify merging by treatment in next step 
  dplyr::filter(Date > 20190807) 
Resp.non.controls.DAYS.15.21$Treatment <- Resp.non.controls.DAYS.15.21$Treatment.ID.SUBSQ # "Treatment" to merge with the ambient means from days 15-21

# data merge and for loop prep Days 1 - 14 data
Resp.non.controls.DAYS.1.14$Date_Treat <- paste(substr(Resp.non.controls.DAYS.1.14$Date, 1,8), substr(Resp.non.controls.DAYS.1.14$Treatment, 1,2)) # create common column based on date and treatment to prep for for loop
ambient.means.DAYS.1.14$Date_Treat  <- paste(substr(ambient.means.DAYS.1.14$Date, 1,8), substr(ambient.means.DAYS.1.14$Treatment, 1,2)) # create common column based on date and treatment to prep for for loop
# data merge and for loop prep Days 15 - 21 data
Resp.non.controls.DAYS.15.21$Date_Treat <- paste(substr(Resp.non.controls.DAYS.15.21$Date, 1,8), substr(Resp.non.controls.DAYS.15.21$Treatment, 1,3)) # create common column based on date and treatment to prep for for loop
ambient.means.DAYS.15.21$Date_Treat  <- paste(substr(ambient.means.DAYS.15.21$Date, 1,8), substr(ambient.means.DAYS.15.21$Treatment, 1,3)) # create common column based on date and treatment to prep for for loop

# merge the datasets
DAYS.1.14.DATA <- merge(Resp.non.controls.DAYS.1.14, ambient.means.DAYS.1.14, by = c('Date_Treat')) # day 1-14 data ready for the for loop
DAYS.15.21.DATA <- merge(Resp.non.controls.DAYS.15.21, ambient.means.DAYS.15.21, by = c('Date_Treat')) # day 15 - 21 data ready for the for loop

######################################## #
######################################## #
########## Days 1 - 14 for loop ######## #
######################################## #
######################################## #

df_Days.1.14 <- data.frame() # start dataframe 
relative.resp.Days.1.14 <- data.frame(matrix(nrow = 1, ncol = 4)) # create dataframe to save cumunalitively during for loop
colnames(relative.resp.Days.1.14)<-c('Date', 'Treatment', 'rel.resp.COUNT', 'rel.resp.MEAN.LENGTH') # names for comuns in the for loop

as.data.frame(colnames(DAYS.1.14.DATA)) # view column order to make for loop

for(i in 1:nrow(DAYS.1.14.DATA)) {
  relative.resp.Days.1.14$Date <- DAYS.1.14.DATA[i,2] # call Date
  relative.resp.Days.1.14$Treatment <- DAYS.1.14.DATA[i,7] # call Treatment.x
  relative.resp.Days.1.14$rel.resp.COUNT <- DAYS.1.14.DATA[i,6] - DAYS.1.14.DATA[i,10] # (ambient mean for resp corrected for # indivs) - each non ambient resp value of the same treatment history
  relative.resp.Days.1.14$rel.resp.MEAN.LENGTH <- DAYS.1.14.DATA[i,5] - DAYS.1.14.DATA[i,9]  # (ambient mean for resp corrected for # indivs) - each non ambient resp value of the same treatment history
  
  df <- data.frame(relative.resp.Days.1.14) # name dataframe for this singl e row
  df_Days.1.14 <- rbind(df_Days.1.14,df) #bind to a cumulative list dataframe
  print(df_Days.1.14) # print to monitor progress
}
df_Days.1.14 # view data

######################################## #
######################################## #
########## Days 15 - 21 for loop ######## #
######################################## #
######################################## #

df_Days.15.21 <- data.frame() # start dataframe 
relative.resp.Days.15.21 <- data.frame(matrix(nrow = 1, ncol = 4)) # create dataframe to save cumunalitively during for loop
colnames(relative.resp.Days.15.21)<-c('Date', 'Treatment', 'rel.resp.COUNT', 'rel.resp.MEAN.LENGTH') # names for comuns in the for loop

as.data.frame(colnames(DAYS.15.21.DATA)) # view column order to make for loop

for(i in 1:nrow(DAYS.15.21.DATA)) {
  relative.resp.Days.15.21$Date <- DAYS.15.21.DATA[i,2] # call Date
  relative.resp.Days.15.21$Treatment <- DAYS.15.21.DATA[i,7] # call Treatment.x
  relative.resp.Days.15.21$rel.resp.COUNT <- DAYS.15.21.DATA[i,6] - DAYS.15.21.DATA[i,10] # (ambient mean for resp corrected for # indivs) - each non ambient resp value of the same treatment history
  relative.resp.Days.15.21$rel.resp.MEAN.LENGTH <-  DAYS.15.21.DATA[i,5] - DAYS.15.21.DATA[i,9]# (ambient mean for resp corrected for # indivs) - each non ambient resp value of the same treatment history
  
  df <- data.frame(relative.resp.Days.15.21) # name dataframe for this singl e row
  df_Days.15.21 <- rbind(df_Days.15.21,df) #bind to a cumulative list dataframe
  print(df_Days.15.21) # print to monitor progress
}
df_Days.15.21 # view data

SMR_relative <- rbind(df_Days.1.14, df_Days.15.21)

write.table(SMR_relative,"Data/SDR_data/Relative_resp_rates.csv",sep=",", row.names=FALSE)  # write out to the path names outputNAME
