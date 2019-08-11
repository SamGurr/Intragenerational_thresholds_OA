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
  print(Final.resp.table) # print to monitor progress
}

Final.resp.table # view table 
Final.resp.table$row.num <- seq.int(nrow(Final.resp.table)) # make new column for run numbers

Final.resp.table.EXP <- Final.resp.table %>% # table for all data on and after 20190723 
  dplyr::filter(row.num > 22) 
########################################### #
# SECTION V - plot the data ############### #
########################################### #

#---------------------- (1) Mean, total shell size, and biovolume against raw respiration rate 

Raw.vs.MEAN <- ggplot(Final.resp.table.EXP, aes(x = mean_length, y = resp.RAW.µg.L.hr)) +
  geom_point() +
  theme_classic() +
  #theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  #ylim(4,13) + 
  labs(y=expression("Resp.rate.RAW"~(~µg~O[2]*hr^{-1})), 
       x=expression("MEAN.shell.length"~(mm)))
Raw.vs.MEAN <- Raw.vs.MEAN + stat_smooth(method="lm", se=FALSE)
Raw.vs.MEAN
 
Raw.vs.MEAN.BIOVOLUME <- ggplot(Final.resp.table.EXP, aes(x = mean_biovolume,y = resp.RAW.µg.L.hr)) +
  geom_point() +
  theme_classic() +
  #theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  #ylim(4,13) + 
  labs(y=expression("Resp.rate.RAW"~(~µg~O[2]*hr^{-1})), 
       x=expression("MEAN.biovolume"~(ml)))
Raw.vs.MEAN.BIOVOLUME <- Raw.vs.MEAN.BIOVOLUME + stat_smooth(method="lm", se=FALSE)
Raw.vs.MEAN.BIOVOLUME


Biovol.vs.MEAN <- ggplot(Final.resp.table.EXP, aes(x = mean_length,y = mean_biovolume)) +
  geom_point() +
  theme_classic() +
  #theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  #ylim(4,13) + 
  labs(y=expression("MEAN.biovolume"~(ml)), 
       x=expression("MEAN.shell.length"~(mm)))
Biovol.vs.MEAN <- Biovol.vs.MEAN + stat_smooth(method="lm", se=FALSE)
Biovol.vs.MEAN

# Raw.vs.TOTAL<- Raw.vs.TOTAL + stat_smooth(method="lm", se=FALSE) # Adds line
Resp.vs.Biometrics.plot <- ggarrange(Raw.vs.MEAN,Raw.vs.MEAN.BIOVOLUME, Biovol.vs.MEAN,ncol = 3, nrow = 1) # combine plots 
Resp.vs.Biometrics.plot # view plots


#--------------------- (2) Respiration rate plots

# pre, initial, and ambient exposure
RESP.pre.initial.ambient <- Final.resp.table.EXP %>% dplyr::filter(Date < 20190808) # whole dataset 

EH.initial.exposure <- subset(RESP.pre.initial.ambient, grepl("^E", Treatment.ID.initial)) # subset from elevated history
AH.initial.exposure <- subset(RESP.pre.initial.ambient, grepl("^A", Treatment.ID.initial))  # subset from ambient history

# subsequent exposure
RESP.subsequent <- Final.resp.table.EXP %>% dplyr::filter(Date > 20190807) # whole dataset

EH.subs.exposure <- subset(RESP.subsequent, grepl("^E", Treatment.ID.SUBSQ)) # subset from elevated history
EHA.subs.exposure <- subset(EH.subs.exposure, grepl("^EHA", Treatment.ID.SUBSQ)) # elev. history initial ambient
EHM.subs.exposure <- subset(EH.subs.exposure, grepl("^EHM", Treatment.ID.SUBSQ)) # elev. history initial moderate
EHS.subs.exposure <- subset(EH.subs.exposure, grepl("^EHS", Treatment.ID.SUBSQ)) # elev. history initial severe

AH.subs.exposure  <- subset(RESP.subsequent, grepl("^A", Treatment.ID.SUBSQ)) # subset from ambient history
AHA.subs.exposure <- subset(AH.subs.exposure, grepl("^AHA", Treatment.ID.SUBSQ)) # ambient history initial ambient
AHM.subs.exposure <- subset(AH.subs.exposure, grepl("^AHM", Treatment.ID.SUBSQ)) # ambient history initial moderate
AHS.subs.exposure <- subset(AH.subs.exposure, grepl("^AHS", Treatment.ID.SUBSQ)) # ambient history initial severe


# --------------------  PRE, INITIAL EXPOSURE, AMBIENT RECOVERY  ------------------------#
# RAW DATA 
RESP.plots.RAW <- ggplot(RESP.pre.initial.ambient, aes(x = factor(Date), y = resp.RAW.µg.L.hr, fill = Treatment.ID.initial)) + 
  theme_classic() +
  scale_fill_manual(values=c("blue", "blue", "orange", "red", "orange", "blue", "orange", "red")) +
  geom_boxplot(alpha = 0.5, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), 
               width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted", 
             color = "black", size=0.5) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  #ylim(0,10) + 
  labs(y=expression("Respiration rate (raw data corrected to blank)"~µg~O[2]*hr^{-1}), x=expression("Date"))

RESP.plots.RAW # view plot 
#RESP.plots.MEAN.shell <- RESP.plots.MEAN.shell + stat_smooth(method="lm", se=FALSE)

# --------------------  (B) ------------------------#
# corrected for MEAN SHELL LENGTH

RESP.plots.MEAN.shell <- ggplot(RESP.pre.initial.ambient, aes(x = factor(Date), y = resp.MEAN.µg.L.hr.mm, fill = Treatment.ID.initial)) + 
  theme_classic() +
  scale_fill_manual(values=c("blue", "blue", "orange", "red", "orange", "blue", "orange", "red")) +
                  #  labels=c("AH","EH", "Severe")) +
  geom_boxplot(alpha = 0.5, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), 
               width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted", 
             color = "black", size=0.5) +
  #scale_x_discrete(labels = c(8,12,16,20,24,32)) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  #ylim(0,10) + 
  labs(y=expression("Respiration rate (mean shell length)"~µg~O[2]*hr^{-1}*mm^{-1}), x=expression("Date"))
  # + geom_errorbar(aes(ymin=resp.MEAN.µg.L.hr.mm-SEM , ymax=resp.MEAN.µg.L.hr.mm+SEM ), width=.1)

RESP.plots.MEAN.shell <- RESP.plots.MEAN.shell + 
  annotate(geom="text", x=1, y=6, 
           label="Pre-exp (trays)", color="black") +
  annotate(geom="text", x=3, y=6, 
           label="Initial treatment", color="black") +
  annotate(geom="text", x=5, y=6, 
           label="Ambient", color="black")
RESP.plots.MEAN.shell # view plot 
#RESP.plots.MEAN.shell <- RESP.plots.MEAN.shell + stat_smooth(method="lm", se=FALSE)

# --------------------  (C) ------------------------#
# corrected for TOTAL SHELL LENGTH

RESP.plots.TOTAL.shell <- ggplot(RESP.pre.initial.ambient, aes(x = factor(Date), y = resp.TOTAL.µg.L.hr.mm, fill = Treatment.ID.initial)) + 
  theme_classic() +
  scale_fill_manual(values=c("blue", "blue", "orange", "red", "orange", "blue", "orange", "red")) +
                  #labels=c("AH","EH", "Severe")) +
  geom_boxplot(alpha = 0.5, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), 
               width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted", 
             color = "black", size=0.5) +
  #scale_x_discrete(labels = c(8,12,16,20,24,32)) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("Respiration rate (total shell length)"~µg~O[2]*hr^{-1}*mm^{-1}), x=expression("Date"))
  #geom_errorbar(aes(ymin=mean_resp-SEM , ymax=mean_resp+SEM ), width=.1)

RESP.plots.TOTAL.shell <- RESP.plots.TOTAL.shell + 
  annotate(geom="text", x=1, y=2.1, 
           label="Pre-exp (trays)", color="black") +
  annotate(geom="text", x=3, y=2.1, 
           label="Initial treatment", color="black") +
  annotate(geom="text", x=5, y=2.1, 
           label="Ambient", color="black")
RESP.plots.TOTAL.shell # view plot 

# --------------------  (D) ------------------------#
# corrected for NUMBER OF INDIVIDUALS

RESP.plots.INDIV <- ggplot(RESP.pre.initial.ambient, aes(x = factor(Date), y = resp.COUNT.µg.L.hr.indiv, fill = Treatment.ID.initial)) + 
  theme_classic() +
  scale_fill_manual(values=c("blue", "blue", "orange", "red", "orange", "blue", "orange", "red")) +
  #labels=c("AH","EH", "Severe")) +
  geom_boxplot(alpha = 0.5, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), 
               width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted", 
             color = "black", size=0.5) +
  #scale_x_discrete(labels = c(8,12,16,20,24,32)) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("Respiration rate (per individual)"~µg~O[2]*hr^{-1}*indiv^{-1}), x=expression("Date"))
#geom_errorbar(aes(ymin=mean_resp-SEM , ymax=mean_resp+SEM ), width=.1)

RESP.plots.INDIV <- RESP.plots.INDIV + 
  annotate(geom="text", x=1, y=16, 
           label="Pre-exp (trays)", color="black") +
  annotate(geom="text", x=3, y=16, 
           label="Initial treatment", color="black") +
  annotate(geom="text", x=5, y=16, 
           label="Ambient", color="black")
RESP.plots.INDIV # view plot 

# --------------------  (E) ------------------------#
# corrected for MEAN BIOVOLUME

RESP.plots.MEAN.biovolume <- ggplot(RESP.pre.initial.ambient, aes(x = factor(Date), y =  resp.MEAN.µg.L.hr.ml, fill = Treatment.ID.initial)) + 
  theme_classic() +
  scale_fill_manual(values=c("blue", "blue", "orange", "red", "orange", "blue", "orange", "red")) +
  #labels=c("AH","EH", "Severe")) +
  geom_boxplot(alpha = 0.5, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), 
               width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted", 
             color = "black", size=0.5) +
  #scale_x_discrete(labels = c(8,12,16,20,24,32)) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("Respiration rate (per mean biovolume)"~µg~O[2]*hr^{-1}*mL^{-1}), x=expression("Date"))
#geom_errorbar(aes(ymin=mean_resp-SEM , ymax=mean_resp+SEM ), width=.1)

RESP.plots.MEAN.biovolume <- RESP.plots.MEAN.biovolume + annotate(geom="text", x=1, y=500, 
                                     label="Pre-exp (trays)", color="black") +
                            annotate(geom="text", x=3, y=500, 
                                     label="Initial treatment", color="black") +
                            annotate(geom="text", x=5, y=500, 
                                     label="Ambient", color="black")
RESP.plots.MEAN.biovolume # view plot

Resp.plots <- ggarrange(RESP.plots.MEAN.shell,
                        RESP.plots.INDIV, RESP.plots.MEAN.biovolume, ncol= 2, nrow = 2)
Resp.plots # view plots

# ------------------------------ TOTAL PLOTS WITH MEAN SHELL LENGTH -------------------------------------------------#

#######################################
# ELEVATED HISTORY
#######################################

EH.pre.init.amb.PLOT <- ggplot(EH.initial.exposure, aes(x = factor(Date), y = resp.MEAN.µg.L.hr.mm, fill = Treatment.ID.initial)) + 
  theme_classic() + scale_fill_manual(values=c("orange", "blue", "orange", "red", "orange", "blue", "orange", "red")) +
  geom_boxplot(alpha = 0.5, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted",  color = "black", size=0.5) +
  ylim(0,7) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("Respiration rate (mean shell length)"~µg~O[2]*hr^{-1}*mm^{-1}), x=expression("Date"))

EH.pre.init.amb.PLOT <- EH.pre.init.amb.PLOT + ggtitle("ELEVATED HISTORY") +
  annotate(geom="text", x=1, y=6,  label="Pre-exp (trays)", color="black") +
  annotate(geom="text", x=3, y=6,   label="Initial treatment", color="black") +
  annotate(geom="text", x=6, y=6,  label="Ambient", color="black")
EH.pre.init.amb.PLOT # view plot 

EHA.subs.PLOT <- ggplot(EHA.subs.exposure, aes(x = factor(Date), y = resp.MEAN.µg.L.hr.mm, fill = Treatment.ID.SUBSQ)) + 
  theme_classic() + scale_fill_manual(values=c("blue1", "blue4")) +
  geom_boxplot(alpha = 0.75, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted",  color = "black", size=0.5) +
  ylim(0,7) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("Respiration rate (mean shell length)"~µg~O[2]*hr^{-1}*mm^{-1}), x=expression("Date"))

EHS.subs.PLOT <- ggplot(EHS.subs.exposure, aes(x = factor(Date), y = resp.MEAN.µg.L.hr.mm, fill = Treatment.ID.SUBSQ)) + 
  theme_classic() + scale_fill_manual(values=c("red1", "red4")) +
  geom_boxplot(alpha = 0.75, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted",  color = "black", size=0.5) +
  ylim(0,7) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("Respiration rate (mean shell length)"~µg~O[2]*hr^{-1}*mm^{-1}), x=expression("Date"))

EHM.subs.PLOT <- ggplot(EHM.subs.exposure, aes(x = factor(Date), y = resp.MEAN.µg.L.hr.mm, fill = Treatment.ID.SUBSQ)) + 
  theme_classic() + scale_fill_manual(values=c("orange1", "orange4")) +
  geom_boxplot(alpha = 0.75, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  ylim(0,7) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted",  color = "black", size=0.5) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("Respiration rate (mean shell length)"~µg~O[2]*hr^{-1}*mm^{-1}), x=expression("Date"))

EH.sub.plots <- grid.arrange(EHA.subs.PLOT, EHM.subs.PLOT,EHS.subs.PLOT, nrow = 1)
Elevated.plots.ALL <- grid.arrange(EH.pre.init.amb.PLOT, EH.sub.plots, nrow = 2)

#######################################
# AMBIENT HISTORY
#######################################
AH.pre.init.amb.PLOT <- ggplot(AH.initial.exposure, aes(x = factor(Date), y = resp.MEAN.µg.L.hr.mm, fill = Treatment.ID.initial)) + 
  theme_classic() + scale_fill_manual(values=c("blue", "blue", "orange", "red", "orange", "blue", "orange", "red")) +
  geom_boxplot(alpha = 0.5, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted",  color = "black", size=0.5) +
  ylim(0,7) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("Respiration rate (mean shell length)"~µg~O[2]*hr^{-1}*mm^{-1}), x=expression("Date"))

AH.pre.init.amb.PLOT <- AH.pre.init.amb.PLOT + ggtitle("AMBIENT HISTORY") +
  annotate(geom="text", x=1, y=6,  label="Pre-exp (trays)", color="black") +
  annotate(geom="text", x=3, y=6,   label="Initial treatment", color="black") +
  annotate(geom="text", x=6, y=6,  label="Ambient", color="black")
AH.pre.init.amb.PLOT # view plot 

AHA.subs.PLOT <- ggplot(AHA.subs.exposure, aes(x = factor(Date), y = resp.MEAN.µg.L.hr.mm, fill = Treatment.ID.SUBSQ)) + 
  theme_classic() + scale_fill_manual(values=c("blue1", "blue4")) +
  geom_boxplot(alpha = 0.75, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted",  color = "black", size=0.5) +
  ylim(0,7) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("Respiration rate (mean shell length)"~µg~O[2]*hr^{-1}*mm^{-1}), x=expression("Date"))

AHS.subs.PLOT <- ggplot(AHS.subs.exposure, aes(x = factor(Date), y = resp.MEAN.µg.L.hr.mm, fill = Treatment.ID.SUBSQ)) + 
  theme_classic() + scale_fill_manual(values=c("red1", "red4")) +
  geom_boxplot(alpha = 0.75, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted",  color = "black", size=0.5) +
  ylim(0,7) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("Respiration rate (mean shell length)"~µg~O[2]*hr^{-1}*mm^{-1}), x=expression("Date"))

AHM.subs.PLOT <- ggplot(AHM.subs.exposure, aes(x = factor(Date), y = resp.MEAN.µg.L.hr.mm, fill = Treatment.ID.SUBSQ)) + 
  theme_classic() + scale_fill_manual(values=c("orange1", "orange4")) +
  geom_boxplot(alpha = 0.75, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  ylim(0,7) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted",  color = "black", size=0.5) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("Respiration rate (mean shell length)"~µg~O[2]*hr^{-1}*mm^{-1}), x=expression("Date"))

AH.sub.plots <- grid.arrange(AHA.subs.PLOT, AHM.subs.PLOT,AHS.subs.PLOT, nrow = 1)
Ambient.plots.ALL <- grid.arrange(AH.pre.init.amb.PLOT, AH.sub.plots, nrow = 2)

#----------------------OUTPUT - save plots and cumulative tables-----------------------------------------

write.table(blanks_total,"Data/SDR_data/All.blank.resp.rates.csv",sep=",", row.names=FALSE)  # write out to the path names outputNAME
write.table(Final.resp.table,"Data/SDR_data/Final.resp.rates.csv",sep=",", row.names=FALSE)  # write out to the path names outputNAME
ggsave(file="Output/Resp.vs.Size.Biovolume.pdf", Resp.vs.Biometrics.plot, width = 12, height = 8, units = c("in")) # respiration rate plots
ggsave(file="Output/Resp.plots.pdf", Resp.plots, width = 12, height = 8, units = c("in")) # respiration rate plots
ggsave(file="Output/Ambient.plots.ALL.pdf", Ambient.plots.ALL, width = 12, height = 8, units = c("in")) # respiration rate plots
ggsave(file="Output/Elevated.plots.ALLs.pdf", Elevated.plots.ALL, width = 12, height = 8, units = c("in")) # respiration rate plots

############################################################################################ #
################################ STATISTICAL ANALYSIS  ##################################### #
############################################################################################ #

Exp.resp <- Final.resp.table.EXP %>% 
  dplyr::filter(Date == c("20190725", "20190728", "20190731", "20190801")) # filter out day one to test anova effects

#---------------------------- raw data ----------------------------#
Exp.resp$Date <- as.character(Exp.resp$Date)
mod.raw <- aov(resp.MEAN.µg.L.hr.mm ~ Treatment.ID*Date, data = Exp.resp)
shapiro.test(residuals(mod.raw)) # shapiro test (model residuals) normal residuals p-value = 0.9886

par(mfrow=c(1,3)) # hist qq residual diagnostic; set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(residuals(mod.raw)) #plot histogram of residuals
boxplot(residuals(mod.raw)) #plot boxplot of residuals
plot(fitted(day1.mod.raw),residuals(mod.raw))
qqnorm(residuals(mod.raw)) # qqplot
plot(lm(mod.raw)) # plot and test model for heteroscedasticity
leveneTest(mod.raw) # Levene's test for homogeneity p = 0.2795
bptest(lm(mod.raw))  # Breusch-Pagan test p-value = 0.2508
# summary of test
summary(mod.raw) # significant effect of treatment - 0.0406 *
library(lsmeans)        # Version: 2.27-62, Date/Publication: 2018-05-11, Depends: methods, R (>= 3.2)
day1.mod.raw.POSTHOC <- lsmeans(mod.raw, pairwise ~  Treatment.ID*Date) # pariwise Tukey Post-hoc test between repeated treatments
day1.mod.raw.POSTHOC # view post hoc summary
day1.mod.raw.POSTHOC.05 <- cld(day1.mod.raw.POSTHOC, alpha=.05, Letters=letters) #list pairwise tests and letter display p < 0.05
day1.mod.raw.POSTHOC.05

#---------------------------- corrected for mean biovolume ----------------------------#
day1.mod.biovolume <- aov(resp.MEAN.µg.L.hr.ml ~ Treatment.ID, data = DAY1)
shapiro.test(residuals(day1.mod.biovolume)) # shapiro test (model residuals) normal residuals p-value = 0.7549

par(mfrow=c(1,3)) # hist qq residual diagnostic; set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(residuals(day1.mod.biovolume)) #plot histogram of residuals
boxplot(residuals(day1.mod.biovolume)) #plot boxplot of residuals
plot(fitted(day1.mod.biovolume),residuals(day1.mod.biovolume))
qqnorm(residuals(day1.mod.biovolume)) # qqplot
plot(lm(day1.mod.biovolume)) # plot and test model for heteroscedasticity
leveneTest(day1.mod.biovolume) # Levene's test for homogeneity p = 0.9717
bptest(lm(day1.mod.biovolume))  # Breusch-Pagan test p-value = 0.7344
# summary of test
summary(day1.mod.biovolume) # no significant effect, p = 0.459 

#---------------------------- corrected for mean shell length ----------------------------#
day1.mod.mean.length <- aov(resp.MEAN.µg.L.hr.mm  ~ Treatment.ID, data = DAY1)
shapiro.test(residuals(day1.mod.mean.length)) # shapiro test (model residuals) normal residuals p-value = 0.9092

par(mfrow=c(1,3)) # hist qq residual diagnostic; set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(residuals(day1.mod.mean.length)) #plot histogram of residuals
boxplot(residuals(day1.mod.mean.length)) #plot boxplot of residuals
plot(fitted(day1.mod.mean.length),residuals(day1.mod.mean.length))
qqnorm(residuals(day1.mod.mean.length)) # qqplot
plot(lm(day1.mod.mean.length)) # plot and test model for heteroscedasticity
leveneTest(day1.mod.mean.length) # Levene's test for homogeneity p = 0.2852
bptest(lm(day1.mod.mean.length))  # Breusch-Pagan test p-value = 0.2352
# summary of test
summary(day1.mod.mean.length) # significant effect, p = 0.048 *
day1.mod.mean.length.POSTHOC <- lsmeans(day1.mod.mean.length, pairwise ~  Treatment.ID) # pariwise Tukey Post-hoc test between repeated treatments
day1.mod.mean.length.POSTHOC # view post hoc summary
day1.mod.mean.length.POSTHOC.05 <- cld(day1.mod.mean.length.POSTHOC, alpha=.05, Letters=letters) #list pairwise tests and letter display p < 0.05
day1.mod.mean.length.POSTHOC.05
