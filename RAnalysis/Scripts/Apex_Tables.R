#Title: Juvenile Repeat Exposure Experiment 2018
#Project: FFAR
#Author: HM Putnam & Sam Gurr
#Edit by: Sam Gurr
#Date Last Modified: 20190625
#See Readme file for details

rm(list=ls())
# Install packages if not already in your library
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("ggplot2" %in% rownames(installed.packages()) == 'FALSE') install.packages('ggplot2') 
if ("ggpubr" %in% rownames(installed.packages()) == 'FALSE') install_github('ggpubr') 
if ("Rmisc" %in% rownames(installed.packages()) == 'FALSE') install.packages('Rmisc') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("lsmeans" %in% rownames(installed.packages()) == 'FALSE') install.packages('lsmeans') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("reshape" %in% rownames(installed.packages()) == 'FALSE') install.packages('reshape') 
if ("multcompView" %in% rownames(installed.packages()) == 'FALSE') install.packages('multcompView') 
if ("lmtest" %in% rownames(installed.packages()) == 'FALSE') install.packages('lmtest') 
if ("car" %in% rownames(installed.packages()) == 'FALSE') install.packages('car') 

# Load packages and pacage version/date/import/depends info
library(dplyr)          # Version 0.7.6, Packaged: 2018-06-27, Depends: R (>= 3.1.2)Imports: assertthat (>= 0.2.0), bindrcpp (>= 0.2.0.9000), glue (>=1.1.1), magrittr (>= 1.5), methods, pkgconfig (>= 2.0.1), R6(>= 2.2.2), Rcpp (>= 0.12.15), rlang (>= 0.2.0), tibble (>=1.3.1), tidyselect (>= 0.2.3), utils
library(ggplot2)        # Version 2.2.1, Packaged: 2016-12-30, Depends: R (>= 3.1)Imports: digest, grid, gtable (>= 0.1.1), MASS, plyr (>= 1.7.1),reshape2, scales (>= 0.4.1), stats, tibble, lazyeval
library(ggpubr)         # Version: 0.1.8 Date: 2018-08-30, Depends: R (>= 3.1.0), ggplot2, magrittrImports: ggrepel, grid, ggsci, stats, utils, tidyr, purrr, dplyr(>=0.7.1), cowplot, ggsignif, scales, gridExtra, glue, polynom
library(Rmisc)          # Version: 1.5 Packaged: 2013-10-21, Depends: lattice, plyr
library(plotrix)        # Version: 3.7-4, Date/Publication: 2018-10-03
library(lsmeans)        # Version: 2.27-62, Date/Publication: 2018-05-11, Depends: methods, R (>= 3.2)
library(gridExtra)      # Version: 2.3, Date/Publication: 2017-09-09, Imports: gtable, grid, grDevices, graphics, utils
library(reshape)        # Version: 0.8.7, Date/Publication: 2017-08-06, Depends: R (>= 2.6.1) Imports: plyr
library(multcompView)   # Version: 0.1-7, Date/Publication: 2015-07-31, Imports: grid
library(Rmisc)
library(lmtest)
library(car)

#set working directory--------------------------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/") #set working

### CONICAL Seawater chemistry Data - Analysis, Graphs, Tables (APEX DATA) ####

#CONTINUOUS EXPERIMENTAL APEX DATA 
#Load Apex Data 
APEX_1<-read.csv("Data/Apex_data/Output/20190520_Apex_Data_Output.data.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_2<-read.csv("Data/Apex_data/Output/20190621_Apex_Data_Output.data.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_3<-read.csv("Data/Apex_data/Output/20190625_Apex_Data_Output.data.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_4<-read.csv("Data/Apex_data/Output/20190722_Apex_Data_Output.data.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_5<-read.csv("Data/Apex_data/Output/20190805_Apex_Data_Output.data.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_6<-read.csv("Data/Apex_data/Output/20190809_Apex_Data_Output.data.csv", header=T, sep=",", na.string="NA", as.is=T) 

APEX_data<- do.call("rbind", list(APEX_1, APEX_2, APEX_3, APEX_4, APEX_5, APEX_6)) # bind all data together
APEX_data$Date.Time <-as.POSIXct(APEX_data$Date.Time, format="%Y-%m-%d %H:%M:%S") #convert date format
APEX_data$X <- 1:nrow(APEX_data) # make new column for cumulative number of rows

##########################################      HEATH STACK CHEMISTRY    #################################################
# ------------------------------------- pediveliger to juvenile in two treatments ---------------------------------------#
# T5 - pH and temperature of head tank of the moderate low pH (elevated pCO2)
# T6 - pH and temperature of H1_T1; a low pH tray in the heath stack
# T7 - pH of H1_T2 of an ambient tray in the heat stack + temp of the ambient head tank conical T7 (senor did not reach H1_T2)

#plot raw data 
plot(APEX_data$Date.Time,APEX_data$pH_T5) # tail end is after exposure experiment
plot(APEX_data$Date.Time,APEX_data$pH_T6) # tail end is after exposure experiment
plot(APEX_data$Date.Time,APEX_data$pH_T7) # start and tail end show before and after exposure experiment
plot(APEX_data$Date.Time,APEX_data$TMP_T5)
plot(APEX_data$Date.Time,APEX_data$TMP_T6)
plot(APEX_data$Date.Time,APEX_data$TMP_T7)


APEX_data$datehour <- cut(as.POSIXct(APEX_data$Date.Time),
                          format="%d-%m-%Y %H:%M:%S", breaks="hour") # create new column for datehour to aggregate data
head(APEX_data) # check this new column - all 10 minute increment data is called for the hour as "datehour"

date <- format(as.POSIXct(APEX_data$Date.Time) ,format = "%Y-%m-%d %H") # call year month and date
days <- as.Date(date) - as.Date(date[1]) # subtract from start to get number of days
days <- as.numeric(days) # convert to numeric

# pH tables and graphs for all exposures 
pH.low_t5 <- do.call(data.frame,aggregate(pH_T5 ~ datehour + days, data = APEX_data, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour
pH.low_t6 <- do.call(data.frame,aggregate(pH_T6 ~ datehour + days, data = APEX_data, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour
pH.amb_t7 <- do.call(data.frame,aggregate(pH_T7 ~ datehour + days, data = APEX_data, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour

pH.low_t5$Treatment <- "MODERATE_conical" #Add treatment Information
colnames(pH.low_t5) <- c("datehour", "days", "mean.ph", "se.ph", "Treatment") #rename columns to generic format
pH.low_t6$Treatment <- "MODERATE_tray" #Add treatment Information
colnames(pH.low_t6) <- c("datehour", "days", "mean.ph", "se.ph", "Treatment") #rename columns to generic format
pH.amb_t7$Treatment <- "AMBIENT_tray" #Add treatment Information
colnames(pH.amb_t7) <- c("datehour", "days", "mean.ph", "se.ph", "Treatment") #rename columns to generic format

hourly.pH.heathstack <- rbind(pH.low_t5, pH.low_t6, pH.amb_t7) #bind treatment data 
hourly.pH.heathstack <- hourly.pH.heathstack[!is.na(hourly.pH.heathstack$se), ] # ommit rows with NA for stand error
hourly.pH.heathstack <- hourly.pH.heathstack[!(hourly.pH.heathstack$se > 0.09),] # ommit rows with high stand error (conical cleaning)
hourly.pH.heathstack #view data

# temp tables and graphs for all exposures 
TMP.low_t5 <- do.call(data.frame,aggregate(TMP_T5 ~ datehour + days, data = APEX_data, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour
TMP.low_t6 <- do.call(data.frame,aggregate(TMP_T6 ~ datehour + days, data = APEX_data, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour
TMP.amb_t7 <- do.call(data.frame,aggregate(TMP_T7 ~ datehour + days, data = APEX_data, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour

TMP.low_t5$Treatment <- "MODERATE_conical" #Add treatment Information
colnames(TMP.low_t5) <- c("datehour", "days", "mean.temp", "se.temp", "Treatment") #rename columns to generic format
TMP.low_t6$Treatment <- "MODERATE_tray" #Add treatment Information
colnames(TMP.low_t6) <- c("datehour", "days", "mean.temp", "se.temp", "Treatment") #rename columns to generic format
TMP.amb_t7$Treatment <- "AMBIENT_tray" #Add treatment Information
colnames(TMP.amb_t7) <- c("datehour", "days", "mean.temp", "se.temp", "Treatment") #rename columns to generic format

hourly.temp.heathstack <- rbind(TMP.low_t5, TMP.low_t6, TMP.amb_t7) #bind treatment data 
hourly.temp.heathstack <- hourly.temp.heathstack[!is.na(hourly.temp.heathstack$se), ] # ommit rows with NA for stand error
hourly.temp.heathstack <- hourly.temp.heathstack[!(hourly.temp.heathstack$se > 0.2),] # ommit rows with high stand error (conical cleaning)
hourly.temp.heathstack #view data

heathstack.chem.FULL <- merge(hourly.temp.heathstack, hourly.pH.heathstack, by = c('datehour', 'Treatment', 'days'))


##########################################      EXPERIMENT CHEMISTRY    #################################################
# ------------------------------------- pediveliger to juvenile in two treatments ---------------------------------------#
# T0 - AMBIENT head tank - feeds tanks 25 - 36 in water bath trays 0_1 and 0_2
# T1 - SEVERE low pH head tank - feed tanks 13 - 24 in water bath trays 1_1 and 1_2
# T2 - MODERATE low pH head tank - feed tanks 1 - 12 in water bath trays 2_1 and 2_2
# T3 - in tank #16; tray 1_1 - pH and temperature of SEVERE treatment from conical T1
# T4 - in tank #4; tray 2_1 - pH and temperature of MODERATE treatment from conical T1
APEX_data[,1]
APEX_data[c(10878),] #  july 6 at 00:00:00 - x = 870; use this to filter data
APEX_data.experiment <- APEX_data %>% 
  filter(X  > 10878)
head(APEX_data.experiment) # data for this analysis/figure output will start with 7/6/19

#plot raw data 
plot(APEX_data.experiment$Date.Time,APEX_data.experiment$pH_T0) # tail end is after exposure experiment
plot(APEX_data.experiment$Date.Time,APEX_data.experiment$pH_T1) # tail end is after exposure experiment
plot(APEX_data.experiment$Date.Time,APEX_data.experiment$pH_T2) # start and tail end show before and after exposure experiment
plot(APEX_data.experiment$Date.Time,APEX_data.experiment$pH_T3) # tail end is after exposure experiment
plot(APEX_data.experiment$Date.Time,APEX_data.experiment$pH_T4) # tail end is after exposure experiment

plot(APEX_data.experiment$Date.Time,APEX_data.experiment$TMP_T0)
plot(APEX_data.experiment$Date.Time,APEX_data.experiment$TMP_T1)
plot(APEX_data.experiment$Date.Time,APEX_data.experiment$TMP_T2)
plot(APEX_data.experiment$Date.Time,APEX_data.experiment$TMP_T4) 

APEX_data.experiment$datehour <- cut(as.POSIXct(APEX_data.experiment$Date.Time),
                          format="%d-%m-%Y %H:%M:%S", breaks="hour") # create new column for datehour to aggregate data
head(APEX_data.experiment) # check this new column - all 10 minute increment data is called for the hour as "datehour"

date <- format(as.POSIXct(APEX_data.experiment$Date.Time) ,format = "%Y-%m-%d %H") # call year month and date
days <- as.Date(date) - as.Date(date[1]) # subtract from start to get number of days
days <- as.numeric(days) # convert to numeric

# pH tables and graphs for all exposures 
pH.amb.conical_T0 <- do.call(data.frame,aggregate(pH_T0 ~ datehour + days, data = APEX_data.experiment, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour
pH.severe.conical_T1 <- do.call(data.frame,aggregate(pH_T1 ~ datehour + days, data = APEX_data.experiment, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour
pH.moderate.conical_T2 <- do.call(data.frame,aggregate(pH_T2 ~ datehour + days, data = APEX_data.experiment, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour
pH.severe.tank_T3 <- do.call(data.frame,aggregate(pH_T3 ~ datehour + days, data = APEX_data.experiment, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour
pH.moderate.tank_T4 <- do.call(data.frame,aggregate(pH_T4 ~ datehour + days, data = APEX_data.experiment, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour

pH.amb.conical_T0$Treatment <- "AMBIENT_conical" #Add treatment Information
colnames(pH.amb.conical_T0) <- c("datehour", "days", "mean.ph", "se.ph", "Treatment") #rename columns to generic format
pH.severe.conical_T1$Treatment <- "SEVERE_conical" #Add treatment Information
colnames(pH.severe.conical_T1) <- c("datehour", "days", "mean.ph", "se.ph", "Treatment") #rename columns to generic format
pH.moderate.conical_T2$Treatment <- "MODERATE_conical" #Add treatment Information
colnames(pH.moderate.conical_T2) <- c("datehour", "days", "mean.ph", "se.ph", "Treatment") #rename columns to generic format
pH.severe.tank_T3$Treatment <- "SEVERE_tank" #Add treatment Information
colnames(pH.severe.tank_T3) <- c("datehour", "days", "mean.ph", "se.ph", "Treatment") #rename columns to generic format
pH.moderate.tank_T4$Treatment <- "MODERATE_tank" #Add treatment Information
colnames(pH.moderate.tank_T4) <- c("datehour", "days", "mean.ph", "se.ph", "Treatment") #rename columns to generic format

hourly.pH.experiment <- rbind(pH.amb.conical_T0, pH.severe.conical_T1, pH.moderate.conical_T2,
                              pH.severe.tank_T3, pH.moderate.tank_T4) #bind treatment data 
hourly.pH.experiment <- hourly.pH.experiment[!is.na(hourly.pH.experiment$se), ] # ommit rows with NA for stand error
hourly.pH.experiment <- hourly.pH.experiment[!(hourly.pH.experiment$se > 0.09),] # ommit rows with high stand error (conical cleaning)
hourly.pH.experiment #view data

# temp tables and graphs for all exposures 
TMP.amb.conical_T0 <- do.call(data.frame,aggregate(TMP_T0 ~ datehour + days, data = APEX_data.experiment, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour
TMP.severe.conical_T1 <- do.call(data.frame,aggregate(TMP_T1 ~ datehour + days, data = APEX_data.experiment, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour
TMP.moderate.conical_T2 <- do.call(data.frame,aggregate(TMP_T2 ~ datehour + days, data = APEX_data.experiment, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour
TMP.severe.tank_T3 <- do.call(data.frame,aggregate(TMP_T3 ~ datehour + days, data = APEX_data.experiment, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour
TMP.moderate.tank_T4 <- do.call(data.frame,aggregate(TMP_T4 ~ datehour + days, data = APEX_data.experiment, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour

TMP.amb.conical_T0$Treatment <- "AMBIENT_conical" #Add treatment Information
colnames(TMP.amb.conical_T0) <- c("datehour", "days", "mean.temp", "se.temp", "Treatment") #rename columns to generic format
TMP.severe.conical_T1$Treatment <- "SEVERE_conical" #Add treatment Information
colnames(TMP.severe.conical_T1) <- c("datehour", "days", "mean.temp", "se.temp", "Treatment") #rename columns to generic format
TMP.moderate.conical_T2$Treatment <- "MODERATE_conical" #Add treatment Information
colnames(TMP.moderate.conical_T2) <- c("datehour", "days", "mean.temp", "se.temp", "Treatment") #rename columns to generic format
TMP.severe.tank_T3$Treatment <- "SEVERE_tank" #Add treatment Information
colnames(TMP.severe.tank_T3) <- c("datehour", "days", "mean.temp", "se.temp", "Treatment") #rename columns to generic format
TMP.moderate.tank_T4$Treatment <- "MODERATE_tank" #Add treatment Information
colnames(TMP.moderate.tank_T4) <- c("datehour", "days", "mean.temp", "se.temp", "Treatment") #rename columns to generic format

hourly.temp.experiment <- rbind(TMP.amb.conical_T0, TMP.severe.conical_T1, TMP.moderate.conical_T2,
                                TMP.severe.tank_T3, TMP.moderate.tank_T4) #bind treatment data 
hourly.temp.experiment <- hourly.temp.experiment[!is.na(hourly.temp.experiment$se), ] # ommit rows with NA for stand error
hourly.temp.experiment <- hourly.temp.experiment[!(hourly.temp.experiment$se > 0.2),] # ommit rows with high stand error (conical cleaning)
hourly.temp.experiment #view data

thresholds.exp.hrly.APEX.FULL <- merge(hourly.temp.experiment, hourly.pH.experiment, by = c('datehour', 'Treatment', 'days'))

################################3
# write tables
write.table(thresholds.exp.hrly.APEX.FULL, "Output/Hourly.Apex.Thresholds.Experiment.csv", sep=",", row.names = FALSE) #save data
write.table(heathstack.chem.FULL, "Output/Hourly.Apex.Heathstack.Rearing.csv", sep=",", row.names = FALSE) #save data

