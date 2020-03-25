#Title: Supplementary Figure - Heathstack Figure
#Project: FFAR
#Author: HM Putnam & Sam Gurr
#Edit by: Sam Gurr
#Date Last Modified: 20200318
#See Readme file for details

rm(list=ls())
# Install packages if not already in your library
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("tidyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('tidyr') 
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
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 

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
library(lubridate)
library(tidyr)

#set working directory--------------------------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/") #set working

#Load Heathstack temp and pH data - hourly; created iin script "Apex_Tables.R"
heathstack.data<-read.csv("Output/Hourly.Apex.Heathstack.Rearing.csv", header=T, sep=",", na.string="NA", as.is=T) 

heathstack.data$datehour2 <- as.POSIXct(heathstack.data$datehour, format="%Y-%m-%d %H:%M:%S")
heathstack.data$datehour2
heathstack.data$date <- substr(heathstack.data$datehour2, 1,10) # sub string for just the date
heathstack.data$time <- substr(heathstack.data$datehour2, 11,18) # sub string for just the time
heathstack.data$day <- day(heathstack.data$date) # use lubridate to assign day
heathstack.data$month <- month(heathstack.data$date) # use lubridate to assign month

Figure_heathstack <- heathstack.data %>% 
  dplyr::filter(Treatment %in% c("AMBIENT_tray", "MODERATE_tray")) %>% 
  tidyr::pivot_longer("mean.temp":"se.ph", # columns pivoting 
      names_to = "measurement", # what is the new column called? names_to = ""
      values_to = "values") %>% 
  dplyr::filter(measurement %in% c("mean.ph", "mean.temp"))  %>% # filter just the mean temp and pH
  ggplot(aes(
      x = datehour2, # call the date formatted
      y = values,  # y as thr av temp
      group = Treatment,  # group by the sensor
      colour=factor(Treatment))) +  # color by treatment
  geom_line() + # make a line
  scale_color_manual(values=c("skyblue3", "tomato1"), 
                     name=NULL,
                     breaks=c("AMBIENT_tray","MODERATE_tray"),
                     labels=c("Ambient","Moderate")) + 
  facet_wrap(~measurement,
      nrow = 2, 
      scales = "free_y", # y axis catered to values discrete
      strip.position = "left",  # create a strip on the left of the graphs
      labeller = as_labeller(c(mean.temp = "Temperature (°C)", mean.ph = "pH (NBS scale)"))) + # label the strip with discrete y axis values
  theme(strip.background = element_blank(),
      strip.placement = "outside") + 
  scale_x_datetime(limits=c(as.POSIXct('2019/04/19 00:00:00'), as.POSIXct('2019/07/23 00:00:00'))) +
  ylab(NULL) + # null the y label beause a strip via 'labeller' created this for discrete labels
  labs(title = "3-Month Conditioning; Initial Exposure",
     x = "Month") # label the title and x axis

Figure_heathstack # view the Figure


#Load 21-day experiment temp and pH data - hourly; created iin script "Apex_Tables.R"
experiment.data<-read.csv("Output/Hourly.Apex.Thresholds.Experiment.csv", header=T, sep=",", na.string="NA", as.is=T) 

experiment.data$datehour2 <- as.POSIXct(experiment.data$datehour, format="%Y-%m-%d %H:%M:%S")
experiment.data$datehour2
experiment.data$date <- substr(experiment.data$datehour2, 1,10) # sub string for just the date
experiment.data$time <- substr(experiment.data$datehour2, 11,18) # sub string for just the time
experiment.data$day <- day(experiment.data$date) # use lubridate to assign day
experiment.data$month <- month(experiment.data$date) # use lubridate to assign month

experiment.data2 <- experiment.data %>% filter(date %in% "2019-07-24":"2019-08-13")

Figure_experiment <- experiment.data %>% 
  dplyr::filter(Treatment %in% c("AMBIENT_conical", "MODERATE_tank", "SEVERE_tank")) %>% 
  tidyr::pivot_longer("mean.temp":"se.ph", # columns pivoting 
                      names_to = "measurement", # what is the new column called? names_to = ""
                      values_to = "values") %>% 
  dplyr::filter(measurement %in% c("mean.ph", "mean.temp"))  %>% # filter just the mean temp and pH
  ggplot(aes(
    x = datehour2, # call the date formatted
    y = values,  # y as thr av temp
    group = Treatment,  # group by the sensor
    colour=factor(Treatment))) +  # color by treatment
  geom_line() + # make a line
  scale_x_datetime(limits=c(as.POSIXct('2019/07/24 00:00:00'), as.POSIXct('2019/08/13 23:59:59'))) +
  scale_color_manual(values=c("skyblue3", "tomato1", "firebrick4"), 
                     name=NULL,
                     breaks=c("AMBIENT_conical", "MODERATE_tank", "SEVERE_tank"),
                     labels=c("Ambient","Moderate", "Severe")) + 
  facet_wrap(~measurement,
             nrow = 2, 
             scales = "free_y", # y axis catered to values discrete
             strip.position = "left",  # create a strip on the left of the graphs
             labeller = as_labeller(c(mean.temp = "Temperature (°C)", mean.ph = "pH (NBS scale)"))) + # label the strip with discrete y axis values
  theme(strip.background = element_blank(),
        strip.placement = "outside") + 
  ylab(NULL) + # null the y label beause a strip via 'labeller' created this for discrete labels
  labs(title = "21-Day Experiment; Secondary and Tertiary Exposures",
       x = "Month") # label the title and x axis

Figure_experiment # view the Figure  
