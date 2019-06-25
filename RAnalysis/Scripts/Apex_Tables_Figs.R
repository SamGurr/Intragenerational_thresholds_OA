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


#Required Data files
# ----Conical Chemistry (APEX data)
#20180724_Apex_Data_Output.csv
#20180805_Apex_Data_Output.csv
#20180814_Apex_Data_Output.csv
# ----Heath Tray Chemistry (discrete probe data)
#Flow.rates.csv
#pH_Calibration_Files (tris data)
#Seawater_chemistry_table_Output_All.csv
# ----Respiration data
#Resp.pre.Exposure.csv
#Resp.30d.Exposure.csv
# ----Shell size data
#Length.30d.Exposure.csv

#set working directory--------------------------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/") #set working

### CONICAL Seawater chemistry Data - Analysis, Graphs, Tables (APEX DATA) ####

#CONTINUOUS EXPERIMENTAL APEX DATA 
#Load Apex Data 
APEX_1<-read.csv("Data/Apex_data/Output/20190520_Apex_Data_Output.data.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_2<-read.csv("Data/Apex_data/Output/20190621_Apex_Data_Output.data.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_3<-read.csv("Data/Apex_data/Output/20190625_Apex_Data_Output.data.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_data<- do.call("rbind", list(APEX_1, APEX_2, APEX_3)) # bind all data together
APEX_data$Date.Time <-as.POSIXct(APEX_data$Date.Time, format="%Y-%m-%d %H:%M:%S") #convert date format

#plot raw data
plot(APEX_data$Date.Time,APEX_data$pH_T5) # tail end is after exposure experiment
plot(APEX_data$Date.Time,APEX_data$pH_T6) # tail end is after exposure experiment
plot(APEX_data$Date.Time,APEX_data$pH_T7) # start and tail end show before and after exposure experiment
plot(APEX_data$Date.Time,APEX_data$pH_T3) # tail end is after exposure experiment
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

pH.low_t5$Treatment <- "Low_headtank" #Add treatment Information
colnames(pH.low_t5) <- c("datehour", "days", "mean", "se", "Treatment") #rename columns to generic format
pH.low_t6$Treatment <- "Low_tray" #Add treatment Information
colnames(pH.low_t6) <- c("datehour", "days", "mean", "se", "Treatment") #rename columns to generic format
pH.amb_t7$Treatment <- "Ambient_tray" #Add treatment Information
colnames(pH.amb_t7) <- c("datehour", "days", "mean", "se", "Treatment") #rename columns to generic format

hourly.pH <- rbind(pH.low_t5, pH.low_t6, pH.amb_t7) #bind treatment data 
hourly.pH <- hourly.pH[!is.na(hourly.pH$se), ] # ommit rows with NA for stand error
hourly.pH <- hourly.pH[!(hourly.pH$se > 0.2),] # ommit rows with high stand error (conical cleaning)
hourly.pH #view data

# temp tables and graphs for all exposures 
TMP.low_t5 <- do.call(data.frame,aggregate(TMP_T5 ~ datehour + days, data = APEX_data, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour
TMP.low_t6 <- do.call(data.frame,aggregate(TMP_T6 ~ datehour + days, data = APEX_data, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour
TMP.amb_t7 <- do.call(data.frame,aggregate(TMP_T7 ~ datehour + days, data = APEX_data, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Hour

TMP.low_t5$Treatment <- "Low_headtank" #Add treatment Information
colnames(TMP.low_t5) <- c("datehour", "days", "mean", "se", "Treatment") #rename columns to generic format
TMP.low_t6$Treatment <- "Low_tray" #Add treatment Information
colnames(TMP.low_t6) <- c("datehour", "days", "mean", "se", "Treatment") #rename columns to generic format
TMP.amb_t7$Treatment <- "Ambient_tray" #Add treatment Information
colnames(TMP.amb_t7) <- c("datehour", "days", "mean", "se", "Treatment") #rename columns to generic format

hourly.temp <- rbind(TMP.low_t5, TMP.low_t6, TMP.amb_t7) #bind treatment data 
hourly.temp <- hourly.temp[!is.na(hourly.pH$se), ] # ommit rows with NA for stand error
hourly.temp <- hourly.temp[!(hourly.pH$se > 0.2),] # ommit rows with high stand error (conical cleaning)
hourly.temp #view data

# Plot daily averages of pH data from 4/20 - 6/25  (continuous APEX data)
hourly.pH$datehour <- as.POSIXct(hourly.pH$datehour, format="%Y-%m-%d %H:%M:%S") #format datehour 
pH.Apex.FIG <- ggplot(hourly.pH, aes(x=datehour, y=mean, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  #geom_line() +
  geom_point(aes(x = datehour, y = mean, group=Treatment, color=Treatment),cex=1) + #Plot points using time as the x axis, light as the Y axis and black dots
  geom_errorbar(aes(x=datehour, ymax=mean+se, ymin=mean-se), 
                position=position_dodge(0.9), data=hourly.pH, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("A) Hourly pH (NBS)") + #Label the graph with the main title
  #scale_x_date(date_minor_breaks = "1 day") +
  #scale_x_date(breaks = APEX.pH.Exp1$datehour[seq(1, length(APEX.pH.Exp1$datehour), by = 24)]) +
  ylim(6.8,8.2) + #Set Y axis limits
  xlab("Time") + #Label the X Axis
  ylab("pH (NBS)") + #Label the Y Axis
  #scale_x_date(date_minor_breaks = "1 day") +
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left

pH.Apex.FIG_2 <- pH.Apex.FIG + scale_color_manual(values=c("#009E73", "#0072B2", "#E69F00")) #colorblindess color theme
pH.Apex.FIG_2 # view figure


# Plot daily averages of temperatyre  data from 4/20 - 6/25  (continuous APEX data)
hourly.temp$datehour <- as.POSIXct(hourly.temp$datehour, format="%Y-%m-%d %H:%M:%S") #format datehour 
TMP.Apex.FIG <- ggplot(hourly.temp, aes(x=datehour, y=mean, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  #geom_line() +
  geom_point(aes(x = datehour, y = mean, group=Treatment, color=Treatment),cex=1) + #Plot points using time as the x axis, light as the Y axis and black dots
  geom_errorbar(aes(x=datehour, ymax=mean+se, ymin=mean-se), 
                position=position_dodge(0.9), data=hourly.pH, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("B) Hourly temperature (C)") + #Label the graph with the main title
  #scale_x_date(date_minor_breaks = "1 day") +
  #scale_x_date(breaks = APEX.pH.Exp1$datehour[seq(1, length(APEX.pH.Exp1$datehour), by = 24)]) +
  ylim(10,20) + #Set Y axis limits
  xlab("Time") + #Label the X Axis
  ylab("temperature (C)") + #Label the Y Axis
  #scale_x_date(date_minor_breaks = "1 day") +
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left

TMP.Apex.FIG_2 <- TMP.Apex.FIG + scale_color_manual(values=c("#009E73", "#0072B2", "#E69F00")) #colorblindess color theme
TMP.Apex.FIG_2 # view figure

Fig.continuous.pH.temp <- grid.arrange(arrangeGrob(pH.Apex.FIG_2, 
                                                   left="pH", ncol=1), 
                                            arrangeGrob(TMP.Apex.FIG_2, 
                                                        left="Temperature", ncol=1), ncol=1)
Fig.continuous.pH.temp # view figure

# output figure
ggsave(file="Output/Fig.continuous.pH.temp.pdf", Fig.continuous.pH.temp, width = 12, height = 8, units = c("in"))
