#Title: Chem Figs
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

# upload data
heathstack.chem.APEX<-read.csv("Output/Hourly.Apex.Heathstack.Rearing.csv", header=T, sep=",", na.string="NA", as.is=T) 
thresholds.exp.APEX<-read.csv("Output/Hourly.Apex.Thresholds.Experiment.csv", header=T, sep=",", na.string="NA", as.is=T) 


# Plot daily averages of pH data from 4/20 - 6/25  (continuous APEX data)
heathstack.TRAY.CHEM <- heathstack.chem.APEX %>% 
  filter(Treatment %in% c("AMBIENT_tray","MODERATE_tray")) # call only sensors in the trays (ommit the conical sensors)

heathstack.TRAY.CHEM$datehour <- as.POSIXct(heathstack.TRAY.CHEM$datehour, format="%Y-%m-%d %H:%M:%S") #format datehour 
heathstack.TRAY.CHEM[4461,] # 7/24 is the first day of the 21 day experiment, do all days 96 and after should be ommitted for this figure 
heathstack.TRAY.CHEM <- heathstack.TRAY.CHEM %>%  filter(heathstack.TRAY.CHEM$days < 96)

FIG.pH.Apex.heathstack <- ggplot(heathstack.TRAY.CHEM, aes(x=datehour, y=mean.ph , group=Treatment, color=Treatment)) + #Plot average diurnal cycle of temperature data
  #geom_line() +
  geom_point(aes(x = datehour, y = mean.ph, group=Treatment, color=Treatment),cex=1) + #Plot points using time as the x axis, light as the Y axis and black dots
  geom_errorbar(aes(x=datehour, ymax=mean.ph + se.ph, ymin=mean.ph - se.ph), 
                position=position_dodge(0.9), data=heathstack.TRAY.CHEM, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("A) Heath stack hourly pH (NBS)") + #Label the graph with the main title
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
FIG.pH.Apex.heathstack_2 <- FIG.pH.Apex.heathstack + scale_color_manual(values=c("#0072B2", "#E69F00")) #colorblindess color theme
FIG.pH.Apex.heathstack_2 # view figure

# Plot daily averages of temperatyre  data from 4/20 - 6/25  (continuous APEX data)
heathstack.TRAY.CHEM$datehour <- as.POSIXct(heathstack.TRAY.CHEM$datehour, format="%Y-%m-%d %H:%M:%S") #format datehour 
FIG.temp.Apex.heathstack<- ggplot(heathstack.TRAY.CHEM, aes(x=datehour, y=mean.temp , group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  #geom_line() +
  geom_point(aes(x = datehour, y = mean.temp , group=Treatment, color=Treatment),cex=1) + #Plot points using time as the x axis, light as the Y axis and black dots
  geom_errorbar(aes(x=datehour, ymax=mean.temp +se.temp, ymin=mean.temp -se.temp), 
                position=position_dodge(0.9), data=heathstack.TRAY.CHEM, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("B) Heath stack hourly temperature (C)") + #Label the graph with the main title
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
FIG.temp.Apex.heathstack_2 <- FIG.temp.Apex.heathstack + scale_color_manual(values=c("#0072B2", "#E69F00")) #colorblindess color theme
FIG.temp.Apex.heathstack_2 # view figure


# ---------------# DAYS 1 - 7 - Amb v. Noderate vs. Severe exposure #-----------------#
thresholds.exp.APEX$datehour <- as.POSIXct(thresholds.exp.APEX$datehour, format="%Y-%m-%d %H:%M:%S") #format datehour 
DAYS.1.to.7.hourly.pH <-  thresholds.exp.APEX %>% 
  dplyr::filter(days %in% (18:25)) %>%  # pH data
  dplyr::filter(Treatment %in% (c("MODERATE_tank","AMBIENT_conical", "SEVERE_tank")))
DAYS.1.to.7.hourly.temp <-  thresholds.exp.APEX %>% 
  dplyr::filter(days %in% (18:25)) %>% # temperature data
  dplyr::filter(Treatment %in% (c("MODERATE_tank","AMBIENT_conical", "SEVERE_tank")))

# PH.FIGURE
FIG.pH.DAYS.1.to.7 <- ggplot(DAYS.1.to.7.hourly.pH, aes(x=datehour, y=mean.ph, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(aes(x = datehour, y = mean.ph, group=Treatment, color=Treatment),cex=1) + #Plot points using time as the x axis, light as the Y axis and black dots
  geom_errorbar(aes(x=datehour, ymax=mean.ph+se.ph, ymin=mean.ph-se.ph), position=position_dodge(0.9), data=DAYS.1.to.7.hourly.pH, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Days 1-7 hourly pH (NBS)") + #Label the graph with the main title
  ylim(6.5,8) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("pH (NBS)") + #Label the Y Axis
  theme_bw() + #Set the background color
  geom_line() + # add line to plot
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
FIG.pH.DAYS.1.to.7 <- FIG.pH.DAYS.1.to.7 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
FIG.pH.DAYS.1.to.7 # view figure

# temp.FIGURE
FIG.temp.DAYS.1.to.7 <- ggplot(DAYS.1.to.7.hourly.temp, aes(x=datehour, y=mean.temp, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(aes(x = datehour, y = mean.temp, group=Treatment, color=Treatment),cex=1) + #Plot points using time as the x axis, light as the Y axis and black dots
  geom_errorbar(aes(x=datehour, ymax=mean.temp+se.temp, ymin=mean.temp-se.temp), position=position_dodge(0.9), data=DAYS.1.to.7.hourly.temp, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Days 1-7 hourly temperature (C)") + #Label the graph with the main title
  ylim(15,20) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("Temperature (C)") + #Label the Y Axis
  theme_bw() + #Set the background color
  geom_line() + # add line to plot
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
FIG.temp.DAYS.1.to.7 <- FIG.temp.DAYS.1.to.7 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
FIG.temp.DAYS.1.to.7 # view figure

# ---------------# DAYS 8 - 14 - AMBIENT PERIOD #-----------------# 
DAYS.8.to.14.hourly.pH <-  thresholds.exp.APEX %>% 
  dplyr::filter(days %in% (26:32))  %>%  # pH data
  dplyr::filter(Treatment %in% (c("MODERATE_tank","AMBIENT_conical", "SEVERE_tank")))
DAYS.8.to.14.hourly.temp <-  thresholds.exp.APEX %>% 
  dplyr::filter(days %in% (26:32)) %>% # temperature data
  dplyr::filter(Treatment %in% (c("MODERATE_tank","AMBIENT_conical", "SEVERE_tank")))

# PH.FIGURE
FIG.pH.DAYS.8.to.14 <- ggplot(DAYS.8.to.14.hourly.pH, aes(x=datehour, y=mean.ph, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(aes(x = datehour, y = mean.ph, group=Treatment, color=Treatment),cex=1) + #Plot points using time as the x axis, light as the Y axis and black dots
  geom_errorbar(aes(x=datehour, ymax=mean.ph+se.ph, ymin=mean.ph-se.ph), position=position_dodge(0.9), data=DAYS.8.to.14.hourly.pH, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Days 7-14hourly pH (NBS)") + #Label the graph with the main title
  ylim(6.5,8) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("pH (NBS)") + #Label the Y Axis
  theme_bw() + #Set the background color
  geom_line() + # add line to plot
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
FIG.pH.DAYS.8.to.14 <- FIG.pH.DAYS.8.to.14 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
FIG.pH.DAYS.8.to.14 # view figure

# temp.FIGURE
FIG.temp.DAYS.8.to.14 <- ggplot(DAYS.8.to.14.hourly.temp, aes(x=datehour, y=mean.temp, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(aes(x = datehour, y = mean.temp, group=Treatment, color=Treatment),cex=1) + #Plot points using time as the x axis, light as the Y axis and black dots
  geom_errorbar(aes(x=datehour, ymax=mean.temp+se.temp, ymin=mean.temp-se.temp), position=position_dodge(0.9), data=DAYS.8.to.14.hourly.temp, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Days 7-14 hourly temperature (C)") + #Label the graph with the main title
  ylim(15,20) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("Temperature (C)") + #Label the Y Axis
  theme_bw() + #Set the background color
  geom_line() + # add line to plot
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
FIG.temp.DAYS.8.to.14 <- FIG.temp.DAYS.8.to.14 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
FIG.temp.DAYS.8.to.14 # view figure

# ---------------# DAYS  15 - 21 - Amb v. Moderate vs. Severe exposure #-----------------# 
DAYS.15.to.21.hourly.pH <-  thresholds.exp.APEX %>% 
  dplyr::filter(days %in% (33:39))  %>%  # pH data
  dplyr::filter(Treatment %in% (c("MODERATE_tank","AMBIENT_conical")))
DAYS.15.to.21.hourly.temp <-  thresholds.exp.APEX %>% 
  dplyr::filter(days %in% (33:39)) %>% # temperature data
  dplyr::filter(Treatment %in% (c("MODERATE_tank","AMBIENT_conical")))
tail(DAYS.15.to.21.hourly.temp)

# PH.FIGURE
FIG.pH.DAYS.15.to.21 <- ggplot(DAYS.15.to.21.hourly.pH, aes(x=datehour, y=mean.ph, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(aes(x = datehour, y = mean.ph, group=Treatment, color=Treatment),cex=1) + #Plot points using time as the x axis, light as the Y axis and black dots
  geom_errorbar(aes(x=datehour, ymax=mean.ph+se.ph, ymin=mean.ph-se.ph), position=position_dodge(0.9), data=DAYS.15.to.21.hourly.pH, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Days 15-21 hourly pH (NBS)") + #Label the graph with the main title
  ylim(6.5,8) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("pH (NBS)") + #Label the Y Axis
  theme_bw() + #Set the background color
  geom_line() + # add line to plot
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
FIG.pH.DAYS.15.to.21 <- FIG.pH.DAYS.15.to.21 + scale_color_manual(values=c("#56B4E9",  "#E69F00")) #colorblindess color theme
FIG.pH.DAYS.15.to.21 # view figure

# temp.FIGURE
FIG.temp.DAYS.15.to.21 <- ggplot(DAYS.15.to.21.hourly.temp, aes(x=datehour, y=mean.temp, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(aes(x = datehour, y = mean.temp, group=Treatment, color=Treatment),cex=1) + #Plot points using time as the x axis, light as the Y axis and black dots
  geom_errorbar(aes(x=datehour, ymax=mean.temp+se.temp, ymin=mean.temp-se.temp), position=position_dodge(0.9), data=DAYS.15.to.21.hourly.temp, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Days 15-21 hourly temperature (C)") + #Label the graph with the main title
  ylim(15,20) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("Temperature (C)") + #Label the Y Axis
  theme_bw() + #Set the background color
  geom_line() + # add line to plot
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
FIG.temp.DAYS.15.to.21 <- FIG.temp.DAYS.15.to.21 + scale_color_manual(values=c("#56B4E9",  "#E69F00")) #colorblindess color theme
FIG.temp.DAYS.15.to.21 # view figure

#compile figures and output 
FIG.FINAL.pH.temp <- grid.arrange(arrangeGrob(FIG.pH.Apex.heathstack_2, FIG.temp.Apex.heathstack_2,FIG.pH.DAYS.1.to.7, FIG.temp.DAYS.1.to.7,FIG.pH.DAYS.8.to.14,
                                              FIG.temp.DAYS.8.to.14, FIG.pH.DAYS.15.to.21, FIG.temp.DAYS.15.to.21, ncol=2, nrow =4))
FIG.FINAL.pH.temp # view figure

# output figure
ggsave(file="Output/Fig.pH.temp.pdf", FIG.FINAL.pH.temp, width = 12, height = 20, units = c("in"))

# ---------------------------------- Discrete chemistry table and figures =---------------------------- #
# call chem tables to get the mean and standard error discrete measurements 
path.p<-"Output/chem.tables/" #the location of all your respirometry files 
file.names.full<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names
df_total <- data.frame() # start dataframe 
chem.final <- data.frame(matrix(nrow = 1, ncol = 12)) # create dataframe to save cumunalitively during for loop
colnames(chem.final)<-c('date', 'Treatment', 'mean.Salinity', 'sem.Salinity', 'mean.Temperature', 'sem.Temperature' 
                        , 'mean.pH', 'sem.pH', 'mean.pCO2',    'sem.pCO2', 'mean.Aragonite.Sat', 'sem.Aragonite.Sat') # names for comuns in the for loop

for(i in 1:length(file.names.full)) { # for every file in list start at the first and run this following function
  chem.data <-read.table(file.path(path.p,file.names.full[i]), header=T, sep=",", na.string="NA", fill = TRUE, as.is=TRUE) #reads in the data files
  chem.data$date <- substr(file.names.full[i], 1,8)
  chem.data <- chem.data %>% 
    dplyr::select(date, Treatment, mean.Salinity, sem.Salinity, mean.Temperature, sem.Temperature,  
                  mean.pH, sem.pH, mean.pCO2, sem.pCO2, mean.Aragonite.Sat,  sem.Aragonite.Sat)
  df <- data.frame(chem.data) # name dataframe for this singl e row
  chem.final <- rbind(chem.final,df) #bind to a cumulative list dataframe
  print(chem.final) # print to monitor progress
}
chem.final # view table
chem.daily.discrete.TANK <- chem.final %>% 
  filter(Treatment %in% c("tank.Ambient", "tank.Severe", "tank.Moderate")) %>% 
  filter(date > 20190723)

# divide data into the three experimental periods
# days 1 - 7
DISCRETE.1.to.7 <-  chem.daily.discrete.TANK %>% 
  dplyr::filter(date %in% (20190724:20190731))  
# days 8 - 14
DISCRETE.8.to.14 <-  chem.daily.discrete.TANK %>% 
  dplyr::filter(date %in% (20190801:20190807)) 
# days 15 - 21
DISCRETE.15.to.21 <-  chem.daily.discrete.TANK %>% 
  dplyr::filter(date %in% (20190808:20190816)) 
tail(DISCRETE.15.to.21)
 
############################################################################################################################# #
# Assemble the mean and st error for the entire exposure periods ############################################################ #
# days 1 - 7 Amb v. Sev. vs. Mod ############################################################################################ #

D.1.7.LONG <- melt(DISCRETE.1.to.7, id.vars=c("date", "Treatment")) # uses tidyr to make a long table from wide
D.1.7.LONG.2 <- D.1.7.LONG %>%  filter(variable %in% c("mean.Salinity", "mean.Temperature", "mean.pH", "mean.pCO2", "mean.Aragonite.Sat")) # filter only  mean values for summary table

D.1.7_MEAN_ST.E <- ddply(D.1.7.LONG.2, c("Treatment","variable"), summarise, #Calculate descriptive stats by Treatment and exposure 
                        N = length(na.omit(value)), #count the sample size removing NA
                        mean = mean(value), #calculate average 
                        sem = sd(value)/sqrt(N)) # calcualte the standard error of the mean

D.1.7_MEAN_ST.E$period <- "DAYS.1.TO.7" # assign column for the exposure period
D.1.7_MEAN_ST.E # view table

# days 8 - 14 - Amb ######################################################################################################### #

D.8.14.LONG <- melt(DISCRETE.8.to.14, id.vars=c("date", "Treatment")) # uses tidyr to make a long table from wide
D.8.14.LONG.2 <- D.8.14.LONG %>%  filter(variable %in% c("mean.Salinity", "mean.Temperature", "mean.pH", "mean.pCO2", "mean.Aragonite.Sat")) # filter only  mean values for summary table

D.8.14_MEAN_ST.E <- ddply(D.8.14.LONG.2, c("Treatment","variable"), summarise, #Calculate descriptive stats by Treatment and exposure 
                        N = length(na.omit(value)), #count the sample size removing NA
                        mean = mean(value), #calculate average 
                        sem = sd(value)/sqrt(N)) # calcualte the standard error of the mean

D.8.14_MEAN_ST.E$period <- "DAYS.8.TO.14" # assign column for the exposure period
D.8.14_MEAN_ST.E # view table

# days 15 - 21 - Amb vs. Mod ################################################################################################ #

D.15.to.21.LONG <- melt(DISCRETE.15.to.21, id.vars=c("date", "Treatment")) # uses tidyr to make a long table from wide
D.15.to.21.LONG.2 <- D.15.to.21.LONG %>%  filter(variable %in% c("mean.Salinity", "mean.Temperature", "mean.pH", "mean.pCO2", "mean.Aragonite.Sat")) # filter only  mean values for summary table

D.15.21_MEAN_ST.E <- ddply(D.15.to.21.LONG.2, c("Treatment","variable"), summarise, #Calculate descriptive stats by Treatment and exposure 
                        N = length(na.omit(value)), #count the sample size removing NA
                        mean = mean(value), #calculate average 
                        sem = sd(value)/sqrt(N)) # calcualte the standard error of the mean

D.15.21_MEAN_ST.E$period <- "DAYS.15.TO.21" # assign column for the exposure period
D.15.21_MEAN_ST.E # view table

#################################### #
###################### DAYS 1 - 7 #################################### #
#################################### #
 
pd <- position_dodge(0.5) # dodge between treatments to ease asthetics and interpretation
                     
# PH.FIGURE DAYS 1 - 7
PH.FIG.1.7 <- ggplot(DISCRETE.1.to.7, aes(x=date, y=mean.pH, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.pH+sem.pH, ymin=mean.pH-sem.pH), data=DISCRETE.1.to.7, col="black",  width=.1, position=pd) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle(" pH (total scale)") + #Label the graph with the main title
  geom_line(position=pd) + # add line to plot
  ylim(6.8,8) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("pH (total scale)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
PH.FIG.1.7 <- PH.FIG.1.7 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
PH.FIG.1.7 # view figure

# pCO2.FIGURE DAYS 1 - 7
pCO2.FIG.1.7 <- ggplot(DISCRETE.1.to.7, aes(x=date, y=mean.pCO2, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.pCO2+sem.pCO2, ymin=mean.pCO2-sem.pCO2), position=pd, data=DISCRETE.1.to.7, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("pCO2") + #Label the grapCO2 with the main title
  geom_line(position=pd) + # add line to plot
  ylim(600,6000) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("pCO2") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
pCO2.FIG.1.7 <- pCO2.FIG.1.7 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
pCO2.FIG.1.7 # view figure

# Aragonite.Sat.FIGURE DAYS 1 - 7
Aragonite.Sat.FIG.1.7 <- ggplot(DISCRETE.1.to.7, aes(x=date, y=mean.Aragonite.Sat, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.Aragonite.Sat+sem.Aragonite.Sat, ymin=mean.Aragonite.Sat-sem.Aragonite.Sat), position=pd, data=DISCRETE.1.to.7, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Aragonite.Sat") + #Label the graAragonite.Sat with the main title
  geom_line(position=pd) + # add line to plot
  ylim(0,1.5) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("Aragonite.Sat") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
Aragonite.Sat.FIG.1.7 <- Aragonite.Sat.FIG.1.7 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
Aragonite.Sat.FIG.1.7 # view figure

# Salinity.FIGURE DAYS 1 - 7
Salinity.FIG.1.7 <- ggplot(DISCRETE.1.to.7, aes(x=date, y=mean.Salinity, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.Salinity+sem.Salinity, ymin=mean.Salinity-sem.Salinity), position=pd, data=DISCRETE.1.to.7, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Salinity") + #Label the graSalinity with the main title
  geom_line(position=pd) + # add line to plot
  ylim(29,29.5) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("Salinity") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
Salinity.FIG.1.7 <- Salinity.FIG.1.7 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
Salinity.FIG.1.7 # view figure

# Temperature.FIGURE DAYS 1 - 7
Temperature.FIG.1.7 <- ggplot(DISCRETE.1.to.7, aes(x=date, y=mean.Temperature, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.Temperature+sem.Temperature, ymin=mean.Temperature-sem.Temperature), position=pd, data=DISCRETE.1.to.7, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Temperature") + #Label the graTemperature with the main title
  geom_line(position=pd) + # add line to plot
  ylim(16,20) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("Temperature") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
Temperature.FIG.1.7 <- Temperature.FIG.1.7 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
Temperature.FIG.1.7 # view figure

#compile figures and output 
ALL.FIGS.1.7 <- grid.arrange(arrangeGrob(PH.FIG.1.7, pCO2.FIG.1.7,Aragonite.Sat.FIG.1.7,
                                         Salinity.FIG.1.7, Temperature.FIG.1.7, ncol=2, nrow =3))
ALL.FIGS.1.7 # view figure


###################### DAYS 8 - 14 #################################### #
#################################### #
# PH.FIGURE DAYS 8 - 14
PH.FIG.7.14 <- ggplot(DISCRETE.8.to.14, aes(x=date, y=mean.pH, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.pH+sem.pH, ymin=mean.pH-sem.pH), position=pd, data=DISCRETE.8.to.14, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle(" pH (total scale)") + #Label the graph with the main title
  geom_line(position=pd) + # add line to plot
  ylim(6.8,8) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("pH (total scale)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
PH.FIG.7.14 <- PH.FIG.7.14 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
PH.FIG.7.14 # view figure

# pCO2.FIGURE DAYS 8 - 14
pCO2.FIG.7.14 <- ggplot(DISCRETE.8.to.14, aes(x=date, y=mean.pCO2, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.pCO2+sem.pCO2, ymin=mean.pCO2-sem.pCO2), position=pd, data=DISCRETE.8.to.14, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("pCO2") + #Label the grapCO2 with the main title
  geom_line(position=pd) + # add line to plot
  ylim(600,6000) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("pCO2") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
pCO2.FIG.7.14 <- pCO2.FIG.7.14 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
pCO2.FIG.7.14 # view figure

# Aragonite.Sat.FIGURE DAYS 8 - 14
Aragonite.Sat.FIG.7.14 <- ggplot(DISCRETE.8.to.14, aes(x=date, y=mean.Aragonite.Sat, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.Aragonite.Sat+sem.Aragonite.Sat, ymin=mean.Aragonite.Sat-sem.Aragonite.Sat), position=pd, data=DISCRETE.8.to.14, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Aragonite.Sat") + #Label the graAragonite.Sat with the main title
  geom_line(position=pd) + # add line to plot
  ylim(0,1.5) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("Aragonite.Sat") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
Aragonite.Sat.FIG.7.14 <- Aragonite.Sat.FIG.7.14 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
Aragonite.Sat.FIG.7.14 # view figure

# Salinity.FIGURE DAYS 8 - 14
Salinity.FIG.7.14 <- ggplot(DISCRETE.8.to.14, aes(x=date, y=mean.Salinity, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
    geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.Salinity+sem.Salinity, ymin=mean.Salinity-sem.Salinity), position=pd, data=DISCRETE.8.to.14, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Salinity") + #Label the graSalinity with the main title
  geom_line(position=pd) + # add line to plot
  ylim(28,29.5) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("Salinity") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
Salinity.FIG.7.14 <- Salinity.FIG.7.14 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
Salinity.FIG.7.14 # view figure

# Temperature.FIGURE DAYS 8 - 14
Temperature.FIG.7.14 <- ggplot(DISCRETE.8.to.14, aes(x=date, y=mean.Temperature, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.Temperature+sem.Temperature, ymin=mean.Temperature-sem.Temperature), position=pd, data=DISCRETE.8.to.14, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Temperature") + #Label the graTemperature with the main title
  geom_line(position=pd) + # add line to plot
  ylim(16,20) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("Temperature") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
Temperature.FIG.7.14 <- Temperature.FIG.7.14 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
Temperature.FIG.7.14 # view figure

#compile figures and output 
ALL.FIGS.7.14 <- grid.arrange(arrangeGrob(PH.FIG.7.14, pCO2.FIG.7.14,Aragonite.Sat.FIG.7.14,
                                         Salinity.FIG.7.14, Temperature.FIG.7.14, ncol=2, nrow =3))
ALL.FIGS.7.14 # view figure

###################### DAYS 15 - 21 #################################### #
#################################### #
# PH.FIGURE DAYS 15 - 21
PH.FIG.14.21 <- ggplot(DISCRETE.15.to.21, aes(x=date, y=mean.pH, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.pH+sem.pH, ymin=mean.pH-sem.pH), position=pd, data=DISCRETE.15.to.21, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle(" pH (total scale)") + #Label the graph with the main title
  geom_line(position=pd) + # add line to plot
  ylim(6.8,8) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("pH (total scale)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
PH.FIG.14.21 <- PH.FIG.14.21 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
PH.FIG.14.21 # view figure

# pCO2.FIGURE DAYS 15 - 21
pCO2.FIG.14.21 <- ggplot(DISCRETE.15.to.21, aes(x=date, y=mean.pCO2, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.pCO2+sem.pCO2, ymin=mean.pCO2-sem.pCO2), position=pd, data=DISCRETE.15.to.21, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("pCO2") + #Label the grapCO2 with the main title
    geom_line(position=pd) + # add line to plot
  ylim(600,6000) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("pCO2") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
pCO2.FIG.14.21 <- pCO2.FIG.14.21 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
pCO2.FIG.14.21 # view figure

# Aragonite.Sat.FIGURE DAYS 15 - 21
Aragonite.Sat.FIG.14.21 <- ggplot(DISCRETE.15.to.21, aes(x=date, y=mean.Aragonite.Sat, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.Aragonite.Sat+sem.Aragonite.Sat, ymin=mean.Aragonite.Sat-sem.Aragonite.Sat), position=pd, data=DISCRETE.15.to.21, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Aragonite.Sat") + #Label the graAragonite.Sat with the main title
    geom_line(position=pd) + # add line to plot
  ylim(0,1.5) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("Aragonite.Sat") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
Aragonite.Sat.FIG.14.21 <- Aragonite.Sat.FIG.14.21 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
Aragonite.Sat.FIG.14.21 # view figure

# Salinity.FIGURE DAYS 15 - 21
Salinity.FIG.14.21 <- ggplot(DISCRETE.15.to.21, aes(x=date, y=mean.Salinity, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.Salinity+sem.Salinity, ymin=mean.Salinity-sem.Salinity), position=pd, data=DISCRETE.15.to.21, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Salinity") + #Label the graSalinity with the main title
    geom_line(position=pd) + # add line to plot
  ylim(28.8,29.5) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("Salinity") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
Salinity.FIG.14.21 <- Salinity.FIG.14.21 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
Salinity.FIG.14.21 # view figure

# Temperature.FIGURE DAYS 15 - 21
Temperature.FIG.14.21 <- ggplot(DISCRETE.15.to.21, aes(x=date, y=mean.Temperature, group=Treatment, color=Treatment)) +#Plot average diurnal cycle of temperature data
  geom_point(position=pd, size=0.5, shape=21, fill="white") + # 21 is filled circle
  geom_errorbar(aes(x=date, ymax=mean.Temperature+sem.Temperature, ymin=mean.Temperature-sem.Temperature), position=pd, data=DISCRETE.15.to.21, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  ggtitle("Temperature") + #Label the graTemperature with the main title
    geom_line(position=pd) + # add line to plot
  ylim(16,20) + #Set Y axis limits
  xlab("Date") + #Label the X Axis
  ylab("Temperature") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 0, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0),
        legend.position="bottom", #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) #Justify the title to the top left
Temperature.FIG.14.21 <- Temperature.FIG.14.21 + scale_color_manual(values=c("#56B4E9",  "#E69F00", "#D55E00")) #colorblindess color theme
Temperature.FIG.14.21 # view figure

#compile figures and output 
ALL.FIGS.14.21 <- grid.arrange(arrangeGrob(PH.FIG.14.21, pCO2.FIG.14.21,Aragonite.Sat.FIG.14.21,
                                          Salinity.FIG.14.21, Temperature.FIG.14.21, ncol=2, nrow =3))
ALL.FIGS.14.21 # view figure


# output figure
ggsave(file="Output/Discrete.Fig.Days.1-7.pdf", ALL.FIGS.1.7, width = 12, height = 8, units = c("in"))
ggsave(file="Output/Discrete.Fig.Days.8-14.pdf", ALL.FIGS.7.14, width = 12, height = 8, units = c("in"))
ggsave(file="Output/Discrete.Fig.Days.15-21.pdf", ALL.FIGS.14.21, width = 12, height = 8, units = c("in"))



