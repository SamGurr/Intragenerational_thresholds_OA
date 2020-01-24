#Title: Shell.size.Analysis 
#Author: Sam Gurr 
#Edited by: Sam Gurr
#Date Last Modified: 20190628
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') install_github('colin-olito/LoLinR') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 


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
library(ggpubr)
#Required Data files

# Set Working Directory:
# setwd("~/MyProjects/Geoduck_Conditioning/RAnalysis/") #set working
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/")
#Load Sample Info
Size.data <- read.csv(file="Data/Shell_length/20190628_shell_size.csv", header=T) #read sample.info data
Size.data.ALL <- read.csv(file="Data/Shell_length/Shell_length_data.csv", header=T) #read sample.info data
ID.reference.all <- read.csv(file="Data/Tank.ID.reference.subsequent.csv", header=T) #read sample.info data

################################################################################################################### #
# DATA CARPENTRY FOR FIGURES AND STATISTICAL ANALYSIS   ########################################################### #
################################################################################################################### #

# Merge ID.ref with shell size; select desired columns; divide into treatment periods; rbind to one table
# ID ref modifications - filter tank IDs for treatment periods
# ID.reference.D.1.14 <- ID.reference.all %>%  dplyr::filter(Tank.ID == 1:36) %>%  dplyr::select(Tank.ID, INITIAL.TREATMENT.ID)
# ID.reference.D.15.21 <- ID.reference.all  %>%  dplyr::select(Tank.ID, TREATMENT.ID.TOTAL)
# filter and select dates of experiment and ommit NAs
Size.data.EXPERIMENT <- Size.data.ALL %>% dplyr::filter(Date > 20190722) %>%  dplyr::select(Date, Tank.ID, Sw.Condition, Length, Notes, Tank.ID.SPLIT) # select data since 20190723
Size.data.EXPERIMENT <- Size.data.EXPERIMENT %>% filter_at(vars(Tank.ID), any_vars(!is.na(.))) # only ommit NA of Tank.ID 
Size.data.EXPERIMENT # view table
# modify ID ref columns and merge with shell size table = shell_size_data
ID.reference.short <- ID.reference.all %>% dplyr::select(Tank.ID,TREATMENT.ID.TOTAL, Tank.ID.SPLIT) # select columns desired from ID.references to merge with shell size below
shell_size_data <- merge(Size.data.EXPERIMENT,ID.reference.short,by="Tank.ID") # merge shell size data and ID reference
shell_size_data <- shell_size_data %>%  dplyr::select(-Tank.ID.SPLIT.x)
colnames(shell_size_data)[7] <- "Tank.ID.SPLIT"
# modify shell_size_data for targetted columns and create new columns for different treatments
shell_size_data <- shell_size_data %>% dplyr::select(Date, Length, TREATMENT.ID.TOTAL,Tank.ID, Tank.ID.SPLIT)
shell_size_data$Treatment_history <- substr(shell_size_data$TREATMENT.ID.TOTAL, 1,1) # new column for treatment history
shell_size_data$Treatment.EXP_1 <- substr(shell_size_data$TREATMENT.ID.TOTAL, 3,3) # new column for d1-7 exposure
shell_size_data$Treatment.EXP_2 <- substr(shell_size_data$TREATMENT.ID.TOTAL, 4,4) # new column for d 15-21 exposure
# divide shell_size_data into the 4 different periods (pre, first 7 days, second 7 days, third 7 days of 21-day experiment)
# (1) Pre-experiment data
#Size.pre <- shell_size_data %>% dplyr::filter(Date %in% c(20190723,20190724))
Size.pre <- shell_size_data %>% dplyr::filter(Date %in% c(20190724))
Size.pre$Treatment.EXP_1 <- "NA" # first exposure NA - not needed in figure or analysis
Size.pre$Treatment.EXP_2 <- "NA" # second exposure NA - not needed in figure or analysis # IMPORTANT! PRE-EXPERIMENT DATA
# (2 and 3) First 14 days of the expeiment - Same treatments for figure and analsis (hisotyr and Exp_1 NOT Exp_2 during last 7-day period)
Size.D.1.14 <- shell_size_data %>% dplyr::filter(Date %in% 20190725:20190807)  # filter data with six treatments
Size.D.1.14$Treatment.EXP_2 <- "NA" # second exposure NA - not needed in figure or analysis
Size.D.1.7 <- Size.D.1.14 %>% dplyr::filter(Date < 20190801 ) # IMPORTANT! DAYS 1 - 7 DATA
Size.D.8.14 <- Size.D.1.14 %>% dplyr::filter(Date > 20190731) # IMPORTANT! DAYS 8 - 14 DATA
# (4) Last 7 days of the experiment 
Size.D.15.21 <- shell_size_data %>% dplyr::filter(Date > 20190807)  # IMPORTANT! DAYS 15 - 21 DATA
# bind all data to one table 
SizeTableFINAL <- rbind(Size.pre, Size.D.1.14, Size.D.15.21)


################################################################################################################### #
# DATA CARPENTRY FOR AVERAGE SHELL LENTH TO INITIAL VALUES BY BIOLOGICAL REPLICATE (CUP)        ################### #
################################################################################################################### #

# make a table of the mean and standard deviation of the pre experiment sizes by Tank.ID
pre_experiment <- SizeTableFINAL %>% filter(Date %in% 20190724)
tapply(pre_experiment$Length, pre_experiment$Treatment_history, mean) # mean value of Ambient and Elevated animals
PRE.Table <- as.table(tapply(pre_experiment$Length, (paste(pre_experiment$Treatment_history, pre_experiment$Tank.ID.SPLIT, sep ="_")), mean)) # mean value of treatment and tray
PRE.Table.melt <- melt(PRE.Table, id.vars=c("Tank.ID.SPLIT"))
PRE.Table.melt$Tank.ID.SPLIT <- substr(PRE.Table.melt$indices, 3,7)
MERGE_averages <- merge(SizeTableFINAL, PRE.Table.melt, by = "Tank.ID.SPLIT")
MERGE_averages$length.DIFF <- MERGE_averages$Length - MERGE_averages$value  # calculate the length difference
# PLAY OF THE DATA A BIT AND PLOT
#x_after <- MERGE_averages %>%  filter(Date > 20190724)
#x_20190731$Treatment_hist_EXP1 <- paste(a=x_20190731$Treatment_history, x_20190731$Treatment.EXP_1, sep = "")
#ggboxplot(x_after, x = "Date", y = "length.DIFF",  ylab = "length.difference (mm)",  fill = "Treatment.EXP_1",
#          palette = c( "rickandmorty"),add = "jitter", title = "Length_difference", xlab = "Shell length")

Size.DIFF.D.1.14 <- MERGE_averages %>% dplyr::filter(Date %in% 20190725:20190807)  # filter data with six treatments
Size.DIFF.D.1.7 <- Size.DIFF.D.1.14 %>% dplyr::filter(Date < 20190801 ) # IMPORTANT! DAYS 1 - 7 DATA
Size.DIFF.D.8.14 <- Size.DIFF.D.1.14 %>% dplyr::filter(Date > 20190731) # IMPORTANT! DAYS 8 - 14 DATA
Size.DIFF.D.15.21 <- MERGE_averages %>% dplyr::filter(Date > 20190807)  # IMPORTANT! DAYS 15 - 21 DATA

################################################################################################################### #
# STATISTICAL ANALYSIS  ########################################################################################### #
################################################################################################################### #

# PRE-EXPERIMENT T-TEST ------------------------------------------------------------------------------------------- #
library(ggpubr)
# t.test of "pre" data prior to the experiment
t.test(Length~Treatment_history, data=Size.pre) # p-value = 4.491e-05; difference between pCO2 treatment
mean.pre.length <- as.table(tapply(Size.pre$Length, Size.pre$Treatment_history, mean)) # mean value of Ambient and Elevated animals
mean.pre.length # view table
((mean.pre.length[2] - mean.pre.length[1]) / mean.pre.length[2])*100 # percent length diff of larger animals in elevated v. ambient

PRE_dens <- ggdensity(Size.pre, x = "Length", fill = "Treatment_history", title = "Pre_Experiment" ,palette = "jco")
PRE_dens
PRE_violin.DATE <- ggviolin(Size.pre, x = "Treatment_history", y = "Length", fill = "Treatment_history", 
                            palette = c("#FC4E07", "#00AFBB"), add = "none", title = "Pre_Experiment")
PRE_violin.DATE <- PRE_violin.DATE %>% ggadd(c("boxplot", "jitter"), fill ="white")  # Add box plot
PRE_violin.DATE <- ggpar(PRE_violin.DATE, ylim = c(4,12))
PRE_violin.DATE


M <- as.table(tapply(Size.pre$Length, (paste(Size.pre$Treatment_history, Size.pre$Tank.ID)), mean)) # mean value of treatment and tray
M.melt <- melt(M, id.vars=c("TREAT_TANK.ID")) # make into a table
M.melt$TREAT <- substr(M.melt$indices, 1,1) # make new column for just treatment
tapply(M.melt$value, M.melt$TREAT, mean) # get the mean of treat from the mean of tank; NOTE: this is the same as the mean W/O the tank ID above

# look at the data for fun
tapply(Size.D.1.7$Length, (paste(Size.D.1.7$Treatment_history, Size.D.1.7$Treatment.EXP_1, Size.D.1.7$Date)), mean)
tapply(Size.D.8.14$Length, (paste(Size.D.8.14$Treatment_history,  Size.D.8.14$Treatment.EXP_1, Size.D.8.14$Date)), mean)
tapply(Size.D.15.21$Length, (paste(Size.D.15.21$Treatment_history, Size.D.15.21$Date)), mean)


# MODELS FOR DAYS 1 - 7 ------------------------------------------------------------------------------------------- #


# interaction plots
d.1.7.A <- interaction.plot(Size.D.1.7$Treatment_history, Size.D.1.7$Date, Size.D.1.7$Length)
d.1.7.B <- interaction.plot(Size.D.1.7$Treatment_history, Size.D.1.7$Treatment.EXP_1, Size.D.1.7$Length)
d.1.7.C <- interaction.plot(Size.D.1.7$Date, Size.D.1.7$Treatment.EXP_1, Size.D.1.7$Length)
d.1.7.ints <- ggarrange(d.1.7.A,  d.1.7.B, d.1.7.C, ncol = 3, nrow = 1) # combine plots 
d.1.7.ints # view interaction plots
# three way ANOVA treatment and date
Size.D.1.7$Date <- as.factor(Size.D.1.7$Date) # make Date a factor
# model on raw data
threewayanova_D1.7 <- aov(Length ~ Treatment_history*Treatment.EXP_1*Date, data=Size.D.1.7) # run the model
shapiro.test(residuals(threewayanova_D1.7)) # shaprio wilk test of model residuals p = 0.0297; non-normal distribution
hist((residuals(threewayanova_D1.7))) # histogram of model - looks normal - slight left skew
boxplot(residuals(threewayanova_D1.7)) #plot boxplot of residuals - some outliers present
plot(fitted(threewayanova_D1.7),residuals(threewayanova_D1.7)) # plot residuals
qqnorm(residuals(threewayanova_D1.7)) # qqplot - looks normal
# summary from the log model
summary(threewayanova_D1.7) # significant effect of time and marginal effect of treatment history (same as the raw data)
TukeyHSD(threewayanova_D1.7,  conf.level=0.95) # tukey test on the effect of treatment with 95% confidence

# model on transformed data
Size.D.1.7$Length.log <- log(Size.D.1.7$Length)
threewayanova_D1.7.LOG <- aov(Length.log ~ Treatment_history*Treatment.EXP_1*Date, data=Size.D.1.7) # run the model
shapiro.test(residuals(threewayanova_D1.7.LOG)) # shaprio wilk test of model residuals p = 0.8429; normal distribution (log worked)
hist((residuals(threewayanova_D1.7.LOG))) # histogram of model - looks normal
boxplot(residuals(threewayanova_D1.7.LOG)) #plot boxplot of residuals - some outliers present
plot(fitted(threewayanova_D1.7.LOG),residuals(threewayanova_D1.7.LOG)) # plot residuals
qqnorm(residuals(threewayanova_D1.7.LOG)) # qqplot - looks normal
# summary from the log model
summary(threewayanova_D1.7.LOG) # significant effect of time and marginal effect of treatment history (same as the raw data)
TukeyHSD(threewayanova_D1.7.LOG, 'Date', conf.level=0.95) # tukey test on the effect of treatment with 95% confidence
# significant difference between:
# 20190731-20190728 p = 0.0235595 (very similar to the raw data model)
Day1.7_dens.DATE <- ggdensity(Size.D.1.7, x = "Length", fill = "Date", title = "Day 1-7_.by.date" ,palette = "jco")
Day1.7_dens.DATE
Day8.14_violin.DATE <- ggviolin(Size.D.1.7, x = "Date", y = "Length", fill = "Date", 
                                palette = "jco", add = "none", title = "Day 1-7_.by.date")
Day8.14_violin.DATE <- Day8.14_violin.DATE %>% ggadd(c("boxplot", "jitter"), fill ="white")  # Add box plot
Day8.14_violin.DATE <- ggpar(Day8.14_violin.DATE, ylim = c(4,12))
Day8.14_violin.DATE

Day1.7_dens.treathist <- ggdensity(Size.D.1.7, x = "Length", fill = "Treatment_history", title = "Day 1-7_by.treatmenthistory" ,palette = "jco")
Day1.7_dens.treathist
Day1.7_violin.treathist <- ggviolin(Size.D.1.7, x = "Treatment_history", y = "Length", fill = "Treatment_history", 
                                    palette = c("#FC4E07", "#00AFBB"), add = "none", title = "Day 1-7_by.treatmenthistory")
Day1.7_violin.treathist <- Day1.7_violin.treathist %>% ggadd(c("boxplot", "jitter"), fill ="white")  # Add box plot
Day1.7_violin.treathist <- ggpar(Day1.7_violin.treathist, ylim = c(4,12))
Day1.7_violin.treathist

PLOTS.D.1.7 <- ggarrange(Day8.14_violin.DATE,Day1.7_violin.treathist, nrow = 1, ncol = 2)
# MODELS FOR DAYS 8 - 14 ------------------------------------------------------------------------------------------ #

# interaction plots
Size.D.1.14$Date <- as.character(Size.D.1.14$Date) # make date as.character for interaction plots
interaction.plot(Size.D.8.14$Treatment_history, Size.D.8.14$Date, Size.D.8.14$Length)
interaction.plot(Size.D.8.14$Treatment_history, Size.D.8.14$Treatment.EXP_1, Size.D.8.14$Length)
interaction.plot(Size.D.8.14$Date, Size.D.8.14$Treatment.EXP_1, Size.D.8.14$Length)
# three way ANOVA treatment and date
Size.D.8.14$Date <- as.factor(Size.D.8.14$Date) # make datae as.factors for ANOVA model
threewayanova_D8.14 <- aov(Length ~ Treatment_history*Treatment.EXP_1*Date, data=Size.D.8.14) # run the model
shapiro.test(residuals(threewayanova_D8.14)) # shaprio wilk test of model residuals p = 0.0003904;  non-normal distribution
hist((residuals(threewayanova_D8.14))) # histogram of model - looks normal
boxplot(residuals(threewayanova_D8.14)) #plot boxplot of residuals - some outliers present
plot(fitted(threewayanova_D8.14),residuals(threewayanova_D8.14)) # plot residuals
qqnorm(residuals(threewayanova_D8.14)) # qqplot - looks normal
summary(threewayanova_D8.14) # significant effect of time
TukeyHSD(threewayanova_D8.14, 'Date', conf.level=0.95) # tukey test on the effect of treatment with 95% confidence

# model on transformed data
Size.D.8.14$Length.log <- log(Size.D.8.14$Length)
threewayanova_D.8.14.LOG <- aov(Length.log ~ Treatment_history*Treatment.EXP_1*Date, data=Size.D.8.14) # run the model
shapiro.test(residuals(threewayanova_D.8.14.LOG)) # shaprio wilk test of model residuals p = 0.3392; normal distribution (log worked)
hist((residuals(threewayanova_D.8.14.LOG))) # histogram of model - looks normal
boxplot(residuals(threewayanova_D.8.14.LOG)) #plot boxplot of residuals - some outliers present
plot(fitted(threewayanova_D.8.14.LOG),residuals(threewayanova_D.8.14.LOG)) # plot residuals
qqnorm(residuals(threewayanova_D.8.14.LOG)) # qqplot - looks normal
# summary of transformed data
summary(threewayanova_D.8.14.LOG) # significant effect of time
TukeyHSD(threewayanova_D.8.14.LOG, 'Date', conf.level=0.95) # tukey test on the effect of treatment with 95% confidence
# significant difference between:
# 20190804-20190801   0.0011880
# 20190807-20190801   0.0132141
Day8.14_dens.DATE <- ggdensity(Size.D.8.14, x = "Length", fill = "Date", title = "Day 8-14_.by.date" ,palette = "jco")
Day8.14_dens.DATE
Day8.14_violin.DATE <- ggviolin(Size.D.8.14, x = "Date", y = "Length", fill = "Date", 
                                 palette = "jco", add = "none", title = "Day 8-14_.by.date")
Day8.14_violin.DATE <- Day8.14_violin.DATE %>% ggadd(c("boxplot", "jitter"), fill ="white")  # Add box plot
Day8.14_violin.DATE <- ggpar(Day8.14_violin.DATE, ylim = c(4,12))
Day8.14_violin.DATE


# MODELS FOR DAYS 15 - 21  --------------------------------------------------------------------------------------- #

# make date a character to address as a factor in the model
Size.D.15.21$Date <- as.character(Size.D.15.21$Date) 
#interaction plots
interaction.plot(Size.D.15.21$Treatment_history, Size.D.15.21$Date, Size.D.15.21$Length)
interaction.plot(Size.D.15.21$Treatment_history, Size.D.15.21$Treatment.EXP_1, Size.D.15.21$Length)
interaction.plot(Size.D.15.21$Date, Size.D.15.21$Treatment.EXP_1, Size.D.15.21$Length)
# four way ANOVA treatment and date
fourwayanova_D15.21 <-aov(Length ~ Treatment_history*Treatment.EXP_1*Treatment.EXP_2*Date, data=Size.D.15.21)
shapiro.test(residuals(fourwayanova_D15.21)) # shaprio wilk test of model residuals p = 0.0005398;  non-normal distribution
hist((residuals(fourwayanova_D15.21)))# histogram of model - looks normal
boxplot(residuals(fourwayanova_D15.21)) #plot boxplot of residuals - some outliers present
plot(fitted(fourwayanova_D15.21),residuals(fourwayanova_D15.21)) # plot residuals
qqnorm(residuals(fourwayanova_D15.21)) # qqplot - looks normal
# summary of transformed data
summary(fourwayanova_D15.21) # significant effect of time and treatment history
TukeyHSD(fourwayanova_D15.21, conf.level=0.95) # tukey test on the effect of treatment with 95% confidence

# model on transformed data
Size.D.15.21$Length.log <- log(Size.D.15.21$Length)
fourwayanova_D15.21.LOG <- aov(Length.log ~ Treatment_history*Treatment.EXP_1*Treatment.EXP_2*Date, data=Size.D.15.21) # run the model
shapiro.test(residuals(fourwayanova_D15.21.LOG)) # shaprio wilk test of model residuals p = 0.377; normal distribution (log worked)
hist((residuals(fourwayanova_D15.21.LOG))) # histogram of model - looks normal
boxplot(residuals(fourwayanova_D15.21.LOG)) #plot boxplot of residuals - some outliers present
plot(fitted(fourwayanova_D15.21.LOG),residuals(fourwayanova_D15.21.LOG)) # plot residuals
qqnorm(residuals(fourwayanova_D15.21.LOG)) # qqplot - looks normal
# summary of transformed data
summary(fourwayanova_D15.21.LOG) # significant effect of time and treatment history
TukeyHSD(fourwayanova_D15.21.LOG, conf.level=0.95) # tukey test on the effect of treatment with 95% confidence
# 20190814-20190808  0.0000000
# 20190814-20190811  0.0000002
#plots
Day15.21_dens.DATE <- ggdensity(Size.D.15.21, x = "Length.log", fill = "Date", title = "Day 15-21_.by.date" ,palette = "jco")
Day15.21_dens.DATE
Day15.21_violin.DATE <- ggviolin(Size.D.15.21, x = "Date", y = "Length", fill = "Date", 
                                   palette = "jco", add = "none", title = "Day 15-21_by.date")
Day15.21_violin.DATE <- Day15.21_violin.DATE %>% ggadd(c("boxplot", "jitter"), fill ="white")  # Add box plot
Day15.21_violin.DATE <- ggpar(Day15.21_violin.DATE, ylim = c(4,12))
Day15.21_violin.DATE
Day15.21_dens.treathist <- ggdensity(Size.D.8.14, x = "Length.log", fill = "Treatment_history", title = "Day 15-21_by.treatmenthistory" ,palette = "jco")
Day15.21_dens.treathist
Day15.21_violin.treathist <- ggviolin(Size.D.8.14, x = "Treatment_history", y = "Length", fill = "Treatment_history", 
                                   palette = c("#FC4E07", "#00AFBB") , add = "none", title = "Day 15-21_by.treatmenthistory")
Day15.21_violin.treathist <- Day15.21_violin.treathist %>% ggadd(c("boxplot", "jitter"), fill ="white")  # Add box plot
Day15.21_violin.treathist <- ggpar(Day15.21_violin.treathist, ylim = c(4,12))
Day15.21_violin.treathist


PLOT.D.15.21 <- ggarrange(Day15.21_violin.DATE,Day15.21_violin.treathist, nrow = 1, ncol = 2)

TOTAL.plots <- ggarrange(PRE_violin.DATE,PLOTS.D.1.7, Day8.14_violin.DATE, PLOT.D.15.21, nrow = 2, ncol = 2)

###############################################################################  #
###############################################################################  #
##### Models same as above but for the length diff data correction          ###  #
###############################################################################  #
###############################################################################  #

# RUN MODEL ON PERIOD 1 DAY 1 - 7
Size.DIFF.D.1.7$Date <- as.factor(Size.DIFF.D.1.7$Date) # make Date a factor
threewayanova_D1.7.DIFF <- aov(length.DIFF ~ Treatment_history*Treatment.EXP_1*Date, data=Size.DIFF.D.1.7) # run the model
shapiro.test(residuals(threewayanova_D1.7.DIFF)) # shaprio wilk test of model residuals p = 0.07836; normal distribution
hist((residuals(threewayanova_D1.7.DIFF))) # histogram of model - looks normal - slight left skew
boxplot(residuals(threewayanova_D1.7.DIFF)) #plot boxplot of residuals - some outliers present
plot(fitted(threewayanova_D1.7.DIFF),residuals(threewayanova_D1.7.DIFF)) # plot residuals
qqnorm(residuals(threewayanova_D1.7.DIFF)) # qqplot - looks normal
leveneTest(threewayanova_D1.7.DIFF) # p = 0.2889; homogeneity of variance
summary(threewayanova_D1.7.DIFF) # summary of the model
TukeyHSD(threewayanova_D1.7.DIFF)

# RUN MODEL ON PERIOD 2 DAY8 - 14
Size.DIFF.D.8.14$Date <- as.factor(Size.DIFF.D.8.14$Date) # make Date a factor
threewayanova_D8.14.DIFF <- aov(length.DIFF ~ Treatment_history*Treatment.EXP_1*Date, data=Size.DIFF.D.8.14) # run the model
shapiro.test(residuals(threewayanova_D8.14.DIFF)) # shaprio wilk test of model residuals p = 0.002007; non-normal distribution
hist((residuals(threewayanova_D8.14.DIFF))) # histogram of model - looks normal - slight left skew
boxplot(residuals(threewayanova_D8.14.DIFF)) #plot boxplot of residuals - some outliers present
plot(fitted(threewayanova_D8.14.DIFF),residuals(threewayanova_D8.14.DIFF)) # plot residuals
qqnorm(residuals(threewayanova_D8.14.DIFF)) # qqplot - looks normal
leveneTest(threewayanova_D8.14.DIFF) # p = 0.215; homogeneity of variance
summary(threewayanova_D8.14.DIFF) # summary of the model
TukeyHSD(threewayanova_D8.14.DIFF)
# transformation
min(Size.DIFF.D.8.14$length.DIFF) # minimum value is -4.599977
Size.DIFF.D.8.14$length.DIFF.sqred <- (Size.DIFF.D.8.14$length.DIFF + 4.599977)^2
# run model on the transformed data
threewayanova_D8.14.DIFF.sq <- aov(length.DIFF.sqred ~ Treatment_history*Treatment.EXP_1*Date, data=Size.DIFF.D.8.14) # run the model
shapiro.test(residuals(threewayanova_D8.14.DIFF.sq)) # shaprio wilk test of model residuals p = 0.0903; non-normal distribution
hist((residuals(threewayanova_D8.14.DIFF.sq))) # histogram of model - looks normal - slight left skew
boxplot(residuals(threewayanova_D8.14.DIFF.sq)) #plot boxplot of residuals - some outliers present
plot(fitted(threewayanova_D8.14.DIFF.sq),residuals(threewayanova_D8.14.DIFF.sq)) # plot residuals
qqnorm(residuals(threewayanova_D8.14.DIFF.sq)) # qqplot - looks normal
leveneTest(threewayanova_D8.14.DIFF.sq) # p = 0.2456; homogeneity of variance
summary(threewayanova_D8.14.DIFF.sq) # summary of the model - no change in statistical outcome
TukeyHSD(threewayanova_D8.14.DIFF.sq)

# RUN MODEL ON PERIOD 3 DAY 15 - 21
Size.DIFF.D.15.21$Date <- as.factor(Size.DIFF.D.15.21$Date) # make Date a factor
fourwayanova_D.15.21.DIFF <- aov(length.DIFF ~ Treatment_history*Treatment.EXP_1*Treatment.EXP_2*Date, data=Size.DIFF.D.15.21) # run the model
shapiro.test(residuals(fourwayanova_D.15.21.DIFF)) # shaprio wilk test of model residuals p = 0.0003508; non-normal distribution
hist((residuals(fourwayanova_D.15.21.DIFF))) # histogram of model - looks normal - slight left skew
boxplot(residuals(fourwayanova_D.15.21.DIFF)) #plot boxplot of residuals - some outliers present
plot(fitted(fourwayanova_D.15.21.DIFF),residuals(fourwayanova_D.15.21.DIFF)) # plot residuals
qqnorm(residuals(fourwayanova_D.15.21.DIFF)) # qqplot - looks normal
leveneTest(fourwayanova_D.15.21.DIFF) # p = 0.5121; homogeneity of variance
summary(fourwayanova_D.15.21.DIFF) # summary of the model
TukeyHSD(fourwayanova_D.15.21.DIFF)
# transformation
#min(Size.DIFF.D.15.21$length.DIFF) # minimum value is -3.88284
#Size.DIFF.D.15.21$length.DIFF.sqred <- (Size.DIFF.D.15.21$length.DIFF + 3.88284)^(1/2)
# run model on the transformed data
#fourwayanova_D.15.21.DIFF.sq <- aov(length.DIFF.sqred ~ Treatment_history*Treatment.EXP_1*Treatment.EXP_2*Date, data=Size.DIFF.D.15.21) # run the model
#shapiro.test(residuals(fourwayanova_D.15.21.DIFF.sq)) # shaprio wilk test of model residuals p = 0.0003508; non-normal distribution


###############################################################################  #
###############################################################################  #
##### Models for Day 7 and Day 21 - aligned with the TAOC and TP data       ###  #
###############################################################################  #
###############################################################################  #

# PREPARE THE DATA FOR THESE MODELS (TWO WAY AND THREE WAY ANOVAS)
# DAY 7  use the Days.1.7 created in previous script (above) and isolate 20190731
Day_7_length<- MERGE_averages %>%  filter(Date %in% 20190731) # filter data on 20190731
# DAY 21  use the DATA_Days.15.21 created in previous script (above) and isolate 20190814
Day_21_length <- MERGE_averages %>%  filter(Date %in% 20190814) # filter data on 20190814

# RUN THE MODELS ON RAW DATA (NOT CORRECTED FOR INITIAL AVERAGE ON 20190724) ############### #
# DAY 7 TWO WAY ANOVA FOR TREATMENT INITIAL (HISTORY) AND TREATMENT SECONDARY (DAYS 1 -7 PERIOD)
SIZE_DAY7_mod <- aov(Length ~ Treatment_history*Treatment.EXP_1, data = Day_7_length)
summary(SIZE_DAY7_mod) # summary of model
shapiro.test(residuals(SIZE_DAY7_mod)) # p-value = 0.1555; normal via shapiro wilk test
hist((residuals(SIZE_DAY7_mod))) # histogram of residuals
qqnorm(residuals(SIZE_DAY7_mod)) # qqplot
leveneTest(SIZE_DAY7_mod) # p = 0.934; homogenity of variance 
# DAY 21 THREE WAY ANOVA FOR TREAMENT INITIAL (HISTORY) × TREATMENT SECONDARY (D 1-7) × TREATMENT TERTIARY (D 14 - 21)
SIZE_DAY21_mod <- aov(Length ~ Treatment_history*Treatment.EXP_1*Treatment.EXP_2, data = Day_21_length)
summary(SIZE_DAY21_mod) # summary of model
shapiro.test(residuals(SIZE_DAY21_mod)) # p-value = 0.005336; NOT normal via shapiro wilk test
hist((residuals(SIZE_DAY21_mod))) # histogram of residuals RIGHT SKEW
qqnorm(residuals(SIZE_DAY21_mod)) # qqplot LOOK NORMAL THOUGH!
leveneTest(SIZE_DAY21_mod) # p = 0.3649; homogenity of variance!
# TRANSFORM  DAY 21 MODEL TO RESOLVE NORMALITY VIA SHAPIRO WILK TEST
Day_21_length$Length.log <- log(Day_21_length$Length)
SIZE_DAY21_mod.TRANS <- aov(Length.log ~ Treatment_history*Treatment.EXP_1*Treatment.EXP_2, data = Day_21_length)
summary(SIZE_DAY21_mod.TRANS) # summary of model
shapiro.test(residuals(SIZE_DAY21_mod.TRANS)) # p-value = 0.1029; normal via shapiro wilk test
hist((residuals(SIZE_DAY21_mod.TRANS))) # histogram of residuals RESOLVED SKEW
qqnorm(residuals(SIZE_DAY21_mod.TRANS)) # qqplot appears MORE normal than before
leveneTest(SIZE_DAY21_mod.TRANS) # p = 0.3123; homogenity of variance!

# RUN THE MODELS ON LENGTH_DIFF CORRECTED FOR THE INTIAL AVERAGE ON 20190724 ############### #
SIZE.DIFF_DAY7_mod <- aov(length.DIFF ~ Treatment_history*Treatment.EXP_1, data = Day_7_length)
summary(SIZE.DIFF_DAY7_mod) # summary of model
shapiro.test(residuals(SIZE.DIFF_DAY7_mod)) # p-value = 0.4901; normal via shapiro wilk test
hist((residuals(SIZE.DIFF_DAY7_mod))) # histogram of residuals
qqnorm(residuals(SIZE.DIFF_DAY7_mod)) # qqplot
leveneTest(SIZE.DIFF_DAY7_mod) # p = 0.9437; homogenity of variance 

# DAY 21 THREE WAY ANOVA FOR TREAMENT INITIAL (HISTORY) × TREATMENT SECONDARY (D 1-7) × TREATMENT TERTIARY (D 14 - 21)
SIZE.DIFF_DAY21_mod <- aov(length.DIFF ~ Treatment_history*Treatment.EXP_1*Treatment.EXP_2, data = Day_21_length)
summary(SIZE.DIFF_DAY21_mod) # summary of model
shapiro.test(residuals(SIZE.DIFF_DAY21_mod)) # p-value = 0.01335; NOT normal via shapiro wilk test
hist((residuals(SIZE.DIFF_DAY21_mod))) # histogram of residuals slight skew
qqnorm(residuals(SIZE.DIFF_DAY21_mod)) # qqplot lookk pretty normal
leveneTest(SIZE.DIFF_DAY21_mod) # p = 0.4578; homogenity of variance 
TukeyHSD(SIZE.DIFF_DAY21_mod, conf.level=0.95) 

Day_21_length$treat_hist_exp1 <- paste(Day_21_length$Treatment_history, Day_21_length$Treatment.EXP_1, sep = "_") 
averages_DAY21_lengthdiff <- Day_21_length %>% 
  select(treat_hist_exp1, length.DIFF) %>% 
  dplyr::group_by(treat_hist_exp1) %>% 
  dplyr::summarise_each(funs(mean,std.error)) # MAX = E_S and min = E_M 
averages_DAY21_lengthdiff

# TRANSFORM  DAY 21 MODEL TO RESOLVE NORMALITY VIA SHAPIRO WILK TEST
# minimum <- min(Day_21_length$length.DIFF) # minimum value is -2.629837 (corrected for avaerge in each biologicla replicate - many values here are negative)
# Day_21_length$length.DIFF.sqrd <- log(Day_21_length$length.DIFF + ((abs(minimum))+1) ) # cuberoot to resolve normality
# SIZE.DIFF_DAY21_mod.TRANS <- aov(length.DIFF.sqrd ~ Treatment_history*Treatment.EXP_1*Treatment.EXP_2, data = Day_21_length)
# summary(SIZE.DIFF_DAY21_mod.TRANS) # summary of model - still same effect as the untransformed data
# shapiro.test(residuals(SIZE.DIFF_DAY21_mod.TRANS)) # p-value = 0.527; normal via shapiro wilk test
# hist((residuals(SIZE.DIFF_DAY21_mod.TRANS))) # histogram of residuals RESOLVED SKEW
# qqnorm(residuals(SIZE.DIFF_DAY21_mod.TRANS)) # qqplot appears MORE normal than before
# leveneTest(SIZE.DIFF_DAY21_mod.TRANS) # p = 0.8428; homogenity of variance!
# TukeyHSD(SIZE.DIFF_DAY21_mod.TRANS, conf.level=0.95) 
# Treatment_history:Treatment.EXP_1


#plot
#Day_21_length$Treatment.EXP_1 <- factor(Day_21_length$Treatment.EXP_1 , levels = c("A", "M", "S")) # arrange levels allows ggpubr to plot correctly
#Day21_boxplot.lengthdiff <- ggboxplot(Day_21_length, x = "Treatment.EXP_1", y = "length.DIFF", 
#                                      fill = "Treatment_history", palette = c( "#00AFBB", "#FC4E07"), 
#                                      add = "jitter", title = "Day 21 Length Diff")
#Day21_boxplot.lengthdiff <- ggpar(Day21_boxplot.lengthdiff, ylim = c(-3,3))
#Day21_boxplot.lengthdiff

Day21_boxplot.lengthdiff <- ggplot(Day_21_length, aes(x=Treatment.EXP_1, y=length.DIFF, fill = TREATMENT.ID.TOTAL, shape=treat_hist_exp1, colour=TREATMENT.ID.TOTAL)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0, fill = "white") + 
  theme_classic() +
  labs(y=expression("Shell Growth"~(mm)), x=expression("Secondary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treatment_history ))), size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_colour_manual(values=c("skyblue1", "skyblue1", "deepskyblue3", "deepskyblue3", 
                               "blue", "blue","tomato1","tomato1", 
                               "red1", "red1", "firebrick4","firebrick4")) +
  scale_shape_manual(values=c(24,1,2,21,
                              24,1,2,21,
                              24,1,2,21,24,2)) + 
  scale_fill_manual(values=c("white", "skyblue1","white", "deepskyblue3",
                             "white", "blue", "white", "tomato1", 
                             "white", "red1", "white", "firebrick4","white", "firebrick4")) +
  ylim(-2,3) + 
  scale_x_discrete(labels = c("Ambient","Moderate", "Severe")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
Day21_boxplot.lengthdiff 


##################################################### #
###### BOX PLOTS FOR THE RAW DATA    ################ #
########## PLOT MEAN ST ERROR  size  ################ #----------------------------------------------------------------------------------------------- #
##################################################### #
##################################################### #
# DIFF has both the raw length and the difference length data
Size.DIFF.D.1.7  # DAYS 1 - 7 DATA
Size.DIFF.D.8.14  #  DAYS 8 - 14 DATA
Size.DIFF.D.15.21 #  DAYS 15 - 21 DATA

Size.DIFF.D.1.7$Treat.hist.second <- paste(Size.DIFF.D.1.7$Treatment_history, Size.DIFF.D.1.7$Treatment.EXP_1, sept = "")
Size.DIFF.D.8.14$Treat.hist.second <- paste(Size.DIFF.D.1.7$Treatment_history, Size.DIFF.D.1.7$Treatment.EXP_1, sept = "")
# secondary exposure day 1 - 7
# Days.1.7.PLOTbox.length <- ggboxplot(Size.DIFF.D.1.7, x = "Treatment.EXP_1", y = "Length",  ylab = "mm shell length (raw)",
#                                       fill = "Treatment_history",add = "jitter",  
#                                       palette = c("#00AFBB", "#FC4E07"), title = "A. Secondary exposure (raw)")
# Days.1.7.PLOTbox.length

Days.1.7.PLOTbox.length <- ggplot(Size.DIFF.D.1.7, aes(x=Treatment.EXP_1, y=Length, fill = Treatment_history, colour=Treat.hist.second)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0) + 
  theme_classic() +
  labs(y=expression("Shell Growth"~(mm)), x=expression("Secondary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treatment_history))), size = 2, position = position_jitterdodge(jitter.width = 0.1))+
  scale_colour_manual(values=c("skyblue1","deepskyblue3","blue", "tomato1",  "red1", "firebrick4")) +
  scale_shape_manual(values=c(24, 21)) + 
  scale_fill_manual(values=c("white", "white")) +
  ylim(5,11) + 
  scale_x_discrete(labels = c("Ambient","Moderate", "Severe")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
Days.1.7.PLOTbox.length 

# Days.1.7.PLOTbox.lengthDIFF <- ggboxplot(Size.DIFF.D.1.7, x = "Treatment.EXP_1", y = "length.DIFF",  ylab = "mm shell length (averrage corrected)",
#                                           fill = "Treatment_history",add = "jitter", 
#                                           palette = c("#00AFBB", "#FC4E07"), title = "A. Secondary exposure (av. corrected)")
# Days.1.7.PLOTbox.lengthDIFF

Days.1.7.PLOTbox.lengthDIFF <- ggplot(Size.DIFF.D.1.7, aes(x=Treatment.EXP_1, y=length.DIFF, fill = Treatment_history, colour=Treat.hist.second)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0) + 
  theme_classic() +
  labs(y=expression("Shell Growth"~(mm)), x=expression("Secondary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treatment_history))), size = 2, position = position_jitterdodge(jitter.width = 0.1))+
  scale_colour_manual(values=c("skyblue1","deepskyblue3","blue", "tomato1",  "red1", "firebrick4")) +
  scale_shape_manual(values=c(24, 21)) + 
  scale_fill_manual(values=c("white", "white")) +
  ylim(-2,5) + 
  scale_x_discrete(labels = c("Ambient","Moderate", "Severe")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
Days.1.7.PLOTbox.lengthDIFF 

# ambient recovery - effect of secondary exposure d 8 - 14

# Days.8.14.PLOTbox.length <- ggboxplot(Size.DIFF.D.8.14, x = "Treatment.EXP_1", y = "Length",  ylab = "mm shell length (raw)",
#                                         fill = "Treatment_history",add = "jitter",  
#                                         palette = c("#00AFBB", "#FC4E07"), title = "B. Amb recovery Secondary exposure (raw)")
# Days.8.14.PLOTbox.length

Days.8.14.PLOTbox.length <- ggplot(Size.DIFF.D.8.14, aes(x=Treatment.EXP_1, y=Length, fill = Treatment_history, colour=Treat.hist.second)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0) + 
  theme_classic() +
  labs(y=expression("Shell Growth"~(mm)), x=expression("Secondary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treatment_history))), size = 2, position = position_jitterdodge(jitter.width = 0.1))+
  scale_colour_manual(values=c("skyblue1","deepskyblue3","blue", "tomato1",  "red1", "firebrick4")) +
  scale_shape_manual(values=c(24, 21)) + 
  scale_fill_manual(values=c("white", "white")) +
  ylim(5,11) + 
  scale_x_discrete(labels = c("Ambient","Moderate", "Severe")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
Days.8.14.PLOTbox.length 

# Days.8.14.PLOTbox.lengthDIFF <- ggboxplot(Size.DIFF.D.8.14, x = "Treatment.EXP_1", y = "length.DIFF",  ylab = "mm shell length (averrage corrected)",
#                                            fill = "Treatment_history",add = "jitter", 
#                                            palette = c("#00AFBB", "#FC4E07"), title = "B. Amb recovery Secondary exposure (av. corrected)")
# Days.8.14.PLOTbox.lengthDIFF

Days.8.14.PLOTbox.lengthDIFF <- ggplot(Size.DIFF.D.8.14, aes(x=Treatment.EXP_1, y=length.DIFF, fill = Treatment_history, colour=Treat.hist.second)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0) + 
  theme_classic() +
  labs(y=expression("Shell Growth"~(mm)), x=expression("Secondary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treatment_history))), size = 2, position = position_jitterdodge(jitter.width = 0.1))+
  scale_colour_manual(values=c("skyblue1","deepskyblue3","blue", "tomato1",  "red1", "firebrick4")) +
  scale_shape_manual(values=c(24, 21)) + 
  scale_fill_manual(values=c("white", "white")) +
  ylim(-2,5) + 
  scale_x_discrete(labels = c("Ambient","Moderate", "Severe")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
Days.8.14.PLOTbox.lengthDIFF 


# tertiary exposure d 15 - 21


 
 Days.15.21.PLOTbox.length <- ggplot(Size.DIFF.D.15.21, aes(x=Treatment.EXP_2, y=Length, fill = TREATMENT.ID.TOTAL, shape=TREATMENT.ID.TOTAL, colour=TREATMENT.ID.TOTAL)) +
   geom_boxplot(position=position_dodge(0.8), outlier.size = 0, fill = "white") + 
   theme_classic() +
   labs(y=expression("Shell Growth"~(mm)), x=expression("Secondary pCO"[2]~"Exposure")) +
   geom_point((aes(shape = factor(Treatment_history ))), size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
   scale_colour_manual(values=c("skyblue1", "skyblue1", "deepskyblue3", "deepskyblue3", 
                                "blue", "blue","tomato1","tomato1", 
                                "red1", "red1", "firebrick4","firebrick4")) +
   scale_shape_manual(values=c(24,1,2,21,
                               24,1,2,21,
                               24,1,2,21,24,2)) + 
   scale_fill_manual(values=c("white", "red1","white", "red1",
                              "white", "red1", "white", "red1", 
                              "white", "red1", "white", "red1","white", "red1")) +
   ylim(5,11) + 
   scale_x_discrete(labels = c("Ambient","Moderate")) +
   # theme(legend.title = element_blank()) just ommit the legend title
   theme(legend.position = "none")
 Days.15.21.PLOTbox.length

# Days.15.21.PLOTbox.lengthDIFF <- ggboxplot(Size.DIFF.D.15.21, x = "Treatment.EXP_1", y = "length.DIFF",  ylab = "mm shell length (averrage corrected)",
#                                 fill = "Treatment_history",add = "jitter", shape = "Treatment.EXP_2", 
#                                 palette = c("#00AFBB", "#FC4E07"), title = "C.Secondary exposure (av. corrected)")
# Days.15.21.PLOTbox.lengthDIFF

Days.15.21.PLOTbox.lengthDIFF <- ggplot(Size.DIFF.D.15.21, aes(x=Treatment.EXP_2, y=length.DIFF, fill = TREATMENT.ID.TOTAL, shape=TREATMENT.ID.TOTAL, colour=TREATMENT.ID.TOTAL)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0, fill = "white") + 
  theme_classic() +
  labs(y=expression("Shell Growth"~(mm)), x=expression("Secondary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treatment_history ))), size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_colour_manual(values=c("skyblue1", "skyblue1", "deepskyblue3", "deepskyblue3", 
                               "blue", "blue","tomato1","tomato1", 
                               "red1", "red1", "firebrick4","firebrick4")) +
  scale_shape_manual(values=c(24,1,2,21,
                              24,1,2,21,
                              24,1,2,21,24,2)) + 
  scale_fill_manual(values=c("white", "red1","white", "red1",
                             "white", "red1", "white", "red1", 
                             "white", "red1", "white", "red1","white", "red1")) +
  ylim(-2,5) + 
  scale_x_discrete(labels = c("Ambient","Moderate", "Severe")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
Days.15.21.PLOTbox.lengthDIFF 


# EXPERIMENTAL DESIGN KEY------------------------------------- #


schematic.data <- read.csv(file="Data/SDR_data/Experimental.design.data_2.csv", header=T) #read Size.info data
schematic.data$x.val
pd <- position_dodge(0)
schematic.data$x.val <- as.numeric(schematic.data$x.val )
FIG.schematic.data <- ggplot(schematic.data, aes(x=x.val, y=y.val, colour=treat, group=treat)) + 
  geom_line() +
  theme_bw() +
  xlab("Experimental periods") +
  ylab("") +
  scale_colour_hue(name="treat",    # Legend label, use darker colors
                   breaks=c("AHAA", "AHAM", "AHMA", "AHMM", 
                            "AHSA","AHSM", "EHAA", "EHAM", 
                            "EHMA", "EHMM", "EHSA", "EHSM"),
                   labels=c("Amb-Amb-Amb", "Amb-Amb-Mod", "Amb-Mod-Amb", "Amb-Mod-Mod", 
                            "Amb-Sev-Amb", "Amb-Sev-Mod", "Elev-Amb-Amb","Elev-Amb-Mod",
                            "Elev-Mod-Amb", "Elev-Mod-Mod", "Elev-Sev-Amb", "Elev-Sev-Mod"), l=40) + # Use darker colors, lightness=40
  ggtitle("Experimental design schematic") +
  geom_point(data = subset(schematic.data, x.val %in% c(1.00,2.00, 3.00, 4.00)),
             position = pd, shape=c(21, 21, 21, 21, 21, 21, 24, 24, 24, 24, 24, 24,
                                    21, 21, 21, 21, 21, 21, 24, 24, 24, 24, 24, 24,
                                    21, 21, 21, 21, 21, 21, 24, 24, 24, 24, 24, 24,
                                    21, 21, 21, 21, 21, 21, 24, 24, 24, 24, 24, 24), 
             size=c(5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4,
                    5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4,
                    5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4,
                    5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4), 
             fill=c("white", "white", "white", "white", "white", "white",
                    "white", "white", "white", "white", "white", "white",
                    "white", "white", "white", "white", "white", "white",
                    "white", "white", "white", "white", "white", "white",
                    "white", "white", "white", "white", "white", "white",
                    "white", "white", "white", "white", "white", "white",
                    "white", "firebrick4", "white", "red1", "white", "tomato1", 
                    "white", "blue", "white", "deepskyblue3", "white", "skyblue1")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype="dotted", 
             color = "black", size=0.5) # all the colors needed
FIG.schematic.data  <- print(FIG.schematic.data +  
                               scale_colour_manual(values = c("skyblue1", "skyblue1", "deepskyblue3", "deepskyblue3","blue", "blue", "tomato1", "tomato1",  "red1", "red1", "firebrick4", "firebrick4")) +
                               scale_x_discrete(labels=c("1" = "Pre", "2" = "Days 1-7","3" = "Days 8-14","4" = "Days 15-21")) +
                               theme(legend.position="none",panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(),panel.border = element_blank(),
                                     axis.ticks = element_blank(),axis.text.y=element_blank()) +
                               annotate(geom="text", x=c(1,1,1,1, # PRE
                                                         2,2,2,2,2,2,2,2, # DAYS 1 -7
                                                         3,3, # DAYS 8 - 14
                                                         4,4,4,4,4,4,4,4,4,4,4,4,4,4), # DAYS 15 - 21
                                        y=c(0.5, 0.2, 9.5,3.5,  # PRE
                                            0.5, 0.2, 11.7,9.7,7.7,5.7,3.7,1.7, # DAYS 1 -7
                                            0.5, 0.2, # DAYS 8 - 14
                                            0.5, 0.2, 12.2,11.2,10.2,9.2,8.2,7.2,6.2,5.2,4.2,3.2,2.2,1.2), # DAYS 15 - 21
                                        label=c("3-month Conditioning", "(postlarval-juvenile)","EH", "AH", # PRE
                                                "Secondary Exposure Period", "(7 Days)", "EHA","EHM","EHS","AHA","AHM","AHS", # DAYS 1 -7
                                                "Ambient Recovery Period","(7 Days)",  # DAYS 8 - 14
                                                "Tertiary Exposure Period","(7 Days)","EHAM","EHAA","EHMM","EHMA","EHSM","EHSA","AHAM","AHAA","AHMM","AHMA","AHSM","AHSA"), # DAYS 15 - 21
                                        size =3, color="black"))


# ASSEMBLE ALL PLOTS AND OUTPUT
# ASSEMBLE ALL PLOTS AND OUTPUT 
shell.GROWTH.plots<- ggarrange(Days.1.7.PLOTbox.lengthDIFF, 
                               Days.8.14.PLOTbox.lengthDIFF, 
                               Days.15.21.PLOTbox.lengthDIFF, nrow = 1, ncol = 3, widths = c(1, 0.5, 1))
shell.GROWTH.plots.with.schematic<-ggarrange(FIG.schematic.data, shell.GROWTH.plots, nrow = 2, ncol = 1, heights = c(1, 2))

shell.LENGTH.plots <- ggarrange(Days.1.7.PLOTbox.length, 
                                   Days.8.14.PLOTbox.length, 
                                   Days.15.21.PLOTbox.length,  nrow = 1, ncol = 3, widths = c(1, 0.5, 1))
shell.LENGTH.plots.with.schematic<-ggarrange(FIG.schematic.data, shell.LENGTH.plots, nrow = 2, ncol = 1, heights = c(1, 2))

ggsave(file="Output/SHELL.GROWTH_plots.pdf", shell.GROWTH.plots.with.schematic, width = 12, height = 8, units = c("in")) # respiration rate plots
ggsave(file="Output/SHELL.LENGTH_plots.pdf", shell.LENGTH.plots.with.schematic, width = 12, height = 8, units = c("in")) # respiration rate plots






















##################################################### #
###### PLOTS FOR THE RAW DATA        ################ #
########## PLOT MEAN ST ERROR  size  ################ #----------------------------------------------------------------------------------------------- #
##################################################### #
##################################################### #

########## FIGURE CALCULATIONS (MEAN ± SE) ######################################################################### #
#################### Raw length data ################### #

# pre
Size_Days.pre.final <- Size.pre %>% 
  dplyr::filter(Date %in% 20190724) %>% 
  dplyr::select(Date, Length, Treatment_history) %>% 
  dplyr::group_by(Date, Treatment_history) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))
# days 1 - 7
Size.D.1.7$TREATMENT <- paste(Size.D.1.7$Treatment_history, "H", Size.D.1.7$Treatment.EXP_1,sep="")
Size_Size.D.1.7.final <- Size.D.1.7 %>% 
  dplyr::select(Date, Length, TREATMENT) %>% 
  dplyr::group_by(Date, TREATMENT) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))

# days 8 - 14
Size.D.8.14$TREATMENT <- paste(Size.D.8.14$Treatment_history, "H", Size.D.8.14$Treatment.EXP_1,sep="")
Size_Size.D.8.14.final <- Size.D.8.14 %>% 
  dplyr::select(Date, Length, TREATMENT) %>% 
  dplyr::group_by(Date, TREATMENT) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))

# days 15 - 21
Size_Days.15.21.final <- Size.D.15.21 %>% 
  dplyr::select(Date, Length, TREATMENT.ID.TOTAL) %>% 
  dplyr::group_by(Date, TREATMENT.ID.TOTAL) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))
pd <- position_dodge(0.25) # dodge between treatments to ease asthetics and interpretation

# FIGURE 1 ######## # 
# shell length (raw mm) #
Size_Days.pre.final$Date <- as.character(Size_Days.pre.final$Date) # call date as a character
FIGURE.size.pre <- ggplot(Size_Days.pre.final, aes(x=factor(Date), y=mean, colour=Treatment_history, group=Treatment_history)) + 
  geom_errorbar(aes(ymin=mean-std.error, 
                    ymax=mean+std.error), width=.1, position=pd) +
  geom_line(position=pd) +
  theme_bw() +
  xlab("Date") +
  ylab("Shell length (mm)") +
  scale_colour_hue(name="Treatment",    # Legend label, use darker colors
                   breaks=c("A", "E"),
                   labels=c("Amb", "Elev"), l=40) + # Use darker colors, lightness=40
  ggtitle("Pre-experiment Shell length") +
  geom_point(position=pd, shape=c(21, 24), size=3,fill=("white")) + # 21 is filled circle
  #expand_limits(y=0) +                        # Expand y range
  ylim(6.0,8.0) +
  #scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
  theme(legend.justification=c(1,1),
        legend.position=c(1,1))               # Position legend in bottom right
FIGURE.size.pre  <- print(FIGURE.size.pre + scale_colour_manual(values = c("skyblue1", "tomato1"))+ theme(legend.position = "none"))

# FIGURE 2
# days 1 - 7 shell length (raw mm) ) ----------------------------------------------------------------------------------------------- #

Size_Size.D.1.7.final$Date <- as.factor(Size_Size.D.1.7.final$Date) # call date as a character
FIGURE.size_Size.D.1.7 <- ggplot(Size_Size.D.1.7.final, aes(x=factor(Date), y=mean, group=TREATMENT, colour=TREATMENT)) + 
  geom_errorbar(aes(ymin=mean-std.error, 
                    ymax=mean+std.error), width=.1, position=pd) +
  geom_line(position=pd) +
  theme_bw() +
  xlab("Date") +
  ylab("Shell length (mm)") +
  scale_colour_hue(name="Treatment",    # Legend label, use darker colors
                   breaks=c("AHA", "AHM", "AHS", "EHA","EHM", "EHS"),
                   labels=c("Amb-Amb", "Amb-Mod", "Amb-Sev", "Elev-Amb", "Elev-Mod", "Elev-Sev"), l=40) + # Use darker colors, lightness=40
  ggtitle("Days 1-7 Shell length") +
  geom_point(position=pd, shape=c(21, 21, 21, 24, 24, 24,21, 21, 21, 24, 24, 24,21, 21, 21, 24, 24, 24), size=3, fill=("white")) + # 21 is filled circle
  #expand_limits(y=0) + 
  ylim(6.0,8.0) +
  #scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
  theme(legend.justification=c(1,1),
        legend.position=c(1,1))               # Position legend in bottom right
FIGURE.size_Size.D.1.7  <- print(FIGURE.size_Size.D.1.7 + scale_colour_manual(values = c("skyblue1", "deepskyblue3", "blue", "tomato1", "red1", "firebrick4"))+ theme(legend.position = "none"))

# FIGURE 3
# days 8-14 shell length (raw mm)  ----------------------------------------------------------------------------------------------- #

Size_Size.D.8.14.final$Date <- as.character(Size_Size.D.8.14.final$Date)
FIGURE.size_Size.D.8.14 <- ggplot(Size_Size.D.8.14.final, aes(x=factor(Date), y=mean, colour=TREATMENT, group=TREATMENT)) + 
  geom_errorbar(aes(ymin=mean-std.error, 
                    ymax=mean+std.error), width=.1, position=pd) +
  geom_line(position=pd) +
  theme_bw() +
  xlab("Date") +
  ylab("Shell length (mm)") +
  scale_colour_hue(name="Treatment",    # Legend label, use darker colors
                   breaks=c("AHA", "AHM", "AHS", "EHA","EHM", "EHS"),
                   labels=c("Amb-Amb", "Amb-Mod", "Amb-Sev", "Elev-Amb", "Elev-Mod", "Elev-Sev"), l=40) + # Use darker colors, lightness=40
  ggtitle("Days 8-14 Shell length") +
  geom_point(position=pd, shape=c(21, 21, 21, 24, 24, 24,21, 21, 21, 24, 24, 24,21, 21, 21, 24, 24, 24), size=3, fill=("white")) + # 21 is filled circle
  #expand_limits(y=0) +   
  ylim(6.0,8.0) +
  #scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
  theme(legend.justification=c(1,1),
        legend.position=c(1,1))            # Position legend in bottom right
FIGURE.size_Size.D.8.14  <- print(FIGURE.size_Size.D.8.14 + scale_colour_manual(values = c("skyblue1", "deepskyblue3", "blue", "tomato1", "red1", "firebrick4"))+ theme(legend.position = "none"))

# FIGURE 4
# days 15-21 shell length (raw mm)   ----------------------------------------------------------------------------------------------- #

Size_Days.15.21.final$Date <- as.character(Size_Days.15.21.final$Date)
FIGURE.size_Days.15.21 <- ggplot(Size_Days.15.21.final, aes(x=factor(Date), y=mean, colour=TREATMENT.ID.TOTAL, group=TREATMENT.ID.TOTAL, shape = TREATMENT.ID.TOTAL)) + 
  geom_errorbar(aes(ymin=mean-std.error, 
                    ymax=mean+std.error), width=.1, position=pd) +
  geom_line(position=pd, width=.3) +
  theme_bw() +
  xlab("Date") +
  ylab("Shell length (mm)") +
  scale_colour_hue(name="Treatment",    # Legend label, use darker colors
                   breaks=c("AHAA", "AHAM", "AHMA", "AHMM", 
                            "AHSA","AHSM", "EHAA", "EHAM", 
                            "EHMA", "EHMM", "EHSA", "EHSM"),
                   labels=c("Amb-Amb-Amb", "Amb-Amb-Mod", "Amb-Mod-Amb", "Amb-Mod-Mod", 
                            "Amb-Sev-Amb", "Amb-Sev-Mod", "Elev-Amb-Amb","Elev-Amb-Mod",
                            "Elev-Mod-Amb", "Elev-Mod-Mod", "Elev-Sev-Amb", "Elev-Sev-Mod"), l=40) + # Use darker colors, lightness=40
  ggtitle("Days 15-21 Shell length") +
  geom_point(position=pd, shape=c(21, 21, 21, 21, 21, 21, 24, 24, 24, 24, 24, 24,  21, 21, 21, 21, 21, 21, 
                                  24, 24, 24, 24, 24, 24, 21, 21, 21, 21, 21, 21,  24, 24, 24, 24, 24, 24), size=3, 
             fill=c("firebrick4", "white", "red1", "white", "tomato1",
                    "white", "blue", "white", "deepskyblue3", "white", "skyblue1",
                    "white", "firebrick4", "white", "red1", "white", "tomato1",
                    "white", "blue", "white", "deepskyblue3", "white", "skyblue1",
                    "white", "firebrick4", "white", "red1", "white", "tomato1", 
                    "white", "blue", "white", "deepskyblue3", "white", "skyblue1","white")) + # all the colors needed
  #expand_limits(y=0) +  
  ylim(6.0,8.0) +
  #scale_y_continuous(breaks=0:20*4) +                          # Set tick every 4
  theme(legend.justification=c(1,1),legend.position=c(1,1))    # Position legend in bottom right
FIGURE.size_Days.15.21  <- print(FIGURE.size_Days.15.21 + scale_colour_manual(values = c("skyblue1", "skyblue1", "deepskyblue3", "deepskyblue3",
                                                                                         "blue", "blue", "tomato1", "tomato1", 
                                                                                         "red1", "red1", "firebrick4", "firebrick4")) + theme(legend.position = "none"))

##################################################### #
###### PLOTS FOR THE corrected data      ################ #
########## PLOT MEAN ST ERROR  size  ################ #----------------------------------------------------------------------------------------------- #
##################################################### #
##################################################### #

########## FIGURE CALCULATIONS (MEAN ± SE) ######################################################################### #
#################### Raw length data ################### #
Size.DIFF.D.1.7 <- Size.DIFF.D.1.14 %>% dplyr::filter(Date < 20190801 ) # IMPORTANT! DAYS 1 - 7 DATA
Size.DIFF.D.8.14 <- Size.DIFF.D.1.14 %>% dplyr::filter(Date > 20190731) # IMPORTANT! DAYS 8 - 14 DATA
Size.DIFF.D.15.21 <- MERGE_averages %>% dplyr::filter(Date > 20190807)  # IMPORTANT! DAYS 15 - 21 DATA

# days 1 - 7
Size.DIFF.D.1.7$TREATMENT <- paste(Size.DIFF.D.1.7$Treatment_history, "H", Size.DIFF.D.1.7$Treatment.EXP_1,sep="")
Size.DIFF.D.1.7.final <- Size.DIFF.D.1.7 %>% 
  dplyr::select(Date, length.DIFF, TREATMENT) %>% 
  dplyr::group_by(Date, TREATMENT) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))

# days 8 - 14
Size.DIFF.D.8.14$TREATMENT <- paste(Size.DIFF.D.8.14$Treatment_history, "H", Size.DIFF.D.8.14$Treatment.EXP_1,sep="")
Size.DIFF.D.8.14.final <- Size.DIFF.D.8.14 %>% 
  dplyr::select(Date, length.DIFF, TREATMENT) %>% 
  dplyr::group_by(Date, TREATMENT) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))

# days 15 - 21
Size.DIFF.D.15.21.final <- Size.DIFF.D.15.21 %>% 
  dplyr::select(Date, length.DIFF, TREATMENT.ID.TOTAL) %>% 
  dplyr::group_by(Date, TREATMENT.ID.TOTAL) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))
pd <- position_dodge(0.25) # dodge between treatments to ease asthetics and interpretation

# FIGURE 1
# days 1 - 7 shell length corrected per indiv  ----------------------------------------------------------------------------------------------- #

Size.DIFF.D.1.7.final$Date <- as.factor(Size.DIFF.D.1.7.final$Date) # call date as a character
FIGURE.size_Size.DIFF.D.1.7 <- ggplot(Size.DIFF.D.1.7.final, aes(x=factor(Date), y=mean, group=TREATMENT, colour=TREATMENT)) + 
  geom_errorbar(aes(ymin=mean-std.error, 
                    ymax=mean+std.error), width=.1, position=pd) +
  geom_line(position=pd) +
  theme_bw() +
  xlab("Date") +
  ylab(" Shell length (mm, av. corrected)") +
  scale_colour_hue(name="Treatment",    # Legend label, use darker colors
                   breaks=c("AHA", "AHM", "AHS", "EHA","EHM", "EHS"),
                   labels=c("Amb-Amb", "Amb-Mod", "Amb-Sev", "Elev-Amb", "Elev-Mod", "Elev-Sev"), l=40) + # Use darker colors, lightness=40
  ggtitle("Days 1-7 Shell length  (mm, av. corrected") +
  geom_point(position=pd, shape=c(21, 21, 21, 24, 24, 24,21, 21, 21, 24, 24, 24,21, 21, 21, 24, 24, 24), size=3, fill=("white")) + # 21 is filled circle
  #expand_limits(y=0) + 
  ylim(-0.2, 2) +
  #scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
  theme(legend.justification=c(1,1),
        legend.position=c(1,1))               # Position legend in bottom right
FIGURE.size_Size.DIFF.D.1.7  <- print(FIGURE.size_Size.DIFF.D.1.7 + scale_colour_manual(values = c("skyblue1", "deepskyblue3", "blue", "tomato1", "red1", "firebrick4"))+ theme(legend.position = "none"))
FIGURE.size_Size.DIFF.D.1.7

# FIGURE 2
# days 8-14 shell length corrected  ----------------------------------------------------------------------------------------------- #

Size.DIFF.D.8.14.final$Date <- as.character(Size.DIFF.D.8.14.final$Date)
FIGURE.size_Size.DIFF.D.8.14 <- ggplot(Size.DIFF.D.8.14.final, aes(x=factor(Date), y=mean, colour=TREATMENT, group=TREATMENT)) + 
  geom_errorbar(aes(ymin=mean-std.error, 
                    ymax=mean+std.error), width=.1, position=pd) +
  geom_line(position=pd) +
  theme_bw() +
  xlab("Date") +
  ylab("Shell length (mm, av. corrected)") +
  scale_colour_hue(name="Treatment",    # Legend label, use darker colors
                   breaks=c("AHA", "AHM", "AHS", "EHA","EHM", "EHS"),
                   labels=c("Amb-Amb", "Amb-Mod", "Amb-Sev", "Elev-Amb", "Elev-Mod", "Elev-Sev"), l=40) + # Use darker colors, lightness=40
  ggtitle("Days 8-14 Shell length  (mm, av. corrected") +
  geom_point(position=pd, shape=c(21, 21, 21, 24, 24, 24,21, 21, 21, 24, 24, 24,21, 21, 21, 24, 24, 24), size=3, fill=("white")) + # 21 is filled circle
  #expand_limits(y=0) +   
  ylim(-0.2, 2) +
  #scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
  theme(legend.justification=c(1,1),
        legend.position=c(1,1))            # Position legend in bottom right
FIGURE.size_Size.DIFF.D.8.14  <- print(FIGURE.size_Size.DIFF.D.8.14 + scale_colour_manual(values = c("skyblue1", "deepskyblue3", "blue", "tomato1", "red1", "firebrick4"))+ theme(legend.position = "none"))

# FIGURE 3
# days 15-21 shell length corrected   ----------------------------------------------------------------------------------------------- #

Size.DIFF.D.15.21.final$Date <- as.character(Size.DIFF.D.15.21.final$Date)
FIGURE.size.DIFF.D.15.21 <- ggplot(Size.DIFF.D.15.21.final, aes(x=factor(Date), y=mean, colour=TREATMENT.ID.TOTAL, group=TREATMENT.ID.TOTAL, shape = TREATMENT.ID.TOTAL)) + 
  geom_errorbar(aes(ymin=mean-std.error, 
                    ymax=mean+std.error), width=.1, position=pd) +
  geom_line(position=pd, width=.3) +
  theme_bw() +
  xlab("Date") +
  ylab("Shell length (mm, av. corrected)") +
  scale_colour_hue(name="Treatment",    # Legend label, use darker colors
                   breaks=c("AHAA", "AHAM", "AHMA", "AHMM", 
                            "AHSA","AHSM", "EHAA", "EHAM", 
                            "EHMA", "EHMM", "EHSA", "EHSM"),
                   labels=c("Amb-Amb-Amb", "Amb-Amb-Mod", "Amb-Mod-Amb", "Amb-Mod-Mod", 
                            "Amb-Sev-Amb", "Amb-Sev-Mod", "Elev-Amb-Amb","Elev-Amb-Mod",
                            "Elev-Mod-Amb", "Elev-Mod-Mod", "Elev-Sev-Amb", "Elev-Sev-Mod"), l=40) + # Use darker colors, lightness=40
  ggtitle("Days 15-21 Shell length (mm, av. corrected)") +
  geom_point(position=pd, shape=c(21, 21, 21, 21, 21, 21, 24, 24, 24, 24, 24, 24,  21, 21, 21, 21, 21, 21, 
                                  24, 24, 24, 24, 24, 24, 21, 21, 21, 21, 21, 21,  24, 24, 24, 24, 24, 24), size=3, 
             fill=c("firebrick4", "white", "red1", "white", "tomato1",
                    "white", "blue", "white", "deepskyblue3", "white", "skyblue1",
                    "white", "firebrick4", "white", "red1", "white", "tomato1",
                    "white", "blue", "white", "deepskyblue3", "white", "skyblue1",
                    "white", "firebrick4", "white", "red1", "white", "tomato1", 
                    "white", "blue", "white", "deepskyblue3", "white", "skyblue1","white")) + # all the colors needed
  #expand_limits(y=0) +  
  ylim(-0.2, 2) +
  #scale_y_continuous(breaks=0:20*4) +                          # Set tick every 4
  theme(legend.justification=c(1,1),legend.position=c(1,1))    # Position legend in bottom right
FIGURE.size.DIFF.D.15.21  <- print(FIGURE.size.DIFF.D.15.21 + scale_colour_manual(values = c("skyblue1", "skyblue1", "deepskyblue3", "deepskyblue3",
                                                                                         "blue", "blue", "tomato1", "tomato1", 
                                                                                         "red1", "red1", "firebrick4", "firebrick4")) + theme(legend.position = "none"))


# EXPERIMENTAL DESIGN KEY/SCHEMATIC FOR SMR PLOT ------------------------------------- #



schematic.data <- read.csv(file="Data/SDR_data/Experimental.design.data_2.csv", header=T) #read Size.info data
schematic.data$x.val
pd <- position_dodge(0)
schematic.data$x.val <- as.numeric(schematic.data$x.val )
FIG.schematic.data <- ggplot(schematic.data, aes(x=x.val, y=y.val, colour=treat, group=treat)) + 
  geom_line() +
  theme_bw() +
  xlab("Experimental periods") +
  ylab("") +
  scale_colour_hue(name="treat",    # Legend label, use darker colors
                   breaks=c("AHAA", "AHAM", "AHMA", "AHMM", 
                            "AHSA","AHSM", "EHAA", "EHAM", 
                            "EHMA", "EHMM", "EHSA", "EHSM"),
                   labels=c("Amb-Amb-Amb", "Amb-Amb-Mod", "Amb-Mod-Amb", "Amb-Mod-Mod", 
                            "Amb-Sev-Amb", "Amb-Sev-Mod", "Elev-Amb-Amb","Elev-Amb-Mod",
                            "Elev-Mod-Amb", "Elev-Mod-Mod", "Elev-Sev-Amb", "Elev-Sev-Mod"), l=40) + # Use darker colors, lightness=40
  ggtitle("Experimental design schematic") +
  geom_point(data = subset(schematic.data, x.val %in% c(1.00,2.00, 3.00, 4.00)),
             position = pd, shape=c(21, 21, 21, 21, 21, 21, 24, 24, 24, 24, 24, 24,
                                    21, 21, 21, 21, 21, 21, 24, 24, 24, 24, 24, 24,
                                    21, 21, 21, 21, 21, 21, 24, 24, 24, 24, 24, 24,
                                    21, 21, 21, 21, 21, 21, 24, 24, 24, 24, 24, 24), size=3, 
             fill=c("white", "white", "white", "white", "white", "white",
                    "white", "white", "white", "white", "white", "white",
                    "white", "white", "white", "white", "white", "white",
                    "white", "white", "white", "white", "white", "white",
                    "white", "white", "white", "white", "white", "white",
                    "white", "white", "white", "white", "white", "white",
                    "white", "firebrick4", "white", "red1", "white", "tomato1", 
                    "white", "blue", "white", "deepskyblue3", "white", "skyblue1")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype="dotted", 
             color = "black", size=0.5) # all the colors needed
FIG.schematic.data  <- print(FIG.schematic.data +  
                               scale_colour_manual(values = c("skyblue1", "skyblue1", "deepskyblue3", "deepskyblue3","blue", "blue", "tomato1", "tomato1",  "red1", "red1", "firebrick4", "firebrick4")) +
                               scale_x_discrete(labels=c("1" = "Pre", "2" = "Days 1-7","3" = "Days 8-14","4" = "Days 15-21")) +
                               theme(legend.position="none",panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(),panel.border = element_blank(),
                                     axis.ticks = element_blank(),axis.text.y=element_blank()) +
                               annotate(geom="text", x=c(1,1,1,1, # PRE
                                                         2,2,2,2,2,2,2,2, # DAYS 1 -7
                                                         3,3, # DAYS 8 - 14
                                                         4,4,4,4,4,4,4,4,4,4,4,4,4,4), # DAYS 15 - 21
                                        y=c(0.5, 0.2, 9.5,3.5,  # PRE
                                            0.5, 0.2, 11.7,9.7,7.7,5.7,3.7,1.7, # DAYS 1 -7
                                            0.5, 0.2, # DAYS 8 - 14
                                            0.5, 0.2, 12.2,11.2,10.2,9.2,8.2,7.2,6.2,5.2,4.2,3.2,2.2,1.2), # DAYS 15 - 21
                                        label=c("3-month Conditioning", "(postlarval-juvenile)","EH", "AH", # PRE
                                                "Secondary Exposure Period", "(7 Days)", "EHA","EHM","EHS","AHA","AHM","AHS", # DAYS 1 -7
                                                "Ambient Recovery Period","(7 Days)",  # DAYS 8 - 14
                                                "Tertiary Exposure Period","(7 Days)","EHAM","EHAA","EHMM","EHMA","EHSM","EHSA","AHAM","AHAA","AHMM","AHMA","AHSM","AHSA"), # DAYS 15 - 21
                                        size =3, color="black"))

FIGURE.size_Size.DIFF.D.1.7
FIGURE.size_Size.DIFF.D.8.14
FIGURE.size.DIFF.D.15.21
# Arrange plots
SIZE.plot <- ggarrange(FIGURE.size.pre,  FIGURE.size_Size.D.1.7, FIGURE.size_Size.D.8.14, FIGURE.size_Days.15.21,ncol = 4, nrow = 1, widths = c(0.7,2,2,2), labels = c("B","C","D","E")) # combine plots 
SIZE.plot # view plots
SIZE.plot.with.schematic <- ggarrange(FIG.schematic.data, SIZE.plot,ncol = 1, nrow = 2, heights = c(1, 2.0), labels = "A") # combine plots

SIZE.plot.DIFF <- ggarrange(FIGURE.size_Size.DIFF.D.1.7,  FIGURE.size_Size.DIFF.D.8.14, FIGURE.size.DIFF.D.15.21,ncol = 3, nrow = 1, widths = c(2,2,2), labels = c("B","C","D","E")) # combine plots 
SIZE.plot.DIFF # view plots
SIZE.DIFF.plot.with.schematic <- ggarrange(FIG.schematic.data, SIZE.plot.DIFF,ncol = 1, nrow = 2, heights = c(1, 2.0), labels = "A") # combine plots

#gsave to Output
ggsave(file="Output/SIZE.plot.pdf", SIZE.plot.with.schematic, width = 12, height = 8, units = c("in")) # SIZE  plots
ggsave(file="Output/SIZE.plot.average.corrected.pdf", SIZE.DIFF.plot.with.schematic, width = 12, height = 8, units = c("in")) # SIZE  plots


############################################################################################################ #
############ PLOTS MADE FOR SHELL SIZE ON 6/28 ############################################################# #
############################################################################################################ #

# Analysis
ttest_size_20190628 <- t.test(Length ~ Treatment, data = Size.data)
ttest_size_20190628 # view t test results - significant difference between treatments

aov.mod_1 <- aov(Length ~ Treatment+ID, data = Size.data)
anova(aov.mod_1)

aov.mod <- aov(Length ~ ID, data = Size.data)
anova(aov.mod)

TukeyHSD(aov.mod, which = "ID")

ID.size.ANOVA <- lsmeans(aov.mod, pairwise ~ ID)# pariwise Tukey Post-hoc test between repeated treatments
ID.size.ANOVA # view post hoc summary
Treatment.pairs.05 <- cld(ID.size.ANOVA, alpha=.05, Letters=letters) #list pairwise tests and letter display p < 0.05
Treatment.pairs.05 #view results

# plot data
size_plot_20190628 <- ggplot(Size.data, aes(x = factor(ID), y = Length, fill = Treatment)) +
  theme_classic() +
  scale_fill_manual(values=c("blue", "orange", "blue", "orange", "blue", "orange", "blue", "orange"), 
                    #labels=c("Ambient","Ambient × Ambient","Ambient × Elevated","Elevated","Elevated × Ambient","Elevated × Elevated")
  ) +
  geom_boxplot(alpha = 0.5, # color hue
               width=0.6, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  geom_point(pch = 19, position = position_jitterdodge(.3), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), 
               width = 0.6, size=0.4, linetype = "dashed", position = position_dodge(preserve = "single")) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  ylim(0,8) +
  labs(y=expression("Shell size"~(mm)), x=expression("Days"))
size_plot_20190628 # view plot

#save file to output
ggsave(file="Output/size_plot_20190628.pdf", size_plot_20190628, width = 12, height = 8, units = c("in"))
