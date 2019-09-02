#Title: Respiration Figs
#Author: Sam Gurr 
#Edited by: Sam Gurr
#Date Last Modified: 20190822
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
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/")

Final.resp.table <- read.csv(file="Data/SDR_data/Final.resp.rates.csv", header=T) #read Size.info data
Final.resp.table.EXP <- Final.resp.table %>% # table for all data on and after 20190723 
  dplyr::filter(row.num > 22) 

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
library(lsmeans)        # Version: 2.27-62, Date/Publication: 2018-05-11, Depends: methods, R (>= 3.2)

model <- aov(lm(rel.resp.COUNT ~ Treatment*Date, data=df_Days.15.21))
summary(model)
shapiro.test(residuals(model))
TukeyHSD(model)

model2 <-aov(lm(rel.resp.COUNT ~ Treatment*Date, data=df_Days.1.14))
summary(model2)
shapiro.test(residuals(model2))
TukeyHSD(model2)

##################################################### #
##################################################### #
########## SUMMARISE FOR MEAN AND STAND ERROR ####### #
##################################################### #
##################################################### #
#df_Days.1.14$Date <- as.character(df_Days.1.14$Date)
#df_Days.1.14$Treatment <- as.character(df_Days.1.14$Treatment)

rel.resp_Days.1.14.final <- df_Days.1.14 %>% 
  dplyr::group_by(Date, Treatment) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,sd,std.error))


rel.resp_Days.15.21.final <- df_Days.15.21 %>% 
  dplyr::group_by(Date, Treatment) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,sd,std.error))
  

##################################################### #
##################################################### #
########## PLOT RELATIVE RESP TO MEAN AMBIENT ####### #
##################################################### #
##################################################### #
names(rel.resp_Days.1.14.final)

pd <- position_dodge(0.5) # dodge between treatments to ease asthetics and interpretation

# days 1 - 14 rel metabolic rate (per individual)
rel.resp_Days.1.14.final$Date <- as.character(rel.resp_Days.1.14.final$Date)
FIGURE.rel.resp.COUNT.Days1.14 <- ggplot(rel.resp_Days.1.14.final, aes(x=Date, y=rel.resp.COUNT_mean, colour=Treatment, group=Treatment)) + 
  geom_errorbar(aes(ymin=rel.resp.COUNT_mean-rel.resp.COUNT_std.error, 
                    ymax=rel.resp.COUNT_mean+rel.resp.COUNT_std.error), colour="black", width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
  xlab("Date") +
  ylab("SMR per individual & relative to Ambient") +
  scale_colour_hue(name="Supplement type",    # Legend label, use darker colors
                   breaks=c("AHM", "AHS", "EHM", "EHS"),
                   labels=c("Amb-Mod", "Amb-Sev", "Elev-Mod", "Elev-Sev"),
                   l=40) +                    # Use darker colors, lightness=40
  ggtitle("Days 1-14 relative metabolic rate") +
  expand_limits(y=0) +                        # Expand y range
  scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
  geom_vline(xintercept = c(3.5), linetype="dotted", color = "black", size=0.5) + # add line between initial exposure nad ambient recovery
  theme_bw() +
  theme(legend.justification=c(1,1),
        legend.position=c(1,1))               # Position legend in bottom right
FIGURE.rel.resp.COUNT.Days1.14  <- FIGURE.rel.resp.COUNT.Days1.14 +  # add text
  annotate(geom="text", x=2, y=6, label="Amb v. Mod. v. Sev (n = 6 treatments)", color="black") +
  annotate(geom="text", x=5, y=6, label="Ambient recovery", color="black")
FIGURE.rel.resp.COUNT.Days1.14 # view plot # view the plot

# days 15- 21 rel metabolic rate (per individual)
rel.resp_Days.15.21.final$Date <- as.character(rel.resp_Days.15.21.final$Date)
FIGURE.rel.resp.COUNT.Days15.21 <- ggplot(rel.resp_Days.15.21.final, aes(x=Date, y=rel.resp.COUNT_mean, colour=Treatment, group=Treatment)) + 
  geom_errorbar(aes(ymin=rel.resp.COUNT_mean-rel.resp.COUNT_std.error, 
                    ymax=rel.resp.COUNT_mean+rel.resp.COUNT_std.error), colour="black", width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
  xlab("Date") +
  ylab("SMR per individual & relative to Ambient") +
  scale_colour_hue(name="Supplement type",    # Legend label, use darker colors
                   breaks=c("AHAM", "AHMM", "AHSM", "EHAM", "EHMM", "EHSM"),
                   labels=c("Amb-Amb-Mod", "Amb-Mod-Mod", "Amb-Sev-Mod", 
                            "Elev-Amb-Mod", "Elev-Mod-Mod", "Elev-Sev-Mod"),
                   l=40) +                    # Use darker colors, lightness=40
  ggtitle("Days 15-21 relative metabolic rate") +
  expand_limits(y=0) +                        # Expand y range
  scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(0,1),
        legend.position=c(0,1))             # Position legend in bottom right
FIGURE.rel.resp.COUNT.Days15.21  <- FIGURE.rel.resp.COUNT.Days15.21 +  # add text
  annotate(geom="text", x=2, y=6, label="Reciprocal Amb v. Mod (n = 12 treatments)", color="black")
FIGURE.rel.resp.COUNT.Days15.21 # view plot # view the plot

# Arrange plots
SMR.relative.plot <- ggarrange(FIGURE.rel.resp.COUNT.Days1.14, FIGURE.rel.resp.COUNT.Days15.21,ncol = 1, nrow = 2) # combine plots 
SMR.relative.plot # view plots
# ggsave to Output
ggsave(file="Output/SMR.relative.plot.pdf", SMR.relative.plot, width = 12, height = 8, units = c("in")) # respiration rate plots

########################################### #
# plot the size and indiv metrics data ############### #
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
  annotate(geom="text", x=6, y=6, 
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
  annotate(geom="text", x=6, y=2.1, 
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
  annotate(geom="text", x=6, y=16, 
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
  annotate(geom="text", x=6, y=500, 
           label="Ambient", color="black")
RESP.plots.MEAN.biovolume # view plot

Resp.plots <- ggarrange(RESP.plots.MEAN.shell,
                        RESP.plots.INDIV, RESP.plots.MEAN.biovolume, ncol= 2, nrow = 2)
Resp.plots # view plots

# ------------------------------ TOTAL PLOTS WITH MEAN SHELL LENGTH -------------------------------------------------#

#######################################
# ELEVATED HISTORY
#######################################

EH.pre.init.amb.PLOT <- ggplot(EH.initial.exposure, aes(x = factor(Date), y = resp.COUNT.µg.L.hr.indiv, fill = Treatment.ID.initial)) + 
  theme_classic() + scale_fill_manual(values=c("orange", "blue", "orange", "red", "orange", "blue", "orange", "red")) +
  geom_boxplot(alpha = 0.5, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted",  color = "black", size=0.5) +
  ylim(0,15) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("SMR (PER INDIVIDUAL)"~µg~O[2]*hr^{-1}*INDIV^{-1}), x=expression("Date"))

EH.pre.init.amb.PLOT <- EH.pre.init.amb.PLOT + ggtitle("ELEVATED HISTORY") +
  annotate(geom="text", x=1, y=6,  label="Pre-exp (trays)", color="black") +
  annotate(geom="text", x=3, y=6,   label="Initial treatment", color="black") +
  annotate(geom="text", x=6, y=6,  label="Ambient", color="black")
EH.pre.init.amb.PLOT # view plot 

EHA.subs.PLOT <- ggplot(EHA.subs.exposure, aes(x = factor(Date), y = resp.COUNT.µg.L.hr.indiv, fill = Treatment.ID.SUBSQ)) + 
  theme_classic() + scale_fill_manual(values=c("blue1", "blue4")) +
  geom_boxplot(alpha = 0.75, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  ylim(0,15) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("SMR (PER INDIVIDUAL)"~µg~O[2]*hr^{-1}*INDIV^{-1}), x=expression("Date"))

EHS.subs.PLOT <- ggplot(EHS.subs.exposure, aes(x = factor(Date), y = resp.COUNT.µg.L.hr.indiv, fill = Treatment.ID.SUBSQ)) + 
  theme_classic() + scale_fill_manual(values=c("red1", "red4")) +
  geom_boxplot(alpha = 0.75, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  ylim(0,15) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("SMR (PER INDIVIDUAL)"~µg~O[2]*hr^{-1}*INDIV^{-1}), x=expression("Date"))

EHM.subs.PLOT <- ggplot(EHM.subs.exposure, aes(x = factor(Date), y = resp.COUNT.µg.L.hr.indiv, fill = Treatment.ID.SUBSQ)) + 
  theme_classic() + scale_fill_manual(values=c("orange1", "orange4")) +
  geom_boxplot(alpha = 0.75, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  ylim(0,15) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("SMR (PER INDIVIDUAL)"~µg~O[2]*hr^{-1}*INDIV^{-1}), x=expression("Date"))

EH.sub.plots <- grid.arrange(EHA.subs.PLOT, EHM.subs.PLOT,EHS.subs.PLOT, nrow = 1)
Elevated.plots.ALL <- grid.arrange(EH.pre.init.amb.PLOT, EH.sub.plots, nrow = 2)

#######################################
# AMBIENT HISTORY
#######################################
AH.pre.init.amb.PLOT <- ggplot(AH.initial.exposure, aes(x = factor(Date), y = resp.COUNT.µg.L.hr.indiv, fill = Treatment.ID.initial)) + 
  theme_classic() + scale_fill_manual(values=c("blue", "blue", "orange", "red", "orange", "blue", "orange", "red")) +
  geom_boxplot(alpha = 0.5, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  geom_vline(xintercept = c(1.5, 4.5), linetype="dotted",  color = "black", size=0.5) +
  ylim(0,15) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("SMR (PER INDIVIDUAL)"~µg~O[2]*hr^{-1}*INDIV^{-1}), x=expression("Date"))

AH.pre.init.amb.PLOT <- AH.pre.init.amb.PLOT + ggtitle("AMBIENT HISTORY") +
  annotate(geom="text", x=1, y=6,  label="Pre-exp (trays)", color="black") +
  annotate(geom="text", x=3, y=6,   label="Initial treatment", color="black") +
  annotate(geom="text", x=6, y=6,  label="Ambient", color="black")
AH.pre.init.amb.PLOT # view plot 

AHA.subs.PLOT <- ggplot(AHA.subs.exposure, aes(x = factor(Date), y = resp.COUNT.µg.L.hr.indiv, fill = Treatment.ID.SUBSQ)) + 
  theme_classic() + scale_fill_manual(values=c("blue1", "blue4")) +
  geom_boxplot(alpha = 0.75, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  ylim(0,15) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("SMR (PER INDIVIDUAL)"~µg~O[2]*hr^{-1}*INDIV^{-1}), x=expression("Date"))

AHS.subs.PLOT <- ggplot(AHS.subs.exposure, aes(x = factor(Date), y = resp.COUNT.µg.L.hr.indiv, fill = Treatment.ID.SUBSQ)) + 
  theme_classic() + scale_fill_manual(values=c("red1", "red4")) +
  geom_boxplot(alpha = 0.75, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  ylim(0,15) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("SMR (PER INDIVIDUAL)"~µg~O[2]*hr^{-1}*INDIV^{-1}), x=expression("Date"))

AHM.subs.PLOT <- ggplot(AHM.subs.exposure, aes(x = factor(Date), y = resp.COUNT.µg.L.hr.indiv, fill = Treatment.ID.SUBSQ)) + 
  theme_classic() + scale_fill_manual(values=c("orange1", "orange4")) +
  geom_boxplot(alpha = 0.75, # color hue
               width=0.5, # boxplot width
               outlier.size=0, # make outliers small
               position = position_dodge(preserve = "single")) + 
  #geom_point(pch = 19, position = position_jitterdodge(0.05), size=1) +
  stat_summary(fun.y=mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),   width = 0.6, size=0.4, linetype = "dashed", 
               position = position_dodge(preserve = "single")) +
  ylim(0,15) +
  theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  labs(y=expression("SMR (PER INDIVIDUAL)"~µg~O[2]*hr^{-1}*INDIV^{-1}), x=expression("Date"))

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
