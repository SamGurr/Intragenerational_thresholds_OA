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
library(ggpubr)
# Set Working Directory:
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/")

# upload the final resp rate data
Final.resp.table <- read.csv(file="Data/SDR_data/Final.resp.rates.csv", header=T) #read Size.info data
DATA <- read.csv(file="Data/SDR_data/Final_table_for_resp_analysis.csv", header=T) #read Size.info data
Final.resp.table.EXP <- Final.resp.table %>% # table for all data on and after 20190723 
  dplyr::filter(row.num > 22) 


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

############################################################################################# #
############################################################################################# #
#########################  RESP RATE DATA ################################################### #
############################################################################################# #
############################################################################################# #
DATA
DATA.pre <- DATA %>% dplyr::filter(Date %in% 20190723) %>% dplyr::select(Date, resp.COUNT.µg.L.hr.indiv,resp.MEAN.µg.L.hr.mm, Treatment.history) # pre experiment 
DATA.Days.1.7 <- DATA %>% dplyr::filter(Date %in% 20190725:20190731) %>% dplyr::select(Date, resp.COUNT.µg.L.hr.indiv,resp.MEAN.µg.L.hr.mm, Treatment.EXP_1, Treatment.history) # d 1 - 7 experiment
DATA_Days.8.14 <- DATA %>% dplyr::filter(Date %in% 20190801:20190807) %>% dplyr::select(Date, resp.COUNT.µg.L.hr.indiv,resp.MEAN.µg.L.hr.mm, Treatment.EXP_1, Treatment.history) # d 8 - 14 experiment
DATA_Days.15.21 <- DATA %>% dplyr::filter(Date %in% 20190808:20190814) %>% dplyr::select(Date, resp.COUNT.µg.L.hr.indiv,resp.MEAN.µg.L.hr.mm, Treatment.EXP_2, Treatment.history, Treatment.EXP_1) # d 15-21 experiment

########## mean and st error calculations ########################### #
# pre
resp_Days.pre.final <- DATA.pre %>% 
  dplyr::group_by(Date, Treatment.history) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))
# days 1 - 7
resp_Days.1.7.final <- DATA.Days.1.7 %>% 
  dplyr::group_by(Date, Treatment.EXP_1,Treatment.history) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))
resp_Days.1.7.final$TREATMENT <- paste(resp_Days.1.7.final$Treatment.history,resp_Days.1.7.final$Treatment.EXP_1, sep="") # IMPORTANT - this combines history and first treatment exposure for figure treatments
# days 8 - 14
resp_Days.8.14.final <- DATA_Days.8.14 %>% 
  dplyr::group_by(Date, Treatment.EXP_1, Treatment.history) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))
resp_Days.8.14.final$TREATMENT <- paste(resp_Days.8.14.final$Treatment.history,resp_Days.8.14.final$Treatment.EXP_1, sep="") # IMPORTANT - this combines history and first treatment exposure for figure treatments
# days 15 - 21
resp_Days.15.21.final <- DATA_Days.15.21 %>% 
  dplyr::group_by(Date, Treatment.EXP_2, Treatment.EXP_1, Treatment.history, Treatment.EXP_1, Treatment.EXP_2) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))
resp_Days.15.21.final$TREATMENT <- paste(resp_Days.15.21.final$Treatment.history,resp_Days.15.21.final$Treatment.EXP_1,resp_Days.15.21.final$Treatment.EXP_2, sep="") # IMPORTANT - this combines history and first treatment exposure for figure treatments


##################################################### #
##################################################### #
########## PLOTS OF THE SIG EFFECTS  ################ #----------------------------------------------------------------------------------------------- #
##################################################### #
##################################################### #

# FIRST 7 DAYS (SECONDARY EXPOSURE) - SIG DIFFERENCE BTWN EHM AND AHM
DATA.pre$Treat <- substr(DATA.pre$Treatment.history, 1,1)

DATA.pre.PLOT <- ggviolin(DATA.pre, x = "Treat", y = "resp.COUNT.µg.L.hr.indiv",  ylab = "µg.L.hr.indiv",  fill = "Treat",
                          palette = c("#FC4E07", "#00AFBB"),add = "none", title = "A.Respiration_pre_experiment")
DATA.pre.PLOT <- DATA.pre.PLOT %>% ggadd(c("boxplot", "jitter"),shape ="Treat", fill = "white") # Add box plot
#DATA.pre.PLOT <- ggpar(DATA.pre.PLOT, ylim = c(-2,35))
DATA.pre.PLOT

DATA.pre.PLOTbox<- ggboxplot(DATA.pre, x = "Treat", y = "resp.COUNT.µg.L.hr.indiv",  ylab = "µg.L.hr.indiv",  fill = "Treat",
                             palette = c("#FC4E07", "#00AFBB"),add = "jitter", title = "Pre-conditioning", xlab = "Initial pCO2 Treatment")
DATA.pre.PLOTbox <- ggpar(DATA.pre.PLOTbox, ylim = c(-2,18))
DATA.pre.PLOTbox

DATA.pre.PLOTbox <- ggplot(DATA.pre, aes(x=Treat, y=resp.COUNT.µg.L.hr.indiv, fill =Treat, colour=Treat)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0, fill = "white") + 
  theme_classic() +
  labs(y=expression("Respiration rate"~(~µg~O[2]*hr^{-1}*individual^{-1})), x=expression("Secondary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treat))), size = 2, position = position_jitterdodge(jitter.width = 0.5))+
  scale_colour_manual(values=c("skyblue1","tomato1")) +
  scale_shape_manual(values=c(24, 21)) + 
  scale_fill_manual(values=c('white', 'white')) +
  ylim(0,20) + 
  scale_x_discrete(labels = c("Ambient","Elevated")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
DATA.pre.PLOTbox # view the plot

# TREATMENT HISTORY (INITIAL EXPOSURE) - NO DIFFERENCE BTWN TREATMENT
DATA.Days.1.7$Treat.initial <- substr(DATA.Days.1.7$Treatment.history, 1,1)
DATA.Days.1.7$treat.hist.secondary <- paste(DATA.Days.1.7$Treatment.history, DATA.Days.1.7$Treatment.EXP_1, sep ="")

Days.1.7.PLOTbox <- ggplot(DATA.Days.1.7, aes(x=Treatment.EXP_1, y=resp.COUNT.µg.L.hr.indiv, fill =factor(Treatment.history), colour=treat.hist.secondary)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0, fill = "white") + 
  stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = 18, size = 3,
               show.legend = FALSE) +
  theme_classic() +
  labs(y=expression("Respiration rate"~(~µg~O[2]*hr^{-1}*individual^{-1})), x=expression("Secondary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treatment.history))), size = 2, position = position_jitterdodge(jitter.width = 0.5))+
  scale_colour_manual(values=c("skyblue1", "deepskyblue3","blue", "tomato1",  "red1", "firebrick4")) +
  scale_shape_manual(values=c(24, 21)) + 
  scale_fill_manual(values=c('white', 'white')) +
  ylim(0,20) + 
  scale_x_discrete(labels = c("Ambient","Moderate", "Severe")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
Days.1.7.PLOTbox # view the plot
# FACET WRAPPED
Days.1.7.PLOTbox.facet <- ggplot(DATA.Days.1.7, aes(Treatment.EXP_1, resp.COUNT.µg.L.hr.indiv, fill = factor(Treatment.EXP_1))) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0) + 
  stat_summary(fun.y = mean, color = "black", position = position_dodge(0.2),
               geom = "point", shape = 15, size = 4,
               show.legend = FALSE) +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.5)) +
  scale_color_manual(values=c("black", "black", "black")) + 
  scale_fill_manual(values=c("white", "grey54", "grey30")) +  
  labs(title = "Respiration rate",
       subtitle = "(Secondary exposure period)",
       y=expression("Respiration rate"~(~µg~O[2]*hr^{-1}*individual^{-1})), 
       x=expression("Secondary pCO"[2]~"Exposure")) + 
  ylim(0,20) + 
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("A","M", "S")) +
  facet_wrap(~ Treat.initial)
Days.1.7.PLOTbox.facet

# TREATMENT HISTORY (INITIAL EXPOSURE) - NO DIFFERENCE BTWN TREATMENT
DATA_Days.15.21$Treat.initial <- substr(DATA_Days.15.21$Treatment.history, 1,1)
DATA_Days.15.21$treat.TOTAL <- paste(DATA_Days.15.21$Treatment.history, DATA_Days.15.21$Treatment.EXP_1,DATA_Days.15.21$Treatment.EXP_2, sep ="")


Days.15.21.PLOTbox <- ggplot(DATA_Days.15.21, aes(x=Treatment.EXP_2, y=resp.COUNT.µg.L.hr.indiv, fill = factor(treat.TOTAL), shape=treat.TOTAL, colour=treat.TOTAL)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0, fill = "white") + 
  stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = 18, size = 3,
               show.legend = FALSE) +
  theme_classic() +
  labs(y=expression("Respiration rate"~(~µg~O[2]*hr^{-1}*individual^{-1})), x=expression("Tertiary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treatment.history))), size = 2, position = position_jitterdodge(jitter.width = 0.15)) +
  scale_colour_manual(values=c("skyblue1", "skyblue1", "deepskyblue3", "deepskyblue3", 
                               "blue", "blue","tomato1","tomato1", 
                               "red1", "red1", "firebrick4","firebrick4")) +
  scale_shape_manual(values=c(24,1,2,21,
                              24,1,2,21,
                              24,1,2,21,24,2)) + 
  scale_fill_manual(values=c("white", "red1","white", "red1",
                             "white", "red1", "white", "red1", 
                             "white", "red1", "white", "red1","white", "red1")) +
  ylim(0,20) + 
  scale_x_discrete(labels = c("Ambient","Moderate")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
Days.15.21.PLOTbox # view the plot
# FACET WRAPPED
Days.15.21.PLOTbox.facet <- ggplot(DATA_Days.15.21, aes(Treatment.EXP_2, resp.COUNT.µg.L.hr.indiv, fill = factor(Treatment.EXP_1))) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0) + 
  stat_summary(fun.y = mean, color = "black", position = position_dodge(0.75),
               geom = "point", shape = 15, size = 4,
               show.legend = FALSE) +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.5))+
  scale_color_manual(values=c("black", "black", "black")) +  
  scale_fill_manual(values=c("white", "grey54", "grey30")) +  
  labs(title = "Respiration rate",
       subtitle = "(Tertiary exposure period)",
       y=expression("Respiration rate"~(~µg~O[2]*hr^{-1}*individual^{-1})), 
       x=expression("Tertiary pCO"[2]~"Exposure")) + 
  ylim(0,20) + 
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("A","M", "S")) +
  facet_wrap(~ Treat.initial)
Days.15.21.PLOTbox.facet

# ggarrange
Resp.plots <- ggarrange(DATA.pre.PLOTbox, Days.1.7.PLOTbox, Days.15.21.PLOTbox, nrow = 1, widths = c(0.5, 1, 1), labels = c("A","B","C"))
Resp.plots.facet <- ggarrange(Days.1.7.PLOTbox.facet, Days.15.21.PLOTbox.facet, nrow = 1, widths = c(0.5, 1, 1), labels = c("A","B"))

# ggsave
RESP.days.plot.with.schematic <- ggarrange(FIG.schematic.data, Resp.plots,ncol = 1, nrow = 2, heights = c(1, 2.0), labels = "A") # combine plots
ggsave(file="Output/Resp.plots.pdf", RESP.days.plot.with.schematic, width = 12, height = 8, units = c("in")) # respiration rate plots
ggsave(file="Output/Resp.plots.facet.pdf", Resp.plots.facet, width = 12, height = 8, units = c("in")) # respiration rate plots



##################################################### #
############PLOTS WITH TIME ######################### #
########## PLOT MEAN ST ERROR  RESP  ################ #----------------------------------------------------------------------------------------------- #
########### SUPPLEMENTARY PLOTS ##################### #
##################################################### #
names(resp_Days.1.7.final)

pd <- position_dodge(0.25) # dodge between treatments to ease asthetics and interpretation


# FIGURE 1
resp_Days.pre.final$Date <- as.character(resp_Days.pre.final$Date) # call date as a character
FIGURE.resp.pre <- ggplot(resp_Days.pre.final, aes(x=factor(Date), y=resp.COUNT.µg.L.hr.indiv_mean, colour=Treatment.history, group=Treatment.history)) + 
  geom_errorbar(aes(ymin=resp.COUNT.µg.L.hr.indiv_mean-resp.COUNT.µg.L.hr.indiv_std.error, 
                    ymax=resp.COUNT.µg.L.hr.indiv_mean+resp.COUNT.µg.L.hr.indiv_std.error), width=.1, position=pd) +
  geom_line(position=pd) +
  theme_bw() +
  xlab("Date") +
  ylab("SMR per individual") +
  scale_colour_hue(name="Treatment",    # Legend label, use darker colors
                   breaks=c("AHA", "AHM", "AHS", "EHA","EHM", "EHS"),
                   labels=c("Amb-Amb", "Amb-Mod", "Amb-Sev", "Elev-Amb", "Elev-Mod", "Elev-Sev"), l=40) + # Use darker colors, lightness=40
  ggtitle("Pre-experiment SMR") +
  geom_point(position=pd, shape=c(21, 24), size=3,fill=("white")) + # 21 is filled circle
  #expand_limits(y=0) +                        # Expand y range
  ylim(3.5, 13) +
  #scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
  theme(legend.justification=c(1,1),
        legend.position=c(1,1))               # Position legend in bottom right
FIGURE.resp.pre  <- print(FIGURE.resp.pre + scale_colour_manual(values = c("skyblue1", "tomato1"))+ theme(legend.position = "none"))


# FIGURE 2
# days 1 - 7 rel metabolic rate (per individual) ----------------------------------------------------------------------------------------------- #


resp_Days.1.7.final$Date <- as.character(resp_Days.1.7.final$Date) # call date as a character
FIGURE.resp_Days.1.7 <- ggplot(resp_Days.1.7.final, aes(x=factor(Date), y=resp.COUNT.µg.L.hr.indiv_mean, group=TREATMENT, colour=TREATMENT)) + 
  geom_errorbar(aes(ymin=resp.COUNT.µg.L.hr.indiv_mean-resp.COUNT.µg.L.hr.indiv_std.error, 
                    ymax=resp.COUNT.µg.L.hr.indiv_mean+resp.COUNT.µg.L.hr.indiv_std.error), width=.3, position=pd) +
  geom_line(position=pd, width=.3) +
  theme_bw() +
  xlab("Date") +
  ylab("SMR per individual") +
  scale_colour_hue(name="Treatment",    # Legend label, use darker colors
                   breaks=c("AHA", "AHM", "AHS", "EHA","EHM", "EHS"),
                   labels=c("Amb-Amb", "Amb-Mod", "Amb-Sev", "Elev-Amb", "Elev-Mod", "Elev-Sev"), l=40) + # Use darker colors, lightness=40
  ggtitle("Days 1-7 SMR") +
  geom_point(position=pd, shape=c(21, 21, 21, 24, 24, 24,21, 21, 21, 24, 24, 24,21, 21, 21, 24, 24, 24), size=3, fill=("white")) + # 21 is filled circle
  #expand_limits(y=0) + 
  ylim(3.5, 13) +
  #scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
  theme(legend.justification=c(1,1),
        legend.position=c(1,1))               # Position legend in bottom right
FIGURE.resp_Days.1.7  <- print(FIGURE.resp_Days.1.7 + scale_colour_manual(values = c("skyblue1", "deepskyblue3", "blue", "tomato1", "red1", "firebrick4"))+ theme(legend.position = "none"))


# FIGURE 3
# days 8-14 rel metabolic rate (per individual)----------------------------------------------------------------------------------------------- #



resp_Days.8.14.final$Date <- as.character(resp_Days.8.14.final$Date)
FIGURE.resp_Days.8.14 <- ggplot(resp_Days.8.14.final, aes(x=factor(Date), y=resp.COUNT.µg.L.hr.indiv_mean, colour=TREATMENT, group=TREATMENT)) + 
  geom_errorbar(aes(ymin=resp.COUNT.µg.L.hr.indiv_mean-resp.COUNT.µg.L.hr.indiv_std.error, 
                    ymax=resp.COUNT.µg.L.hr.indiv_mean+resp.COUNT.µg.L.hr.indiv_std.error), width=.3, position=pd) +
  geom_line(position=pd, width=.3) +
  theme_bw() +
  xlab("Date") +
  ylab("SMR per individual") +
  scale_colour_hue(name="Treatment",    # Legend label, use darker colors
                   breaks=c("AHA", "AHM", "AHS", "EHA","EHM", "EHS"),
                   labels=c("Amb-Amb", "Amb-Mod", "Amb-Sev", "Elev-Amb", "Elev-Mod", "Elev-Sev"), l=40) + # Use darker colors, lightness=40
  ggtitle("Days 8-14 SMR") +
  geom_point(position=pd, shape=c(21, 21, 21, 24, 24, 24,21, 21, 21, 24, 24, 24,21, 21, 21, 24, 24, 24), size=3, fill=("white")) + # 21 is filled circle
  #expand_limits(y=0) +   
  ylim(3.5, 13) +
  #scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
  theme(legend.justification=c(1,1),
        legend.position=c(1,1))            # Position legend in bottom right
FIGURE.resp_Days.8.14  <- print(FIGURE.resp_Days.8.14 + scale_colour_manual(values = c("skyblue1", "deepskyblue3", "blue", "tomato1", "red1", "firebrick4"))+ theme(legend.position = "none"))



# FIGURE 4
# days 15-21 rel metabolic rate (per individual)    ----------------------------------------------------------------------------------------------- #



resp_Days.15.21.final$Date <- as.character(resp_Days.15.21.final$Date)
FIGURE.resp_Days.15.21 <- ggplot(resp_Days.15.21.final, aes(x=factor(Date), y=resp.COUNT.µg.L.hr.indiv_mean, colour=TREATMENT, group=TREATMENT, shape = Treatment.EXP_2)) + 
  geom_errorbar(aes(ymin=resp.COUNT.µg.L.hr.indiv_mean-resp.COUNT.µg.L.hr.indiv_std.error, 
                    ymax=resp.COUNT.µg.L.hr.indiv_mean+resp.COUNT.µg.L.hr.indiv_std.error), width=.3, position=pd) +
  geom_line(position=pd, width=.3) +
  theme_bw() +
  xlab("Date") +
  ylab("SMR per individual") +
  scale_colour_hue(name="Treatment",    # Legend label, use darker colors
                   breaks=c("AHAA", "AHAM", "AHMA", "AHMM", 
                            "AHSA","AHSM", "EHAA", "EHAM", 
                            "EHMA", "EHMM", "EHSA", "EHSM"),
                   labels=c("Amb-Amb-Amb", "Amb-Amb-Mod", "Amb-Mod-Amb", "Amb-Mod-Mod", 
                            "Amb-Sev-Amb", "Amb-Sev-Mod", "Elev-Amb-Amb","Elev-Amb-Mod",
                            "Elev-Mod-Amb", "Elev-Mod-Mod", "Elev-Sev-Amb", "Elev-Sev-Mod"), l=40) + # Use darker colors, lightness=40
  ggtitle("Days 15-21 SMR") +
  geom_point(position=pd, shape=c(21, 21, 21, 21, 21, 21, 24, 24, 24, 24, 24, 24,  21, 21, 21, 21, 21, 21, 
                                  24, 24, 24, 24, 24, 24, 21, 21, 21, 21, 21, 21,  24, 24, 24, 24, 24, 24), size=3, 
             fill=c("firebrick4", "white", "red1", "white", "tomato1",
                    "white", "blue", "white", "deepskyblue3", "white", "skyblue1",
                    "white", "firebrick4", "white", "red1", "white", "tomato1",
                    "white", "blue", "white", "deepskyblue3", "white", "skyblue1",
                    "white", "firebrick4", "white", "red1", "white", "tomato1", 
                    "white", "blue", "white", "deepskyblue3", "white", "skyblue1","white")) + # all the colors needed
  #expand_limits(y=0) +  
  ylim(3.5, 13) +
  #scale_y_continuous(breaks=0:20*4) +                          # Set tick every 4
  theme(legend.justification=c(1,1),legend.position=c(1,1))    # Position legend in bottom right
FIGURE.resp_Days.15.21  <- print(FIGURE.resp_Days.15.21 + scale_colour_manual(values = c("skyblue1", "skyblue1", "deepskyblue3", "deepskyblue3",
                                                                  "blue", "blue", "tomato1", "tomato1", 
                                                                  "red1", "red1", "firebrick4", "firebrick4")) + theme(legend.position = "none"))

                               

# Arrange plots
SMR.plot <- ggarrange(FIGURE.resp.pre,  FIGURE.resp_Days.1.7, FIGURE.resp_Days.8.14, FIGURE.resp_Days.15.21,ncol = 4, nrow = 1, widths = c(0.7,2,2,2), labels = c("B","C","D","E")) # combine plots 
SMR.plot # view plots
SMR.plot.with.schematic <- ggarrange(FIG.schematic.data, SMR.plot,ncol = 1, nrow = 2, heights = c(1, 2.0), labels = "A") # combine plots
#gsave to Output
ggsave(file="Output/SMR.plot.pdf", SMR.plot.with.schematic, width = 12, height = 8, units = c("in")) # respiration rate plots

























##################################################### #
##################################################### #
########## RELATIVE RESP RATE DATA  ################# #
##################################################### #
##################################################### #

# upload relative SMR data (calculation for this table in "Resp.Tables.R")
df <- read.csv(file="Data/SDR_data/Relative_resp_rates.csv", header=T) #read Size.info data
df.experiment <- df %>% 
  dplyr::filter(Date >20190723) # 21-day experiment began on 20190725
df_Days.1.14 <- df.experiment %>%  # divide dataset into first 14 days (same treatment) 
  dplyr::filter(Date < 20190808)
df_Days.15.21 <- df.experiment %>%  # and last 7 days (same treatment)
  dplyr::filter(Date > 20190807)

#df_Days.1.14$Date <- as.character(df_Days.1.14$Date)
#df_Days.1.14$Treatment <- as.character(df_Days.1.14$Treatment)

########## mean and st error calculations ########################### #
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
  annotate(geom="text", x=2, y=6, label="Reciprocal Amb v. Mod (n = 12 treatments)", color="black") +
  annotate(geom="text", x=1.25, y=4, label="a", color="black") +
  annotate(geom="text", x=1.15, y=-0.2, label="ab", color="black") +
  annotate(geom="text", x=1.10, y=-0.2, label="ab", color="black") +
  annotate(geom="text", x=1, y=-0.2, label="ab", color="black") +
  annotate(geom="text", x=0.9, y=-0.2, label="b", color="black") +
  annotate(geom="text", x=0.82, y=-0.2, label="ab", color="black") +
  annotate(geom="text", x=2.25, y=0.2, label="ab", color="black") +
  annotate(geom="text", x=2.15, y=-0.2, label="ab", color="black") +
  annotate(geom="text", x=2.1, y=1.8, label="ab", color="black") +
  annotate(geom="text", x=2, y=-2.6, label="b", color="black") +
  annotate(geom="text", x=1.9, y=0.3, label="ab", color="black") +
  annotate(geom="text", x=1.82, y=-0.3, label="ab", color="black") +
  annotate(geom="text", x=3.25, y=0.2, label="ab", color="black") +
  annotate(geom="text", x=3.15, y=-1.7, label="b", color="black") +
  annotate(geom="text", x=3.1, y=-0.7, label="ab", color="black") +
  annotate(geom="text", x=3, y=-0.3, label="ab", color="black") +
  annotate(geom="text", x=2.9, y=-0.2, label="ab", color="black") +
  annotate(geom="text", x=2.82, y=0.6, label="ab", color="black") 
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
summary(lm(resp.RAW.µg.L.hr~mean_length, data=Final.resp.table.EXP)) # Adjusted R-squared:  0.1401 


Raw.vs.MEAN.BIOVOLUME <- ggplot(Final.resp.table.EXP, aes(x = mean_biovolume,y = resp.RAW.µg.L.hr)) +
  geom_point() +
  theme_classic() +
  #theme(legend.position = c(0.55,0.96), legend.direction="horizontal", legend.title=element_blank()) +
  #ylim(4,13) + 
  labs(y=expression("Resp.rate.RAW"~(~µg~O[2]*hr^{-1})), 
       x=expression("MEAN.biovolume"~(ml)))
Raw.vs.MEAN.BIOVOLUME <- Raw.vs.MEAN.BIOVOLUME + stat_smooth(method="lm", se=FALSE)
Raw.vs.MEAN.BIOVOLUME
summary(lm(resp.RAW.µg.L.hr~mean_biovolume, data=Final.resp.table.EXP)) # Adjusted R-squared:  0.08455 

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
