#Title: Shell.size 
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

#Required Data files

# Set Working Directory:
# setwd("~/MyProjects/Geoduck_Conditioning/RAnalysis/") #set working
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/")
#Load Sample Info
Size.data <- read.csv(file="Data/Shell_length/20190628_shell_size.csv", header=T) #read sample.info data

# Analysis
ttest_size_20190628 <- t.test(Length ~ Treatment, data = Size.data)
ttest_size_20190628 # view t test results - significant difference between treatments

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
