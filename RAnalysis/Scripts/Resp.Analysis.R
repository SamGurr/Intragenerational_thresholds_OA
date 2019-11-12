#Title: Respiration Analysis
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
library(lsmeans)        # Version: 2.27-62, Date/Publication: 2018-05-11, Depends: methods, R (>= 3.2)

# Set Working Directory:
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/")


############################################################################################# #
############################################################################################# #
#########################  RESP RATE DATA ################################################### #
############################################################################################# #
############################################################################################# #

# upload relative SMR data (calculation for this table in "Resp.Tables.R")
DATA <- read.csv(file="Data/SDR_data/Final_table_for_resp_analysis.csv", header=T) #read Size.info data
DATA.pre <- DATA %>% 
  dplyr::filter(Date %in% 20190723) # 21-day experiment began on 20190723
DATA.experiment <- DATA %>% 
  dplyr::filter(Date > 20190723) # 21-day experiment began on 20190725
DATA_Days.1.14 <- DATA.experiment %>%  # divide dataset into first 14 days (same treatment) 
  dplyr::filter(Date < 20190808)
DATA_Days.15.21 <- DATA.experiment %>%  # and last 7 days (same treatment)
  dplyr::filter(Date > 20190807)

# t.test of "pre" data prior to the experiment
ttest_summ <-t.test(resp.COUNT.µg.L.hr.indiv~Treatment.history, data=DATA.pre) # p-value = 0.5516; no difference between pCO2 treatment
ttest_summ # view tet results

# MODELS FOR DAYS 1 - 7 -------------------------------------------------------- #
Days.1.7 <- DATA_Days.1.14 %>% 
  dplyr::filter(Date < 20190801 ) # filter dataframe for the correct dates
# make date a character to address as a factor in the model
Days.1.7$Date <- as.factor(Days.1.7$Date) 
# interaction plots
interaction.plot(Days.1.7$Treatment.history, Days.1.7$Date, Days.1.7$resp.COUNT.µg.L.hr.indiv)
interaction.plot(Days.1.7$Treatment.history, Days.1.7$Treatment.EXP_1, Days.1.7$resp.COUNT.µg.L.hr.indiv)
interaction.plot(Days.1.7$Date, Days.1.7$Treatment.EXP_1, Days.1.7$resp.COUNT.µg.L.hr.indiv)
# two way ANOVA treatment and date
# resp.MEAN.µg.L.hr.mm
# resp.COUNT.µg.L.hr.indiv
threewayanova_D1.7 <- aov(resp.COUNT.µg.L.hr.indiv ~ Treatment.history*Treatment.EXP_1*Date, data=Days.1.7) # run the model
shapiro.test(residuals(threewayanova_D1.7)) # shaprio wilk test of model residuals p = 0.162; normal distribution
hist((residuals(threewayanova_D1.7)))
summary(threewayanova_D1.7) # marginal effect of treatment history
TukeyHSD(threewayanova_D1.7, 'Treatment.history', conf.level=0.95) # tukey test on the effect of treatment with 95% confidence
TukeyHSD(threewayanova_D1.7, 'Treatment.history:Treatment.EXP_1', conf.level=0.95) # tukey test on the effect of treatment with 95% confidence
# EH and AH to get percent difference; effect of pCO2 history (marginal)
EHandAH.D.1.7 <- Days.1.7 %>% 
  dplyr::group_by(Treatment.history) %>% # group by treatment
  dplyr::summarise(mean.resp = mean(resp.COUNT.µg.L.hr.indiv)) # get the mean values
percent.diff.EHandAH <- ((EHandAH.D.1.7[2,2] - EHandAH.D.1.7[1,2])/ EHandAH.D.1.7[2,2])*100 # 12.372% difference
# EHM and AHM to get percent difference; effect of pCO2 history × treatment (marginal)
Days.1.7$treat.inital <- paste(Days.1.7$Treatment.history,Days.1.7$Treatment.EXP_1, sep ="")
EHMandAHM.D.1.7 <- Days.1.7 %>% dplyr::filter(treat.inital %in% c('EHM', 'AHM'))
EHMandAHM.D.1.7.MEANS <- EHMandAHM.D.1.7 %>% 
  dplyr::group_by(treat.inital) %>% # group by treatment
  dplyr::summarise(mean.resp = mean(resp.COUNT.µg.L.hr.indiv)) # get the mean values
percent.diff.EHMandAHM <- ((EHMandAHM.D.1.7.MEANS[2,2] - EHMandAHM.D.1.7.MEANS[1,2])/ EHMandAHM.D.1.7.MEANS[2,2])*100 # 31.1947% difference
# marginal diff : Treatment.history:Treatment.EXP_1
# EH:M-AH:M p = 0.0408163


# MODELS FOR DAYS 8 - 14 -------------------------------------------------------- #

Days.8.14 <- DATA_Days.1.14 %>% 
  dplyr::filter(Date > 20190731) # make dataframe 
# make date a character to address as a factor in the model
# interaction plots
interaction.plot(Days.8.14$Treatment.history, Days.8.14$Date, Days.8.14$resp.COUNT.µg.L.hr.indiv)
interaction.plot(Days.8.14$Treatment.history, Days.8.14$Treatment.EXP_1, Days.8.14$resp.COUNT.µg.L.hr.indiv)
interaction.plot(Days.8.14$Date, Days.8.14$Treatment.EXP_1, Days.8.14$resp.COUNT.µg.L.hr.indiv)

Days.8.14$Date <- as.factor(Days.8.14$Date) 
# two way ANOVA treatment and date
threewayanova_D8.14 <- aov(resp.COUNT.µg.L.hr.indiv ~ Treatment.history*Treatment.EXP_1*Date, data=Days.8.14) # run the model
shapiro.test(residuals(threewayanova_D8.14)) # shaprio wilk test of model residuals p = 0.8904; normal distribution
hist((residuals(threewayanova_D8.14)))
summary(threewayanova_D8.14) # significant effect of treatment
TukeyHSD(threewayanova_D8.14, 'Date', conf.level=0.95) # tukey test on the effect of treatment with 95% confidence
# significant difference between:
# 20190807-20190801 p = 0.0254093
# D8.14_posthoc <- lsmeans(twowayanova_D8.14, pairwise ~ Date)# pariwise Tukey Post-hoc test between repeated treatments
# D8.14_posthoc.05 <- cld(D8.14_posthoc, alpha=.05, Letters=letters) #letters
# D8.14_posthoc.05


# MODELS FOR DAYS 15 - 21 -------------------------------------------------------- #

# make date a character to address as a factor in the model
DATA_Days.15.21$Date <- as.factor(DATA_Days.15.21$Date) 
#interaction plots
interaction.plot(DATA_Days.15.21$Treatment.history, DATA_Days.15.21$Date, DATA_Days.15.21$resp.COUNT.µg.L.hr.indiv)
interaction.plot(DATA_Days.15.21$Treatment.history, DATA_Days.15.21$Treatment.EXP_1, DATA_Days.15.21$resp.COUNT.µg.L.hr.indiv)
interaction.plot(DATA_Days.15.21$Date, DATA_Days.15.21$Treatment.EXP_1, DATA_Days.15.21$resp.COUNT.µg.L.hr.indiv)

# two way ANOVA treatment and date
fourwayanova_D15.21 <-aov(resp.COUNT.µg.L.hr.indiv ~ Treatment.history*Treatment.EXP_1*Treatment.EXP_2*Date, data=DATA_Days.15.21)
shapiro.test(residuals(fourwayanova_D15.21)) # shaprio wilk test of model residuals p = 0.15
hist((residuals(fourwayanova_D15.21)))
summary(fourwayanova_D15.21) # significant interaction between date and treatment
TukeyHSD(fourwayanova_D15.21, 'Treatment.history:Treatment.EXP_1', conf.level=0.95) # tukey test on the effect of treatment with 95% confidence
# sig effect of : Treatment.history:Treatment.EXP_1
# EH:S-EH:A p = 0.0459861
RESP_Days.15.21.final <- DATA_Days.15.21 %>% 
  dplyr::group_by(Treatment.EXP_1, Treatment.history) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))




############################################################################################# #
############################################################################################# #
######################### RELATIVE RESP RATE DATA ########################################### #
############################################################################################# #
############################################################################################# #

# upload relative SMR data (calculation for this table in "Resp.Tables.R")
df <- read.csv(file="Data/SDR_data/Relative_resp_rates.csv", header=T) #read Size.info data
df.experiment <- df %>% 
  dplyr::filter(Date >20190723) # 21-day experiment began on 20190725
df_Days.1.14 <- df.experiment %>%  # divide dataset into first 14 days (same treatment) 
  dplyr::filter(Date < 20190808)
df_Days.15.21 <- df.experiment %>%  # and last 7 days (same treatment)
  dplyr::filter(Date > 20190807)

######################################################### #
##############   STATISTICAL TESTS ###################### #
######################################################### #

# MODELS FOR DAYS 1 - 7 -------------------------------------------------------- #

df_Days.1.7 <- df_Days.1.14 %>% 
  dplyr::filter(Date < 20190801 ) # filter dataframe for the correct dates
# make a subjects column
df_Days.1.7$subject <- paste(df_Days.1.7$Date, df_Days.1.7$Treatment, sep = "_")
# make date a character to address as a factor in the model
df_Days.1.7$Date <- as.character(df_Days.1.7$Date) 

# two way ANOVA treatment and date
twowayanova_D1.7 <- aov(rel.resp.COUNT ~ Treatment*Date, data=df_Days.1.7) # run the model
shapiro.test(residuals(twowayanova_D1.7)) # shaprio wilk test of model residuals p = 0.162; normal distribution
hist((residuals(twowayanova_D1.7)))
summary(twowayanova_D1.7) # significant effect of treatment
#posthoc
TukeyHSD(twowayanova_D1.7, 'Treatment', conf.level=0.95) # tukey test on the effect of treatment with 95% confidence
# significant difference between:
# EHM-AHM
# EHS-EHM
D1.7_posthoc <- lsmeans(twowayanova_D1.7, pairwise ~  Treatment)# pariwise Tukey Post-hoc test between repeated treatments
D1.7_posthoc.05 <- cld(D1.7_posthoc, alpha=.05, Letters=letters) #letters

# linear mixed effect model with random effect of subject
model_lmer.FULL.D1.14 = lmer(rel.resp.COUNT ~ Treatment  +
                               (1|subject), data=df_Days.1.14, REML=FALSE)
model_lmer.NULL.D1.14 = lmer(rel.resp.COUNT ~ 1 +
                               (1|subject), data=df_Days.1.14, REML=FALSE)
anova(model_lmer.FULL.D1.14, model_lmer.NULL.D1.14)
coef(model_lmer.FULL.D1.14)


# MODELS FOR DAYS 8 - 14 -------------------------------------------------------- #

df_Days.8.14 <- df_Days.1.14 %>% 
  dplyr::filter(Date > 20190731) # make dataframe 
# make a subjects column
df_Days.8.14$subject <- paste(df_Days.8.14$Date, df_Days.8.14$Treatment, sep = "_")
# make date a character to address as a factor in the model
df_Days.8.14$Date <- as.character(df_Days.8.14$Date) 
# two way ANOVA treatment and date
twowayanova_D8.14 <- aov(rel.resp.COUNT ~ Treatment*Date, data=df_Days.8.14) # run the model
shapiro.test(residuals(twowayanova_D8.14)) # shaprio wilk test of model residuals p = 0.8904; normal distribution
hist((residuals(twowayanova_D8.14)))
summary(twowayanova_D8.14) # significant effect of treatment
TukeyHSD(twowayanova_D8.14, 'Treatment:Date', conf.level=0.95) # tukey test on the effect of treatment with 95% confidence
# significant difference between:
# EHM-AHM
D8.14_posthoc <- lsmeans(twowayanova_D8.14, pairwise ~  Treatment:Date)# pariwise Tukey Post-hoc test between repeated treatments
D8.14_posthoc.05 <- cld(D8.14_posthoc, alpha=.05, Letters=letters) #letters

# MODELS FOR DAYS 15 - 21 -------------------------------------------------------- #

# make a subjects column
df_Days.15.21$subject <- paste(df_Days.15.21$Date, df_Days.15.21$Treatment, sep = "_")
# make date a character to address as a factor in the model
df_Days.15.21$Date <- as.character(df_Days.15.21$Date) 
# two way ANOVA treatment and date
twowayanova_D15.21 <-aov(rel.resp.MEAN.LENGTH ~ Treatment*Date, data=df_Days.15.21)
shapiro.test(residuals(twowayanova_D15.21)) # shaprio wilk test of model residuals p = 0.15
hist((residuals(twowayanova_D15.21)))
summary(twowayanova_D15.21) # significant interaction between date and treatment
TukeyHSD(twowayanova_D15.21, 'Treatment', conf.level=0.95) # tukey test on the effect of treatment with 95% confidence
# moderate difference (p<0.1) between:
# EHSM-AHSM 
TukeyHSD(twowayanova_D15.21, 'Treatment:Date', conf.level=0.95) # tukey test on the effect of treatment:date with 95% confidence
# significant difference between:
# EHMM:20190814-EHSM:20190808 - 0.0388647
# AHSM:20190811-EHSM:20190808 - 0.0027614
# EHSM:20190808-AHMM:20190808 - 0.0500898
D15.21_posthoc <- lsmeans(twowayanova_D15.21, pairwise ~  Treatment:Date)# pariwise Tukey Post-hoc test between repeated treatments
D15.21_posthoc.05 <- cld(D15.21_posthoc, alpha=.05, Letters=letters) #letters

# linear mixed effect model with random effect of subject
model_lmer.FULL.D15.21 = lmer(rel.resp.COUNT ~ Treatment +
                                (1|subject), data=df_Days.15.21, REML=FALSE)

model_lmer.NULL.D15.21 = lmer(rel.resp.COUNT ~ 1 +
                                (1|subject), data=df_Days.15.21, REML=FALSE)

anova(model_lmer.FULL.D15.21, model_lmer.NULL.D15.21)

model <- aov(lm(rel.resp.COUNT ~ Treatment*Date, data=df_Days.15.21))
summary(model)
shapiro.test(residuals(model))
TukeyHSD(model)

model2 <-aov(lm(rel.resp.COUNT ~ Treatment*Date, data=df_Days.1.14))
summary(model2)
shapiro.test(residuals(model2))
TukeyHSD(model2)