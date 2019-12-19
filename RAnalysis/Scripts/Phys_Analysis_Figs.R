#Author: Sam Gurr 
#Edited by: Sam Gurr
#Date Last Modified: 20191107
#Purpose: Calls the Table (From Phys_Calc.R, output as "Phys.Assay.Table.csv")
#See Readme file for details

rm(list=ls())

# Install packages if not already in your library-----------------------------------------------------------------------------------------------
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("lsmeans" %in% rownames(installed.packages()) == 'FALSE') install.packages('lsmeans') 
if ("ggpubr" %in% rownames(installed.packages()) == 'FALSE') install.packages('ggpubr') 
# Load packages and pacage version/date/import/depends info
library(dplyr) # Version 0.7.6, Packaged: 2018-06-27, Depends: R (>= 3.1.2)Imports: assertthat (>= 0.2.0), bindrcpp (>= 0.2.0.9000), glue (>=1.1.1), magrittr (>= 1.5), methods, pkgconfig (>= 2.0.1), R6(>= 2.2.2), Rcpp (>= 0.12.15), rlang (>= 0.2.0), tibble (>=1.3.1), tidyselect (>= 0.2.3), utils
library(lsmeans) # Version: 2.27-62, Date/Publication: 2018-05-11, Depends: methods, R (>= 3.2)
library(ggpubr)
#set working directory--------------------------------------------------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis") #set working

# NOTE: This script analyzes a subsample of the significant effect(s) found for Resp rate (review /Scripts/Resp.Analysis.R) 
# and analyzes the first run of phys assays of 84 chosen samples TARGETTED to these significant effects of treatment on respiration rate
# Resp rate effects on....
# Days 1 - 7: Effect of initial exposure × history - EHM and AHM -  resp rate marginal difference with greater rates
# Days 15 - 21: Effect of initial exposure × history;  EHS and EHA - resp rate significant difference with greater rates in EHS
# Review the following 

#  ASSEMBLE ALL DATA FOR GRAPHS AND ANALYSIS ------------------------------------- #
# Total antioxidant capacity, total protein and ARDW data
PHYS.DATA <- read.csv(file="Output/Phys.Assay.Table.csv", header=T) #read PPhys assay data output from Phys_Calc.R (TAC, TP, AFDW of targetted samples)
# Respiration rate data
resp.data <- read.csv(file="Data/SDR_data/Final_table_for_resp_analysis.csv", header=T) #read Size.info resp.data
resp.data.Days.1.7 <- resp.data %>% dplyr::filter(Date %in% 20190725:20190731) %>% dplyr::select(Date, resp.COUNT.µg.L.hr.indiv,resp.MEAN.µg.L.hr.mm, Treatment.EXP_1, Treatment.history) # d 1 - 7 experiment
resp.data_Days.8.14 <- resp.data %>% dplyr::filter(Date %in% 20190801:20190807) %>% dplyr::select(Date, resp.COUNT.µg.L.hr.indiv,resp.MEAN.µg.L.hr.mm, Treatment.EXP_1, Treatment.history) # d 8 - 14 experiment
resp.data_Days.15.21 <- resp.data %>% dplyr::filter(Date %in% 20190808:20190814) %>% dplyr::select(Date, resp.COUNT.µg.L.hr.indiv,resp.MEAN.µg.L.hr.mm, Treatment.EXP_2, Treatment.history, Treatment.EXP_1) # d 15-21 experiment

################################################################ #
################################################################ #
################################################################ #
######     DAY 1 -7 initial exposure × history       ###########
##########       RESP RATES EHM > AHM  ######################### #
################################################################ #
################################################################ #

#initial effect
Model_1_RESP.rate <- aov(resp.MEAN.µg.L.hr.mm~Treatment.history*Treatment.EXP_1*Date, data = resp.data.Days.1.7) # run model 
summary(Model_1_RESP.rate)
shapiro.test(residuals(Model_1_RESP.rate)) # p-value = 0.2751; normal via shapiro wilk test
hist((residuals(Model_1_RESP.rate))) # histogram of residuals
qqnorm(residuals(Model_1_RESP.rate)) # qqplot
# post-hoc 
Model_1_RESP.rate.MEANS <- lsmeans(Model_1_RESP.rate, pairwise ~  Treatment.history:Treatment.EXP_1)# pariwise Tukey Post-hoc test between repeated treatments
# AH,M - EH,M -1.07408508 0.3831895 98  -2.803  0.0654 MARGINAL DIFF BETWEEN AHM AND EHM
resp.data.Days.1.7$Treatment.initial <- paste(resp.data.Days.1.7$Treatment.history, resp.data.Days.1.7$Treatment.EXP_1, sep="")
resp.data.Days.1.7.EFFECT <- resp.data.Days.1.7 %>%  dplyr::filter(Treatment.initial %in% c("EHM", "AHM"))
plot_initial.resp <- ggviolin(resp.data.Days.1.7.EFFECT, ylab =expression("Respiration rate"~(~µg~O[2]*hr^{-1}*indiv^{-1})),  
                              x = "Treatment.initial", y = "resp.COUNT.µg.L.hr.indiv",  fill = "Treatment.initial",  
                              palette = c("#FC4E07", "#00AFBB"), add = "none", title = "Respiration rate")
plot_initial.resp

# stat.test.RESP<- compare_means(
# resp.MEAN.µg.L.hr.mm~Treatment.initial, data = resp.data.Days.1.7.EFFECT,
# method = "t.test")
# stat.test.RESP <- stat.test.RESP %>%# Add manually p-values from stat.test data
# Insert the results of the tukey test in the Resp.Analysis.R Script
# result from three way anova history × exposure = p = 0.0647; tukey test EH:M > AH:M p = 0.0408
threewayanova.resp.table <- data.frame(matrix(nrow = 1, ncol = 4))
colnames(threewayanova.resp.table)<-c('test', 'group1', 'group2', 'p-val_tukey') # column names, matches all data
threewayanova.resp.table[1,] <- c("initial", "EHM", "AHM", 0.0408)
threewayanova.resp.table <- threewayanova.resp.table %>%# Add manually p-values from stat.test data
  mutate(y.position = c(22))
plot_initial.resp <- plot_initial.resp %>% ggadd(c("boxplot", "jitter"), fill = "white")  + 
  stat_pvalue_manual(threewayanova.resp.table, label ="p-val_tukey") # Add box plot
plot_initial.resp <- ggpar(plot_initial.resp, ylim = c(0,35))
plot_initial.resp

##################################################################################################### #
##################################################################################################### #
# INITIAL EFFECT - Day 7 ############################################################################ #
##################################################################################################### #
##################################################################################################### #

# data prep
Data.D.7 <- PHYS.DATA %>% dplyr::filter(Date.fixed == '20190731')  # call the data
Data.D.7$Treat_history <- substr(Data.D.7$Treatment, 1,1) # Treat_history
Data.D.7$secondary <- substr(Data.D.7$Treatment, 3,3) # initial
Data.D.7 <- Data.D.7 %>% filter(µM.CRE.mg.protein > 0) # ommit values below 0
Data.D.7$secondary <- factor(Data.D.7$secondary, levels = c("A", "M", "S")) # arrange levels allows ggpubr to plot correctly

# TAOC ##################################################################################################### #
# run the model on raw data
model_Day7_TAC <- aov(µM.CRE.mg.protein ~ Treat_history * secondary, data = Data.D.7)
summary(model_Day7_TAC) # no sig difference
shapiro.test(residuals(model_Day7_TAC)) # p-value =  5.55e-07; non normal via shapiro wilk test - likely due to single outlier
hist((residuals(model_Day7_TAC))) # histogram of residuals - left skew appears to have an outlier
qqnorm(residuals(model_Day7_TAC)) # qqplot
# find the outlier and ommit to run model with the ommission
outliers <- boxplot(Data.D.7$µM.CRE.mg.protein , plot=FALSE)$out # id the outlier
print(outliers) # view the outlier
Data.D.7.OM <- Data.D.7[-which(Data.D.7$µM.CRE.mg.protein %in% outliers),] # remove outlier from new data frame
# run the model with ommitted outlier
model_Day7_TAC.OM<- aov(µM.CRE.mg.protein ~ Treat_history * initial, data = Data.D.7.OM) # run the model without the outlier
summary(model_Day7_TAC.OM) # view model - no significant diffs
shapiro.test(residuals(model_Day7_TAC.OM)) # p-value = 0.8191; normal via shapiro wilk test
hist((residuals(model_Day7_TAC.OM))) # histogram of residuals
qqnorm(residuals(model_Day7_TAC.OM)) # qqplot
# try to transform raw data LOG without removing the outlier
model_Day7_TAC.LOG <- aov((log(µM.CRE.mg.protein)) ~ Treat_history * secondary, data = Data.D.7)
summary(model_Day7_TAC.LOG) # no sig difference
shapiro.test(residuals(model_Day7_TAC.LOG)) # p-value = 0.01333; non normal via shapiro wilk test, tranformation does NOT resolve 
hist((residuals(model_Day7_TAC.LOG))) # histogram of residuals - left skew appears to have an outlier
qqnorm(residuals(model_Day7_TAC)) # qqplot
# plot the data
D.7.PLOTbox<- ggboxplot(Data.D.7, x = "secondary", y = "µM.CRE.mg.protein",  ylab = "µM.CRE.mg.protein",  fill = "Treat_history",
                             palette = c( "#00AFBB", "#FC4E07"),add = "jitter", title = "Secondary TAOC Day 7", xlab = "Secondary pCO2 Treatment")
D.7.PLOTbox <- ggpar(D.7.PLOTbox, ylim = c(0,200))
D.7.PLOTbox

# Protein ##################################################################################################### #
# run the model on raw data
model_Day7_AFDW <- aov(mgProtein.mgAFDW  ~ Treat_history * secondary, data = Data.D.7)
summary(model_Day7_AFDW) # no significant difference
shapiro.test(residuals(model_Day7_AFDW)) # p-value = 0.5455; normal via shapiro wilk test
hist((residuals(model_Day7_AFDW))) # histogram of residuals
qqnorm(residuals(model_Day7_AFDW)) # qqplot
D.7.PLOTbox.Protein<- ggboxplot(Data.D.7, x = "secondary", y = "mgProtein.mgAFDW",  ylab = "mgProtein.mgAFDW",  fill = "Treat_history",
                        palette = c( "#00AFBB", "#FC4E07"),add = "jitter", title = "Secondary TP Day 7",xlab = "Secondary pCO2 Treatment")
D.7.PLOTbox.Protein <- ggpar(D.7.PLOTbox.Protein, ylim = c(0,35))
D.7.PLOTbox.Protein

##################################################################################################### #
#DAY 21 TEST ######################################################################################## #
##################################################################################################### #
##################################################################################################### #
Data.D.21 <- PHYS.DATA %>% dplyr::filter(Date.fixed == '20190814') # call the data
Data.D.21$Treat_history <- substr(Data.D.21$Treatment, 1,1) # Treat_history
Data.D.21$secondary <- substr(Data.D.21$Treatment, 3,3) # secondary
Data.D.21$tertiary <- substr(Data.D.21$Treatment, 4,4) # tertiary
Data.D.21$secondary <- factor(Data.D.21$secondary , levels = c("A", "M", "S")) # arrange levels allows ggpubr to plot correctly

# TAOC ##################################################################################################### #
# data prep
Data.D.21 <- Data.D.21 %>% filter(µM.CRE.mg.protein > 0) # ommit values below 0
Data.D.21$µM.CRE.mg.protein <- ifelse(Data.D.21$µM.CRE.mg.protein  < 0.0000, 0, Data.D.21$µM.CRE.mg.protein ) # values < 0 make 0
# diagnostic
count.treat <- Data.D.21 %>%  dplyr::group_by(Treatment) %>% # diagnostic of coutns for TAOC - still have more samples to run
  dplyr::summarise(count_replicated.by.day =n())
count.treat # STILL NEED TO ANALYZE SOME SAMPLES
# mdoel with raw data
model_Day21_TAC <- aov(µM.CRE.mg.protein ~ Treat_history * secondary * tertiary, data = Data.D.21) # run model for TAOC raw data
summary(model_Day21_TAC) # marginal difference due to treatment history
shapiro.test(residuals(model_Day21_TAC)) # p-value = 0.001127; non normal via shapiro wilk test
hist((residuals(model_Day21_TAC))) # histogram of residuals - LEFT SKEW - LOG TRANSFORM
qqnorm(residuals(model_Day21_TAC)) # qqplot
# odel with transfomed data (LOG) - sig diff from treatment history
model_Day21_TAC <- aov((log(µM.CRE.mg.protein)) ~ Treat_history * secondary * tertiary, data = Data.D.21) # run model for Log transformed to resolce left skew
summary(model_Day21_TAC) # transformed data shows a significant difference due to treatment history
TukeyHSD(model_Day21_TAC) # E-A -0.1991035 -0.3751744 -0.02303262 0.0278336
shapiro.test(residuals(model_Day21_TAC)) # p-value = 0.7787; normal via shapiro wilk test resolved left skew and normality assumptions
hist((residuals(model_Day21_TAC))) # histogram of residuals - GOOD
qqnorm(residuals(model_Day21_TAC)) # qqplot - GOOD
# plot the data
plot.DAY21.TAC <- ggboxplot(Data.D.21, x = "Treat_history", y = "µM.CRE.mg.protein",  ylab = "µM.CRE.mg.protein", fill = "Treat_history",
                               palette = c("#00AFBB", "#FC4E07"), add = "jitter",  title ="Tertiary TAOC Day 21", xlab = "Initial pCO2 Treatment")
plot.DAY21.TAC <- ggpar(plot.DAY21.TAC, ylim = c(0,200)) + geom_bracket(xmin = "A", xmax = "E", y.position = 220, label = "*", tip.length = 0.01)
D.21.PLOTbox<- ggboxplot(Data.D.21, x = "secondary", y = "µM.CRE.mg.protein",  ylab = "µM.CRE.mg.protein",  fill = "Treat_history",
                         palette = c( "#00AFBB", "#FC4E07"),add = "jitter", shape = "tertiary",
                        xlab = "Secondary pCO2 Treatment", title = "Tertiary TAOC Day 21")
D.21.PLOTbox <- ggpar(D.21.PLOTbox, ylim = c(0,200)) 
PLOTS.TAC.D21 <- ggarrange(plot.DAY21.TAC,D.21.PLOTbox, nrow = 1, widths = c(0.5, 1))
PLOTS.TAC.D21 # view arranged plots

# Protein ##################################################################################################### #
# run the model on raw data
model_Day21_Protein <- aov(mgProtein.mgAFDW  ~ Treat_history * secondary * tertiary, data = Data.D.21)
summary(model_Day21_Protein)
shapiro.test(residuals(model_Day21_Protein)) # p-value = 1e-06; non normal via shapiro wilk test
hist((residuals(model_Day21_Protein))) # histogram of residuals left skew
qqnorm(residuals(model_Day21_Protein)) # qqplot
# look for outliers (seems like there are in qqnorm plot)
outliers.d21_protein <- boxplot(Data.D.21$mgProtein.mgAFDW , plot=FALSE)$out # id the outlier
print(outliers.d21_protein) # view the outlier
# run model with LOG transformation
model_Day21_Protein.LOG <- aov((log(mgProtein.mgAFDW))  ~ Treat_history * secondary * tertiary, data = Data.D.21)
summary(model_Day21_Protein.LOG)
shapiro.test(residuals(model_Day21_Protein.LOG)) # p-value = 0.04893; non normal via shapiro wilk test
hist((residuals(model_Day21_Protein.LOG))) # histogram of residuals left skew
qqnorm(residuals(model_Day21_Protein.LOG)) # qqplot
# plot the data
D.21.PLOTbox.Protein<- ggboxplot(Data.D.21, x = "secondary", y = "mgProtein.mgAFDW",  ylab = "mgProtein.mgAFDW",  fill = "Treat_history",
                                palette = c( "#00AFBB", "#FC4E07"),add = "jitter",shape = "tertiary", title = "Tertiary TP Day 21",
                                xlab = "Secondary pCO2 Treatment")
D.21.PLOTbox.Protein <- ggpar(D.21.PLOTbox.Protein, ylim = c(0,35))
D.21.PLOTbox.Protein


TAC.TP.D7.21_PLOTS <- ggarrange(D.7.PLOTbox, D.7.PLOTbox.Protein, PLOTS.TAC.D21, D.21.PLOTbox.Protein, nrow = 2, ncol = 2, widths = c(1, 1, 0.5, 1), labels = c("A","B", "C", "D"))
ggsave(file="Output/TAOC_plots_d7.and.d21.pdf", TAC.TP.D7.21_PLOTS, width = 12, height = 8, units = c("in")) # respiration rate plots







# INITIAL EFFECT - Marginal diff in metabolic response EHM > AHM during initial subseq exposure (Days 1 4 7)
# prep the data
Data.D.1.4.7 <- PHYS.DATA %>% dplyr::filter(Date.fixed < 20190801) # call data by DATE for the correct timeline
Data.D.1.4.7$Treatment <- factor(Data.D.1.4.7$Treatment, levels = c("EHM","AHM"))
count <- Data.D.1.4.7 %>%  dplyr::group_by(Treatment) %>% dplyr::summarise(count = n()) # 18 (AHM) and 17 (EHM) look at the count of replicates per day (removed NAs due to testing, under pipetting, and spec errors)
# ASH FREE DRY WEIGHT
Model_1_AFDW <- aov(mgTOTAL_AFDW~Treatment*Exprmt.Day, data = Data.D.1.4.7) # run model 
summary(Model_1_AFDW)
shapiro.test(residuals(Model_1_AFDW)) # p-value = 0.0003713; non normal via shapiro wilk test
hist((residuals(Model_1_AFDW))) # histogram of residuals
qqnorm(residuals(Model_1_AFDW)) # qqplot
# transform to resolve normality assumptions
Data.D.1.4.7$mgASFW.LOG.TRANS <- log(Data.D.1.4.7$mgTOTAL_AFDW)
Model_1_AFDW.trans <- aov(mgASFW.LOG.TRANS~Treatment*Exprmt.Day, data = Data.D.1.4.7) # run model 
summary(Model_1_AFDW.trans)
shapiro.test(residuals(Model_1_AFDW.trans)) # p-value =  0.2881; non normal via shapiro wilk test
hist((residuals(Model_1_AFDW.trans))) # histogram of residuals
qqnorm(residuals(Model_1_AFDW.trans)) # qqplot
# post-hoc on transformed data - sig effect of treatment
Model_1_AFDW.trans.MEANS <- lsmeans(Model_1_AFDW.trans, pairwise ~  Treatment)# pariwise Tukey Post-hoc test between repeated treatments
#plots
plot_intitial.AFDW <- ggviolin(Data.D.1.4.7, x = "Treatment", y = "mgTOTAL_AFDW",  ylab = "mg AFDW", fill = "Treatment",
                               palette = c("#FC4E07", "#00AFBB"), add = "none",  title ="mgAFDW")
plot_intitial.AFDW %>% ggadd(c("boxplot", "jitter"), color = "Treatment") 
stat.test.AFDW <- compare_means(
  mgASFW.LOG.TRANS~Treatment, data = Data.D.1.4.7,
  method = "t.test")
stat.test.AFDW <- stat.test.AFDW %>%# Add manually p-values from stat.test data
  mutate(y.position = c(19))
plot_intitial.AFDW <- plot_intitial.AFDW %>% ggadd(c("boxplot", "jitter"), fill = "white")  + stat_pvalue_manual(stat.test.AFDW, label = "p.adj") # Add box plot
plot_intitial.AFDW <- ggpar(plot_intitial.AFDW, ylim =c(-2,20))
plot_intitial.AFDW

# TOTAL PROTEIN (corrected for AFDW)
Model_1_TP <- aov(mgProtein.mgAFDW~Treatment*Exprmt.Day, data = Data.D.1.4.7) # run model 
summary(Model_1_TP) # experiment and day 0.0470 *
shapiro.test(residuals(Model_1_TP)) # p-value = 0.7212; normal via shapiro wilk test
hist((residuals(Model_1_TP))) # histogram of residuals
qqnorm(residuals(Model_1_TP)) # qqplot
TukeyHSD(Model_1_TP) # AHM:Day4-EHM:Day1
#plots
plot_intitial.TP <- ggviolin(Data.D.1.4.7, x = "Treatment", y = "mgProtein.mgAFDW", fill = "Treatment",
                             palette = c("#FC4E07", "#00AFBB"), ylab =expression("Total Protein"~(~mg_protein_mg_AFDW^{-1})),
                             add = "none", title ="mg.Total Protein")
plot_intitial.TP <- plot_intitial.TP %>% ggadd(c("boxplot", "jitter"), fill = "white")
plot_intitial.TP <- ggpar(plot_intitial.TP, ylim = c(0,40))
plot_intitial.TP

# TOTAL ANTIOXIDANT CAPACIITY
Model_2_TAC <- aov(µM.CRE.mg.protein ~Treatment*Exprmt.Day, data = Data.D.1.4.7) # run model 
summary(Model_2_TAC) # no effect
# test residuals for normality assumptions
shapiro.test(residuals(Model_2_TAC)) # p-value = 0.6034; normal via shapiro wilk test
hist((residuals(Model_2_TAC))) # histogram of residuals
qqnorm(residuals(Model_2_TAC)) # qqplot
# REMOVE OUTLIERS MODEL
#outliers <- boxplot(Data.D.1.4.7$µM.CRE.µg.protein  , plot=FALSE)$out # call outliers
#Data.D.1.4.7[which(Data.D.1.4.7$µM.CRE.µg.protein  %in% outliers),] # call ouliers
#Data.D.1.4.7.OM <- Data.D.1.4.7[-which(Data.D.1.4.7$µM.CRE.µg.protein  %in% outliers),] # ommit ouliers new dataframe ouliers
#Model_1_TAC.OM <- aov(µM.CRE.µg.protein ~Treatment*Exprmt.Day, data = Data.D.1.4.7.OM) # run model WITH NEGATIVE VALUES AS ZERO
#summary(Model_1_TAC.OM) # P = 0.0338 Treatment
#shapiro.test(residuals(Model_1_TAC.OM)) # p-value = 0.1127; normal via shapiro wilk test - removal of ouliers resolved
#hist((residuals(Model_1_TAC.OM))) # histogram of residuals
#qqnorm(residuals(Model_1_TAC.OM)) # qqplot
#Data.D.1.4.7.OM.MEANS <- lsmeans(Model_1_TAC.OM, pairwise ~  Treatment)# pariwise Tukey Post-hoc test between repeated treatments
#Data.D.1.4.7.OM.MEANS$lsmeans
# Figure
plot_intitial.TAC <- ggviolin(Data.D.1.4.7, x = "Treatment", y = "µM.CRE.mg.protein", fill = "Treatment",
                              palette = c("#FC4E07", "#00AFBB"), ylab =expression("µM Cu Equivalents "~~mgProtein^{-1}), add = "none", title = "Total Antioxidant Capacity")
plot_intitial.TAC %>% ggadd(c("boxplot", "jitter"), color = "Treatment") 
stat.test.TAC <- compare_means(
  µM.CRE.mg.protein ~Treatment, data = Data.D.1.4.7,
  method = "t.test")
stat.test.TAC <- stat.test.TAC %>%# Add manually p-values from stat.test data
  mutate(y.position = c(50000))
plot_intitial.TAC <- plot_intitial.TAC %>% ggadd(c("boxplot", "jitter"), fill = "white") + stat_pvalue_manual(stat.test.TAC, label = "p.adj") # Add box plot
plot_intitial.TAC <- ggpar(plot_intitial.TAC, ylim = c(-0,1200))
plot_intitial.TAC

# HIGHER ANTIOXIDANTS EXPRESSED IN THE ELEVATED HISTORY ANIMALS
ALL.initial.PLOTS <- ggarrange(plot_initial.resp, plot_intitial.AFDW, plot_intitial.TP, plot_intitial.TAC, ncol = 4, nrow = 1)
ALL.initial.PLOTS

################################################################ #
################################################################ #
################################################################ #
######     DAY 15-21 subseq recip effect; exposure × history  ####
##########       RESP RATES EHS > EHA  ######################### #
################################################################ #
################################################################ #

# SUBSEQUENT RECIPROCAL EFFECT
#secondary subseq effect
resp.data_Days.15.21$Treatment <- paste(resp.data_Days.15.21$Treatment.initial, resp.data_Days.15.21$Treatment.EXP_2, sep = "")
resp.data_Days.15.21$Date <- as.factor(resp.data_Days.15.21$Date)
Model_2_RESP.rate <- aov(resp.MEAN.µg.L.hr.mm~Treatment.EXP_2*Treatment.history*Treatment.EXP_1*Date, data = resp.data_Days.15.21) # run model 
summary(Model_2_RESP.rate) # p = 0.0122 Treatment.history:Treatment.EXP_1
shapiro.test(residuals(Model_2_RESP.rate)) # p-value = 0.4021; normal via shapiro wilk test
hist((residuals(Model_2_RESP.rate))) # histogram of residuals
qqnorm(residuals(Model_2_RESP.rate)) # qqplot
# post-hoc
TukeyHSD(Model_2_RESP.rate)
Model_2_resp.POSTHOC <- lsmeans(Model_2_RESP.rate, pairwise ~  Treatment.history:Treatment.EXP_1)# pariwise Tukey Post-hoc test between repeated treatments
Model_2_resp.POSTHOC # view post hoc summary 0.0830 EH,A - EH,S
resp.data_Days.15.21$Treatment.initial <- paste(resp.data_Days.15.21$Treatment.history, resp.data_Days.15.21$Treatment.EXP_1, sep="")
resp.data_Days.15.21.EFFECT <- resp.data_Days.15.21 %>%  dplyr::filter(Treatment.initial %in% c("EHS", "EHA"))
plot_second.subes.resp <- ggviolin(resp.data_Days.15.21.EFFECT, x = "Treatment.initial", y = "resp.MEAN.µg.L.hr.mm", 
                                   ,ylab =expression("Respiration rate"~(~µg~O[2]*hr^{-1}*indiv^{-1})), fill = "Treatment.initial",
                                   palette = c("#FC4E07", "#00AFBB"), add = "none", title = "Respiration rate")
fourwayanova.resp.table <- data.frame(matrix(nrow = 1, ncol = 4))
colnames(fourwayanova.resp.table)<-c('test', 'group1', 'group2', 'p-val_tukey') # column names, matches all data
fourwayanova.resp.table[1,] <- c("initial", "EHS", "EHA", "aov 0.0122 / post-hoc 0.0830")
fourwayanova.resp.table <- fourwayanova.resp.table %>%# Add manually p-values from stat.test data
  mutate(y.position = c(9))
plot_second.subes.resp <- plot_second.subes.resp %>% ggadd(c("boxplot", "jitter"), shape ="Treatment",  fill = "white")  +
  stat_pvalue_manual(fourwayanova.resp.table, label ="p-val_tukey") # Add box plot
plot_second.subes.resp <- ggpar(plot_second.subes.resp, ylim = c(0,10))
plot_second.subes.resp

# prep the PHYS ASSAY DATA
Data.D.15.18.21 <- PHYS.DATA %>% filter(Date.fixed > 20190807) 
count <- Data.D.15.18.21 %>%  dplyr::group_by(Treatment) %>% dplyr::summarise(count = n())
Data.D.15.18.21$Treatment.initial <- substr(Data.D.15.18.21$Treatment, 1,3)
Data.D.15.18.21$Treatment.initial <- factor(Data.D.15.18.21$Treatment.initial, levels = c("EHS", "EHA"))
Data.D.15.18.21$Treatment.subseq <- substr(Data.D.15.18.21$Treatment, 4,4)
# ASH FREE DRY WEIGHT
Model_2_AFDW <- aov(mgTOTAL_AFDW~Treatment.initial*Treatment.subseq*Exprmt.Day, data = Data.D.15.18.21)
summary(Model_2_AFDW) # Treatment.initial:Exprmt.Day  0.0382 *
# test residuals for normality assumptions
shapiro.test(residuals(Model_2_AFDW)) # p-value = 0.09782; NORMAL via shapiro wilk test
hist((residuals(Model_2_AFDW))) # histogram of residuals; left skew
qqnorm(residuals(Model_2_AFDW)) # qqplot
# post-hoc
TukeyHSD(Model_2_AFDW)
Model_2_AFDW.tr.POSTHOC <- lsmeans(Model_2_AFDW, pairwise ~  Treatment.initial:Exprmt.Day)# pariwise Tukey Post-hoc test between repeated treatments
Model_2_AFDW.tr.POSTHOC # view post hoc summary 
# EHS,DAY21 - EHA,DAY21  7.12424421 2.559908 33   2.783  0.0857
# EHS,DAY18 - EHS,DAY21 -8.04663149 2.656542 33  -3.029  0.0494
# EHA,DAY15 - EHS,DAY21 -7.03363265 2.559908 33  -2.748  0.0924
# EHS,DAY15 - EHS,DAY21 -7.60591851 2.656542 33  -2.863  0.0719
# plot the data
plot_subsequent.AFDW <- ggviolin(Data.D.15.18.21, x = "Treatment.initial", y = "mgTOTAL_AFDW",  ylab = "mg AFDW",  fill = "Treatment.initial",
                                 palette = c("#FC4E07", "#00AFBB"),add = "none", title = "mgAFDW")
plot_subsequent.AFDW <- plot_subsequent.AFDW %>% ggadd(c("boxplot", "jitter"),shape ="Treatment", fill = "white") # Add box plot
plot_subsequent.AFDW <- ggpar(plot_subsequent.AFDW, ylim = c(-2,35))
plot_subsequent.AFDW

# TOTAL PROTEIN (corrected for AFDW)
Model_2_TP <- aov(mgProtein.mgAFDW~Treatment.initial*Treatment.subseq*Exprmt.Day, data = Data.D.15.18.21)
summary(Model_2_TP) # Treatment.initial:Exprmt.Day  
# test residuals for normality assumptions
shapiro.test(residuals(Model_2_TP)) # p-value = 0.1411; NORMAL via shapiro wilk test
hist((residuals(Model_2_TP))) # histogram of residuals; left skew
qqnorm(residuals(Model_2_TP)) # qqplot
# plot the data
plot_subsequent.TP <- ggviolin(Data.D.15.18.21, x = "Treatment.initial", y = "mgProtein.mgAFDW", fill = "Treatment.initial",
                               ylab =expression("Total Protein"~(~mg_protein_mg_AFDW^{-1})), palette = c("#FC4E07", "#00AFBB"),add = "none",  title ="mgProtein.mgAFDW")
plot_subsequent.TP <- plot_subsequent.TP %>% ggadd(c("boxplot", "jitter"), shape ="Treatment", fill = "white") # Add box plot
plot_subsequent.TP <- ggpar(plot_subsequent.TP, ylim = c(0,40))
plot_subsequent.TP

# TOTAL ANTIOXIDANT CAPACIITY
# Run the model (3-way anova Treatment*Treatment*time)
Model_2_TAC <- aov(µM.CRE.mg.protein ~Treatment.initial*Treatment.subseq*Exprmt.Day, data = Data.D.15.18.21)
summary(Model_2_TAC) # p 0.0112 *Treatment.initial:Exprmt.Day          
# test residuals for normality assumptions
shapiro.test(residuals(Model_2_TAC)) # p-value = 0.4925; NON NORMAL
hist((residuals(Model_2_TAC))) # histogram of residuals; left skew
qqnorm(residuals(Model_2_TAC)) # qqplot
# POSTHOC
TukeyHSD(Model_2_TAC) # EHA:M:DAY18-EHS:A:DAY18 0.0657267
# Model_2_TAC.tr.POSTHOC <- lsmeans(Model_2_TAC, pairwise ~  Treatment.initial:Exprmt.Day)# pariwise Tukey Post-hoc test between repeated treatments
# Model_2_TAC.tr.POSTHOC #  EHS,DAY18 - EHS,DAY21 0.0753


# plot the data
plot_subsequent.TAC <- ggviolin(Data.D.15.18.21, x = "Date.fixed", y = "µM.CRE.mg.protein", fill = "Treatment.initial",
                                palette = c("#FC4E07", "#00AFBB"), add = "none", ylab =expression("µM Cu Equivalents"~(~mgProtein^{-1})),
                                title = "Total Antioxidant Capacity")
#plot_subsequent.TAC <- plot_subsequent.TAC %>% ggadd(c("jitter"), color ="Treatment.initial") # Add box plot
plot_subsequent.TAC <- plot_subsequent.TAC %>% ggadd(c("boxplot", "jitter"), fill ="Treatment.initial") # Add box plot
plot_subsequent.TAC <- ggpar(plot_subsequent.TAC, ylim = c(0,1200))
plot_subsequent.TAC

plot_subsequent.TAC.2 <- ggviolin(Data.D.15.18.21, x = "Treatment.initial", y = "µM.CRE.mg.protein", fill = "Treatment.initial",
                                palette = c("#FC4E07", "#00AFBB"), add = "none", ylab =expression("log(µM Cu Equivalents)"~(~mgProtein^{-1}*mgAFDW^{-1})),
                                title = "Total Antioxidant Capacity")
#plot_subsequent.TAC <- plot_subsequent.TAC %>% ggadd(c("jitter"), color ="Treatment.initial") # Add box plot
plot_subsequent.TAC.2 <- plot_subsequent.TAC.2 %>% ggadd(c("boxplot", "jitter"), shape ="Treatment", fill ="white") # Add box plot
plot_subsequent.TAC.2 <- ggpar(plot_subsequent.TAC.2, ylim = c(0,1200))
plot_subsequent.TAC.2

# arrange plots
ALL.second.subseq.PLOTS <- ggarrange(plot_second.subes.resp, plot_subsequent.AFDW, plot_subsequent.TP, plot_subsequent.TAC.2, ncol = 4, nrow = 1)
ALL.second.subseq.PLOTS # view plot

# both plots
FINAL.PLOTS <- ggarrange(ALL.initial.PLOTS, ALL.second.subseq.PLOTS, ncol = 1, nrow = 2)
FINAL.PLOTS # view plot

# Saving output plots
ggsave(file="Output/Phys.Data_2.pdf", FINAL.PLOTS, width = 14, height = 8, units = c("in"))













