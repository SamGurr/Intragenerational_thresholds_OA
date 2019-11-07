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
plot_initial.resp <- ggviolin(resp.data.Days.1.7.EFFECT,  ylim(1, 10) ,ylab =expression("Respiration rate"~(~µg~O[2]*hr^{-1}*indiv^{-1})), x = "Treatment.initial", y = "resp.MEAN.µg.L.hr.mm",  fill = "Treatment.initial",
                              palette = c("#FC4E07", "#00AFBB"), add = "none", title = "Respiration rate")
plot_initial.resp
?ggviolin
# stat.test.RESP<- compare_means(
# resp.MEAN.µg.L.hr.mm~Treatment.initial, data = resp.data.Days.1.7.EFFECT,
# method = "t.test")
# stat.test.RESP <- stat.test.RESP %>%# Add manually p-values from stat.test data
  mutate(y.position = c(9))

# Insert the results of the tukey test in the Resp.Analysis.R Script
# result from three way anova history × exposure = p = 0.0647; tukey test EH:M > AH:M p = 0.0408
threewayanova.resp.table <- data.frame(matrix(nrow = 1, ncol = 4))
colnames(threewayanova.resp.table)<-c('test', 'group1', 'group2', 'p-val_tukey') # column names, matches all data
threewayanova.resp.table[1,] <- c("initial", "EHM", "AHM", 0.0408)
threewayanova.resp.table <- threewayanova.resp.table %>%# Add manually p-values from stat.test data
  mutate(y.position = c(9))
plot_initial.resp <- plot_initial.resp %>% ggadd(c("boxplot", "jitter"), fill = "white")  + 
  stat_pvalue_manual(threewayanova.resp.table, label ="p-val_tukey") # Add box plot
plot_initial.resp <- ggpar(plot_initial.resp, ylim = c(0,9))

# INITIAL EFFECT - Marginal diff in metabolic response EHM > AHM during initial subseq exposure (Days 1 4 7)
# prep the data
Data.D.1.4.7 <- PHYS.DATA %>% dplyr::filter(Date.fixed < 20190801) # call data by DATE for the correct timeline
Data.D.1.4.7$Treatment <- factor(Data.D.1.4.7$Treatment, levels = c("EHM","AHM"))
count <- Data.D.1.4.7 %>%  dplyr::group_by(Treatment) %>% dplyr::summarise(count = n()) # 18 (AHM) and 17 (EHM) look at the count of replicates per day (removed NAs due to testing, under pipetting, and spec errors)
# ASH FREE DRY WEIGHT
Data.D.1.4.7$mgTOTAL_AFDW <- Data.D.1.4.7$TOTAL_AFDW*100
Model_1_AFDW <- aov(mgTOTAL_AFDW~Treatment*Exprmt.Day, data = Data.D.1.4.7) # run model 
summary(Model_1_AFDW)
shapiro.test(residuals(Model_1_AFDW)) # p-value = 0.0004303; non normal via shapiro wilk test
hist((residuals(Model_1_AFDW))) # histogram of residuals
qqnorm(residuals(Model_1_AFDW)) # qqplot
# transform to resolve normality assumptions
Data.D.1.4.7$mgASFW.LOG.TRANS <- log(Data.D.1.4.7$mgTOTAL_AFDW)
Model_1_AFDW.trans <- aov(mgASFW.LOG.TRANS~Treatment*Exprmt.Day, data = Data.D.1.4.7) # run model 
summary(Model_1_AFDW.trans)
shapiro.test(residuals(Model_1_AFDW.trans)) # p-value =  0.3181; non normal via shapiro wilk test
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
  mutate(y.position = c(2))
plot_intitial.AFDW <- plot_intitial.AFDW %>% ggadd(c("boxplot", "jitter"), fill = "white")  + stat_pvalue_manual(stat.test.AFDW, label = "p.adj") # Add box plot
plot_intitial.AFDW <- ggpar(plot_intitial.AFDW, ylim = c(0,3))
plot_intitial.AFDW

# TOTAL PROTEIN (corrected for AFDW)
Model_1_TP <- aov(mgProtein.mgAFDW~Treatment*Exprmt.Day, data = Data.D.1.4.7) # run model 
summary(Model_1_TP)
shapiro.test(residuals(Model_1_TP)) # p-value = 0.289; normal via shapiro wilk test
hist((residuals(Model_1_TP))) # histogram of residuals
qqnorm(residuals(Model_1_TP)) # qqplot
#plots
plot_intitial.TP <- ggviolin(Data.D.1.4.7, x = "Treatment", y = "mgProtein.mgAFDW", fill = "Treatment",
                             palette = c("#FC4E07", "#00AFBB"), ylab =expression("Total Protein"~(~mg_protein_mg_AFDW^{-1})),
                             add = "none", title ="mg.Total Protein")
plot_intitial.TP <- plot_intitial.TP %>% ggadd(c("boxplot", "jitter"), fill = "white")
plot_intitial.TP

# TOTAL ANTIOXIDANT CAPACIITY
Model_1_TAC <- aov(mMCuRedEqu.mgPro.mgAFDW ~Treatment*Exprmt.Day, data = Data.D.1.4.7) # run model 
summary(Model_1_TAC)
# test residuals for normality assumptions
shapiro.test(residuals(Model_1_TAC)) # p-value = 0.001308; non normal via shapiro wilk test
hist((residuals(Model_1_TAC))) # histogram of residuals
qqnorm(residuals(Model_1_TAC)) # qqplot
# call ouliers and remove
outliers <- boxplot(Data.D.1.4.7$mMCuRedEqu.mgPro.mgAFDW  , plot=FALSE)$out # call outliers
Data.D.1.4.7[which(Data.D.1.4.7$mMCuRedEqu.mgPro.mgAFDW  %in% outliers),] # call ouliers
Data.D.1.4.7.OM <- Data.D.1.4.7[-which(Data.D.1.4.7$mMCuRedEqu.mgPro.mgAFDW  %in% outliers),] # ommit ouliers new dataframe ouliers
# run model without ouliers
Model_1_TAC.OM <- aov(mMCuRedEqu.mgPro.mgAFDW ~Treatment*Exprmt.Day, data = Data.D.1.4.7.OM) # run model WITH NEGATIVE VALUES AS ZERO
summary(Model_1_TAC.OM)
# test residuals for normality assumptions
shapiro.test(residuals(Model_1_TAC.OM)) # p-value = 0.2552; normal via shapiro wilk test - removal of ouliers resolved
hist((residuals(Model_1_TAC.OM))) # histogram of residuals
qqnorm(residuals(Model_1_TAC.OM)) # qqplot
# post-hoc
Data.D.1.4.7.OM.MEANS <- lsmeans(Model_1_TAC.OM, pairwise ~  Treatment)# pariwise Tukey Post-hoc test between repeated treatments
# D.1.4.7.TAC # view post hoc summary
# Figure
plot_intitial.TAC <- ggviolin(Data.D.1.4.7.OM, x = "Treatment", y = "mMCuRedEqu.mgPro.mgAFDW ", fill = "Treatment",
                              palette = c("#FC4E07", "#00AFBB"), ylab =expression("µM Cu Equivalents"~(~mgProtein^{-1}*mgAFDW^{-1})), add = "none", title = "Total Antioxidant Capacity")
plot_intitial.TAC %>% ggadd(c("boxplot", "jitter"), color = "Treatment") 
stat.test.TAC <- compare_means(
  mMCuRedEqu.mgPro.mgAFDW ~Treatment, data = Data.D.1.4.7.OM,
  method = "t.test")
stat.test.TAC <- stat.test.TAC %>%# Add manually p-values from stat.test data
  mutate(y.position = c(3))
plot_intitial.TAC <- plot_intitial.TAC %>% ggadd(c("boxplot", "jitter"), fill = "white") + stat_pvalue_manual(stat.test.TAC, label = "p.adj") # Add box plot
plot_intitial.TAC <- ggpar(plot_intitial.TAC, ylim = c(-1,5))
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
Model_2_RESP.rate <- aov(resp.MEAN.µg.L.hr.mm~Treatment.EXP_2*Treatment.history*Treatment.EXP_1*Date, data = resp.data_Days.15.21) # run model 
summary(Model_2_RESP.rate)
shapiro.test(residuals(Model_2_RESP.rate)) # p-value = 0.4021; normal via shapiro wilk test
hist((residuals(Model_2_RESP.rate))) # histogram of residuals
qqnorm(residuals(Model_2_RESP.rate)) # qqplot
resp.data_Days.15.21$Treatment.initial <- paste(resp.data_Days.15.21$Treatment.history, resp.data_Days.15.21$Treatment.EXP_1, sep="")
resp.data_Days.15.21.EFFECT <- resp.data_Days.15.21 %>%  dplyr::filter(Treatment.initial %in% c("EHS", "EHA"))
plot_second.subes.resp <- ggviolin(resp.data_Days.15.21.EFFECT, x = "Treatment.initial", y = "resp.MEAN.µg.L.hr.mm", 
                                   ,ylab =expression("Respiration rate"~(~µg~O[2]*hr^{-1}*indiv^{-1})), fill = "Treatment.initial",
                                   palette = c("#FC4E07", "#00AFBB"), add = "none", title = "Respiration rate")
fourwayanova.resp.table <- data.frame(matrix(nrow = 1, ncol = 4))
colnames(fourwayanova.resp.table)<-c('test', 'group1', 'group2', 'p-val_tukey') # column names, matches all data
fourwayanova.resp.table[1,] <- c("initial", "EHS", "EHA", 0.0685)
fourwayanova.resp.table <- fourwayanova.resp.table %>%# Add manually p-values from stat.test data
  mutate(y.position = c(9))
plot_second.subes.resp <- plot_second.subes.resp %>% ggadd(c("boxplot", "jitter"), fill = "white")  +
  stat_pvalue_manual(fourwayanova.resp.table, label ="p-val_tukey") # Add box plot
plot_second.subes.resp 

# prep the PHYS ASSAY DATA
Data.D.15.18.21 <- PHYS.DATA %>% filter(Date.fixed > 20190807) 
count <- Data.D.15.18.21 %>%  dplyr::group_by(Treatment) %>% dplyr::summarise(count = n())
Data.D.15.18.21$Treatment.initial <- substr(Data.D.15.18.21$Treatment, 1,3)
Data.D.15.18.21$Treatment.initial <- factor(Data.D.15.18.21$Treatment.initial, levels = c("EHS", "EHA"))

# ASH FREE DRY WEIGHT
Data.D.15.18.21$mgTOTAL_AFDW <- Data.D.15.18.21$TOTAL_AFDW*100
Model_2_AFDW <- aov(mgTOTAL_AFDW~Treatment.initial*Treatment*Exprmt.Day, data = Data.D.15.18.21)
summary(Model_2_AFDW) # Treatment.initial:Exprmt.Day  2 0.0001755 8.775e-05   3.627 0.0377 *
# test residuals for normality assumptions
shapiro.test(residuals(Model_2_AFDW)) # p-value = 0.09616; NORMAL via shapiro wilk test
hist((residuals(Model_2_AFDW))) # histogram of residuals; left skew
qqnorm(residuals(Model_2_AFDW)) # qqplot
# post-hoc
TukeyHSD(Model_2_AFDW)
Model_2_AFDW.tr.POSTHOC <- lsmeans(Model_2_AFDW, pairwise ~  Treatment.initial:Exprmt.Day)# pariwise Tukey Post-hoc test between repeated treatments
Model_2_AFDW.tr.POSTHOC # view post hoc summary
# plot the data
plot_subsequent.AFDW <- ggviolin(Data.D.15.18.21, x = "Treatment.initial", y = "mgTOTAL_AFDW",  ylab = "mg AFDW",  fill = "Treatment.initial",
                                 palette = c("#FC4E07", "#00AFBB"),add = "none", title = "mgAFDW")
plot_subsequent.AFDW <- plot_subsequent.AFDW %>% ggadd(c("boxplot", "jitter"), fill = "white") # Add box plot
plot_subsequent.AFDW

# TOTAL PROTEIN (corrected for AFDW)
Model_2_TP <- aov(mgProtein.mgAFDW  ~Treatment.initial*Treatment*Exprmt.Day, data = Data.D.15.18.21)
summary(Model_2_TP) # Treatment.initial:Exprmt.Day  2 0.0001755 8.775e-05   3.627 0.0377 *
# test residuals for normality assumptions
shapiro.test(residuals(Model_2_TP)) # p-value = 0.1617; NORMAL via shapiro wilk test
hist((residuals(Model_2_TP))) # histogram of residuals; left skew
qqnorm(residuals(Model_2_TP)) # qqplot
# plot the data
plot_subsequent.TP <- ggviolin(Data.D.15.18.21, x = "Treatment.initial", y = "mgProtein.mgAFDW", fill = "Treatment.initial",
                               ylab =expression("Total Protein"~(~mg_protein_mg_AFDW^{-1})), palette = c("#FC4E07", "#00AFBB"),add = "none", title ="mgProtein.mgAFDW")
plot_subsequent.TP <- plot_subsequent.TP %>% ggadd(c("boxplot", "jitter"), fill = "white") # Add box plot
plot_subsequent.TP

# TOTAL ANTIOXIDANT CAPACIITY
# Run the model (3-way anova Treatment*Treatment*time)
Model_2_TAC <- aov((mMCuRedEqu.mgPro.mgAFDW ^(1/3))~Treatment.initial*Treatment*Exprmt.Day, data = Data.D.15.18.21)
summary(Model_2_TAC)
# test residuals for normality assumptions
shapiro.test(residuals(Model_2_TAC)) # p-value = 0.9556; NORMAL via shapiro wilk test
hist((residuals(Model_2_TAC))) # histogram of residuals; left skew
qqnorm(residuals(Model_2_TAC)) # qqplot
# post-hoc
TukeyHSD(Model_2_TAC)
Model_2_TAC.tr.POSTHOC <- lsmeans(Model_2_TAC, pairwise ~  Treatment.initial:Exprmt.Day)# pariwise Tukey Post-hoc test between repeated treatments
Model_2_TAC.tr.POSTHOC # view post hoc summary
# plot the data
plot_subsequent.TAC <- ggviolin(Data.D.15.18.21, x = "Treatment.initial", y = "mMCuRedEqu.mgPro.mgAFDW", fill = "Treatment.initial",
                                palette = c("#FC4E07", "#00AFBB"), add = "none", ylab =expression("µM Cu Equivalents"~(~mgProtein^{-1}*mgAFDW^{-1})),
                                title = "Total Antioxidant Capacity")
plot_subsequent.TAC <- plot_subsequent.TAC %>% ggadd(c("boxplot", "jitter"), fill = "white") # Add box plot
plot_subsequent.TAC


# arrange plots
ALL.second.subseq.PLOTS <- ggarrange(plot_second.subes.resp, plot_subsequent.AFDW, plot_subsequent.TP, plot_subsequent.TAC, ncol = 4, nrow = 1)
ALL.second.subseq.PLOTS # view plot

# both plots
FINAL.PLOTS <- ggarrange(ALL.initial.PLOTS, ALL.second.subseq.PLOTS, ncol = 1, nrow = 2)
FINAL.PLOTS # view plot

# Saving output plots
ggsave(file="Output/Phys.Data.pdf", FINAL.PLOTS, width = 14, height = 8, units = c("in"))













