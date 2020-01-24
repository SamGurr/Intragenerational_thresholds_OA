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
library(car)
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
###################  IDENTIFY NEGATIVE TAOC #################### #
###################### AND OUTLIERS ############################ #
########## TO REPEAT WITH BIOLOGICAL REPLICATES ################ #
################ (~6 CUPS PER TREATMENT) ####################### #
################################################################ #

# negative values for all datasets
negative_TAOC_rows <- PHYS.DATA %>% dplyr::filter(µM.CRE.mg.protein < 0) # ID rows with TAOC < 0
negative_TAOC_rows # view rows with negative TAOC values - note these and assign more samples to homogenize
# NOTE: negative values for Day 21 and Day  data was repeated with biological replicates (animals fixed from same Tank.ID)
# ID 1483 (negative) replaced with 1636 -  day 21
# ID 1588 (negative) replaced with ID 1587 -  day 21
# ID 337 (negative) replaced with 338 -  day 7
# ID 450 (negative) replaced with 449 -  day 7


# find the outliers FOR TAOC on days 7 and 21
PHYS.OM <- dplyr::filter(PHYS.DATA, µM.CRE.mg.protein > 0) # ommit the negative values (above)
Day.7.PHYS <- PHYS.OM %>% dplyr::filter(Date.fixed  %in% 20190731) # day 7 data - negatives ommitted
boxplot(Day.7.PHYS$µM.CRE.mg.protein) # view boxplot
D.7.OUTLIERS <- boxplot(Day.7.PHYS$µM.CRE.mg.protein , plot=FALSE)$out # id the outlier(s)
D.7.OUTLIERS # view outlier(s) #  9.739665   0.263706 385.208675
min.values_TAOC_D.7 <- dplyr::filter(Day.7.PHYS, µM.CRE.mg.protein == (min(D.7.OUTLIERS)))
min.values_TAOC_D.7 # view row(s) of the outlier(s) # 349 replaced with 351
high.values_TAOC_D.7 <- dplyr::filter(Day.7.PHYS, µM.CRE.mg.protein > (max(D.7.OUTLIERS-1)))
high.values_TAOC_D.7 # view row(s) of the outlier(s) # 357 replaced with 358

Day.21.PHYS <- PHYS.OM %>% dplyr::filter(Date.fixed  %in% 20190814) # day 21 data - negatives ommitted
boxplot(Day.21.PHYS$µM.CRE.mg.protein) # view boxplot
D.21.OUTLIERS <- boxplot(Day.21.PHYS$µM.CRE.mg.protein , plot=FALSE)$out # id the outlier(s)
D.21.OUTLIERS # view outlier(s)
high.values_TAOC_D.21 <- dplyr::filter(Day.21.PHYS, µM.CRE.mg.protein > 170) # two other values that can be re done
high.values_TAOC_D.21 # view rows of outlier(s) - contains the two outliers in D.21.OUTLIERS and two other values ID as outliers in the three way anova
# ID 1438 replaced with 1437
# ID 1642 was not replaced
min.values_TAOC_D.21 <- dplyr::filter(Day.21.PHYS, µM.CRE.mg.protein == (min(D.21.OUTLIERS))) # two other values that can be re done
min.values_TAOC_D.21 # view rows of outlier(s) - contains the two outliers in D.21.OUTLIERS and two other values ID as outliers in the three way anova
# ID 1477 replaced with 1476

# CLEAN THE DATASETS OF THE OUTLIER VALUES (negatives already ommited)
Day.7.PHYS <- Day.7.PHYS %>%  dplyr::filter(!ID %in% c('349', '357')) # ommit IDs 349, 357 (337 and 450 are negatives already ommited)
Day.21.PHYS <- Day.21.PHYS %>%  dplyr::filter(!ID %in% c('1477', '1438', '1594')) # ommit IDs 1477, 1438 (1588 and 1483 are negatives already ommited)

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
Data.D.7.Resp <- resp.data.Days.1.7 %>%  dplyr::filter(Date %in% 20190731)
Data.D.7.Resp$Treatment.EXP_1 <- factor(Data.D.7.Resp$Treatment.EXP_1, levels = c("A", "M", "S"))
Data.D.7.Resp$Treat_hist <- substr(Data.D.7.Resp$Treatment.history, 1,1)
Data.D.7 <- Day.7.PHYS # assign the day 7 data to a new title
Data.D.7$Treat_history <- substr(Data.D.7$Treatment, 1,1) # Treat_history
Data.D.7$secondary <- substr(Data.D.7$Treatment, 3,3) # initial
Data.D.7$secondary <- factor(Data.D.7$secondary, levels = c("A", "M", "S")) # arrange levels allows ggpubr to plot correctly

# Resp  ##################################################################################################### #
# run the model on raw data
model_Day7_RESP<- aov(resp.COUNT.µg.L.hr.indiv  ~ Treat_hist * Treatment.EXP_1, data = Data.D.7.Resp)
summary(model_Day7_RESP) # no sig difference
shapiro.test(residuals(model_Day7_RESP)) # p-value = 0.4121; non normal via shapiro wilk test - likely due to single outlier
hist((residuals(model_Day7_RESP))) # histogram of residuals - left skew appears to have an outlier
qqnorm(residuals(model_Day7_RESP)) # qqplot
leveneTest(model_Day7_RESP) # p = 0.9121; homogenity of variance 
# plot the data
# D.7.PLOTboxRESP<- ggboxplot(Data.D.7.Resp, x = "Treatment.EXP_1", y = "resp.COUNT.µg.L.hr.indiv",  ylab = "µg.L.hr.indiv",  fill = "Treat_hist",
#                         palette = c( "#00AFBB", "#FC4E07"),add = "jitter", title ="Respiration rate Day 7", xlab = "Secondary pCO2 Treatment")
# D.7.PLOTboxRESP <- ggpar(D.7.PLOTboxRESP, ylim = c(0,12))
# D.7.PLOTboxRESP

D.7.PLOTboxRESP<-ggplot(Data.D.7.Resp, aes(x=Treatment.EXP_1, y=resp.COUNT.µg.L.hr.indiv, fill =Treatment.history, colour=Treatment.initial)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0, fill = "white") + 
  theme_classic() +
  labs(y=expression("Respiration rate"~(~µg~O[2]*hr^{-1}*individual^{-1})), x=expression("Secondary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treatment.history))), size = 2, position = position_jitterdodge(jitter.width = 0.5))+
  scale_colour_manual(values=c("skyblue1", "deepskyblue3","blue", "tomato1",  "red1", "firebrick4")) +
  scale_shape_manual(values=c(24, 21)) + 
  scale_fill_manual(values=c('white', 'white')) +
  ylim(0,14) + 
  scale_x_discrete(labels = c("Ambient","Moderate", "Severe")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
D.7.PLOTboxRESP


# TAOC ##################################################################################################### #
# data prep
# diagnostic
count.treat_1 <- Data.D.7 %>%  dplyr::group_by(Treatment) %>% # diagnostic of coutns for TAOC - still have more samples to run
  dplyr::summarise(count_replicated.by.day =n())
count.treat_1 # 6 for each 
count.treat_1.TANK <- Data.D.7 %>%  dplyr::group_by(Tank.ID) %>% # diagnostic of coutns for TAOC - still have more samples to run
  dplyr::summarise(count_replicated.by.day =n())
count.treat_1.TANK # 1 fir each tank replicate
# run the model on raw data
model_Day7_TAC <- aov(µM.CRE.mg.protein ~ Treat_history * secondary, data = Data.D.7)
summary(model_Day7_TAC) # no sig difference
shapiro.test(residuals(model_Day7_TAC)) # p-value =  0.5338; normal via shapiro wilk test
hist((residuals(model_Day7_TAC))) # histogram of residuals - left skew appears to have an outlier
qqnorm(residuals(model_Day7_TAC)) # qqplot
leveneTest(model_Day7_TAC) # p = 0.3634; homogenity of variance 
# plot the data
# D.7.PLOTbox<- ggboxplot(Data.D.7, x = "secondary", y = "µM.CRE.mg.protein",  ylab = "µM.CRE.mg.protein",  fill = "Treat_history",
#                              palette = c( "#00AFBB", "#FC4E07"),add = "jitter", title = "TAOC Day 7", xlab = "Secondary pCO2 Treatment")
# D.7.PLOTbox <- ggpar(D.7.PLOTbox, ylim = c(0,200))
# D.7.PLOTbox

D.7.PLOTbox <- ggplot(Data.D.7, aes(x=secondary, y=µM.CRE.mg.protein, fill = Treat_history, colour=Treatment)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0) + 
  theme_classic() +
  labs(y=expression("Total Antioxidant Capacity"~(~µM~CRE~Reducing~Equivalents^{-1}*mg~total~protein^{-1})), x=expression("Secondary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treat_history))), size = 2, position = position_jitterdodge(jitter.width = 1))+
  scale_colour_manual(values=c("skyblue1", "deepskyblue3","blue", "tomato1",  "red1", "firebrick4")) +
  scale_shape_manual(values=c(24, 21)) + 
  scale_fill_manual(values=c("white", "white")) +
  ylim(0,200) + 
  scale_x_discrete(labels = c("Ambient","Moderate", "Severe")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
D.7.PLOTbox

# Protein ##################################################################################################### #
# run the model on raw data
model_Day7_TotalProtein <- aov(mgProtein.mgAFDW  ~ Treat_history * secondary, data = Data.D.7)
summary(model_Day7_TotalProtein) # no significant difference
shapiro.test(residuals(model_Day7_TotalProtein)) # p-value = 0.7918; normal via shapiro wilk test
hist((residuals(model_Day7_TotalProtein))) # histogram of residuals
qqnorm(residuals(model_Day7_TotalProtein)) # qqplot
leveneTest(model_Day7_TotalProtein) # p = 0.1243; homogenity of variance 
# plot data
# D.7.PLOTbox.Protein<- ggboxplot(Data.D.7, x = "secondary", y = "mgProtein.mgAFDW",  ylab = "mgProtein.mgAFDW",  fill = "Treat_history",
#                         palette = c( "#00AFBB", "#FC4E07"),add = "jitter", title = "TP Day 7",xlab = "Secondary pCO2 Treatment")
# D.7.PLOTbox.Protein <- ggpar(D.7.PLOTbox.Protein, ylim = c(0,35))
# D.7.PLOTbox.Protein

D.7.PLOTbox.Protein<- ggplot(Data.D.7, aes(x=secondary, y=mgProtein.mgAFDW, fill = Treat_history, colour=Treatment)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0) + 
  theme_classic() +
  labs(y=expression("Total Protein"~(~µg~protein^{-1}*mg~AFDW^{-1})), x=expression("Secondary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treat_history))), size = 2, position = position_jitterdodge(jitter.width = 1))+
  scale_colour_manual(values=c("skyblue1", "deepskyblue3","blue", "tomato1",  "red1", "firebrick4")) +
  scale_shape_manual(values=c(24, 21)) + 
  scale_fill_manual(values=c("white", "white")) +
  ylim(0,35) + 
  scale_x_discrete(labels = c("Ambient","Moderate", "Severe")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
D.7.PLOTbox.Protein

# AFDW ##################################################################################################### #
# run the model on raw data
model_Day7_AFDW <- aov(mgTOTAL_AFDW  ~ Treat_history * secondary, data = Data.D.7)
summary(model_Day7_AFDW) # significant difference
shapiro.test(residuals(model_Day7_AFDW)) # p-value = 0.4076; normal via shapiro wilk test
hist((residuals(model_Day7_AFDW))) # histogram of residuals
qqnorm(residuals(model_Day7_AFDW)) # qqplot
leveneTest(model_Day7_AFDW) # p = 0.6027; homogenity of variance 
# POST HOC
TukeyHSD(model_Day7_AFDW, conf.level=0.95) # Treat_history E-A av diff = 2.608801 (E > A); p = 0.0047298 # S < A -2.349025 p = 0.0800328
# PLOTS
# D.7.PLOTbox.AFDW<- ggboxplot(Data.D.7, x = "secondary", y = "mgTOTAL_AFDW",  ylab = "mg AFDW",  fill = "Treat_history",
#                                 palette = c( "#00AFBB", "#FC4E07"),add = "jitter", title = "AFDW Day 7",xlab = "Secondary pCO2 Treatment")
# D.7.PLOTbox.AFDW <- ggpar(D.7.PLOTbox.AFDW, ylim = c(0,15))
# D.7.PLOTbox.AFDW

D.7.PLOTbox.AFDW<- ggplot(Data.D.7, aes(x=secondary, y=mgTOTAL_AFDW, fill = Treat_history, colour=Treatment)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0) + 
  theme_classic() +
  labs(y=expression("Ash Free Dry Weight"~(mg)), x=expression("Secondary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treat_history))), size = 2, position = position_jitterdodge(jitter.width = 1))+
  scale_colour_manual(values=c("skyblue1", "deepskyblue3","blue", "tomato1",  "red1", "firebrick4")) +
  scale_shape_manual(values=c(24, 21)) + 
  scale_fill_manual(values=c("white", "white")) +
  ylim(0,20) + 
  scale_x_discrete(labels = c("Ambient","Moderate", "Severe")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
D.7.PLOTbox.AFDW

##################################################################################################### #
#DAY 21 TEST ######################################################################################## #
##################################################################################################### #
##################################################################################################### #
# data prep
Data.D.21.Resp <- resp.data_Days.15.21 %>%  dplyr::filter(Date %in% 20190814)
Data.D.21.Resp$Treatment.EXP_1 <- factor(Data.D.21.Resp$Treatment.EXP_1, levels = c("A", "M", "S"))
Data.D.21.Resp$Treat_hist <- substr(Data.D.21.Resp$Treatment.history, 1,1)
Data.D.21.Resp$TOTAL <-paste(Data.D.21.Resp$Treatment.history, Data.D.21.Resp$Treatment.EXP_1, Data.D.21.Resp$Treatment.EXP_2 , sep="")
# PHYS data prep
# diagnostic
Data.D.21 <- Day.21.PHYS # assign the day 7 data to a new title
count.treat_2 <- Data.D.21 %>%  dplyr::group_by(Treatment) %>% # diagnostic of coutns for TAOC - still have more samples to run
  dplyr::summarise(count_replicated=n())
count.treat_2 # 6 for each 
count.treat_2.TANK <- Data.D.21 %>%  dplyr::group_by(Tank.ID) %>% # diagnostic of coutns for TAOC - still have more samples to run
  dplyr::summarise(count_replicated =n())
count.treat_2.TANK # 1  each tank replicate
count.treat_2.TANK %>% dplyr::filter(count_replicated > 1) #re there duplicates?
Day.21.PHYS %>%  dplyr::filter(Tank.ID %in% c('13','21', '27','55')) # look at the duplicated values
Day.21.PHYS.2 <- Day.21.PHYS %>%  dplyr::filter(!ID %in% c('1570', '1606', '1551', '1563')) # duplicates and are ommited based on duplication for later assay

Day.21.PHYS.2$Treat_history <- substr(Day.21.PHYS.2$Treatment, 1,1) # Treat_history
Day.21.PHYS.2$secondary <- substr(Day.21.PHYS.2$Treatment, 3,3) # secondary
Day.21.PHYS.2$tertiary <- substr(Day.21.PHYS.2$Treatment, 4,4) # tertiary
Day.21.PHYS.2$secondary <- factor(Day.21.PHYS.2$secondary , levels = c("A", "M", "S")) # arrange levels allows ggpubr to plot correctly

# Resp  ##################################################################################################### #
# run the model on raw data
model_Day21_RESP<- aov(resp.COUNT.µg.L.hr.indiv  ~ Treat_hist * Treatment.EXP_1 * Treatment.EXP_2, data = Data.D.21.Resp)
summary(model_Day21_RESP) # no sig difference
shapiro.test(residuals(model_Day21_RESP)) # p-value = 0.2533; non normal via shapiro wilk test - likely due to single outlier
hist((residuals(model_Day21_RESP))) # histogram of residuals - left skew appears to have an outlier
qqnorm(residuals(model_Day21_RESP)) # qqplot
leveneTest(model_Day21_RESP) # p = 0.5439; homogenity of variance 
# plot the data
# D.21.PLOTboxRESP<- ggboxplot(Data.D.21.Resp, x = "Treatment.EXP_1", y = "resp.COUNT.µg.L.hr.indiv",  ylab = "µg.L.hr.indiv",  fill = "Treat_hist",
#                             palette = c( "#00AFBB", "#FC4E07"),add = "jitter", shape = 'Treatment.EXP_2',  title ="Respiration rate Day 21", 
#                             xlab = "Secondary pCO2 Treatment")
# D.21.PLOTboxRESP <- ggpar(D.21.PLOTboxRESP, ylim = c(0,12))
# D.21.PLOTboxRESP

D.21.PLOTboxRESP<-ggplot(Data.D.21.Resp, aes(x=Treatment.EXP_2, y=resp.COUNT.µg.L.hr.indiv, fill = TOTAL, shape=TOTAL, colour=TOTAL)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0, fill = "white") + 
  theme_classic() +
  labs(y=expression("Respiration rate"~(~µg~O[2]*hr^{-1}*individual^{-1})), x=expression("Tertiary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treatment.history))), size = 2, position = position_jitterdodge(jitter.width = 0.1)) +
  scale_colour_manual(values=c("skyblue1", "skyblue1", "deepskyblue3", "deepskyblue3", 
                               "blue", "blue","tomato1","tomato1", 
                               "red1", "red1", "firebrick4","firebrick4")) +
  scale_shape_manual(values=c(24,1,2,21,
                              24,1,2,21,
                              24,1,2,21,24,2)) + 
  scale_fill_manual(values=c("white", "red1","white", "red1",
                             "white", "red1", "white", "red1", 
                             "white", "red1", "white", "red1","white", "red1")) +
  ylim(0,14) + 
  scale_x_discrete(labels = c("Ambient","Moderate")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
D.21.PLOTboxRESP

# TAOC ##################################################################################################### #
# model with raw data
model_Day21_TAC <- aov(µM.CRE.mg.protein ~ Treat_history * secondary * tertiary, data = Day.21.PHYS.2) # run model for TAOC raw data
summary(model_Day21_TAC) # marginal difference due to treatment history
shapiro.test(residuals(model_Day21_TAC)) # p-value = 0.03158;  normal via shapiro wilk test
hist((residuals(model_Day21_TAC))) # histogram of residuals - LEFT SKEW
qqnorm(residuals(model_Day21_TAC)) # qqplot
leveneTest(model_Day21_TAC) # p = 0.9197; homogenity of variance 
TukeyHSD(model_Day21_TAC) # E-A av. diff = -22.92171; p = 0.0062659
mean_TAC_D21.initial_effect <- Day.21.PHYS.2 %>% 
  dplyr::select(µM.CRE.mg.protein,Treat_history)  %>% 
  dplyr::group_by(Treat_history) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))
mean_TAC_D21.initial_effect # intital
# A = 106 ± 6.04
# E = 82.9 ± 5.19
((106-82.9)/106)*100 # percent difference 21.79245 %
# plot the data
# plot.DAY21.TAC <- ggboxplot(Day.21.PHYS.2, x = "Treat_history", y = "µM.CRE.mg.protein",  ylab = "µM.CRE.mg.protein", fill = "Treat_history",
#                                palette = c("#00AFBB", "#FC4E07"), add = "jitter",  title ="TAOC Day 21", xlab = "Initial pCO2 Treatment")
# plot.DAY21.TAC <- ggpar(plot.DAY21.TAC, ylim = c(0,200)) + geom_bracket(xmin = "A", xmax = "E", y.position = 220, label = "*", tip.length = 0.01)
# plot.DAY21.TAC # main effect plot
# D.21.PLOTbox<- ggboxplot(Day.21.PHYS.2, x = "secondary", y = "µM.CRE.mg.protein",  ylab = "µM.CRE.mg.protein",  fill = "Treat_history",
#                          palette = c( "#00AFBB", "#FC4E07"),add = "jitter", shape = "tertiary",
#                         xlab = "Secondary pCO2 Treatment", title = "TAOC Day 21")
# D.21.PLOTbox <- ggpar(D.21.PLOTbox, ylim = c(0,200)) 
# D.21.PLOTbox
# PLOTS.TAC.D21 <- ggarrange(plot.DAY21.TAC,D.21.PLOTbox, nrow = 1, widths = c(0.5, 1))
# PLOTS.TAC.D21 # view arranged plots

plot.DAY21.TAC <- ggplot(Day.21.PHYS.2, aes(x=tertiary, y=µM.CRE.mg.protein, fill = Treatment, shape=Treatment, colour=Treatment)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0, fill = "white") + 
  theme_classic() +
  labs(y=expression("Total Antioxidant Capacity"~(~µM~CRE~Reducing~Equivalents^{-1}*mg~total~protein^{-1})), x=expression("Tertiary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treat_history))), size = 2, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_colour_manual(values=c("skyblue1", "skyblue1", "deepskyblue3", "deepskyblue3", 
                               "blue", "blue","tomato1","tomato1", 
                               "red1", "red1", "firebrick4","firebrick4")) +
  scale_shape_manual(values=c(24,1,2,21,
                              24,1,2,21,
                              24,1,2,21,24,2)) + 
  scale_fill_manual(values=c("white", "red1","white", "red1",
                             "white", "red1", "white", "red1", 
                             "white", "red1", "white", "red1","white", "red1")) +
  ylim(0,200) + 
  scale_x_discrete(labels = c("Ambient","Moderate")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
plot.DAY21.TAC


# Protein ##################################################################################################### #
# run the model on raw data
model_Day21_Protein <- aov(mgProtein.mgAFDW  ~ Treat_history * secondary * tertiary, data = Day.21.PHYS.2)
summary(model_Day21_Protein)
shapiro.test(residuals(model_Day21_Protein)) # p-value = 3.147e-08; non normal via shapiro wilk test
hist((residuals(model_Day21_Protein))) # histogram of residuals left skew
qqnorm(residuals(model_Day21_Protein)) # qqplot - very bad outliers might be the reason for non normality
leveneTest(model_Day21_Protein) # p = 0.5965; homogenity of variance 
# look for outliers (seems like there are in qqnorm plot)
outliers.d21_protein <- boxplot(Data.D.21$mgProtein.mgAFDW , plot=FALSE)$out # id the outlier
print(outliers.d21_protein) # view the outlier(s)  26.62622 24.47385 50.82861 27.01706 62.79815
Data.D.21.OM <- Day.21.PHYS.2[-which(Day.21.PHYS.2$mgProtein.mgAFDW %in% outliers.d21_protein),] # remove outlier from new data frame
# run model with ommitted outliers
model_Day21_Protein.OM <- aov(mgProtein.mgAFDW  ~ Treat_history * secondary * tertiary, data = Data.D.21.OM)
summary(model_Day21_Protein.OM)
shapiro.test(residuals(model_Day21_Protein.OM)) # p-value = 0.7229; non normal via shapiro wilk test
hist((residuals(model_Day21_Protein.OM))) # histogram of residuals left skew
qqnorm(residuals(model_Day21_Protein.OM)) # qqplot
leveneTest(model_Day21_Protein.OM) # p = 0.4608; homogenity of variance 
# run model with LOG transformation
Day.21.PHYS.2$mgProtein.mgAFDW.LOG <- log(Day.21.PHYS.2$mgProtein.mgAFDW)
model_Day21_Protein.LOG <- aov(mgProtein.mgAFDW.LOG  ~ Treat_history * secondary * tertiary, data = Day.21.PHYS.2)
summary(model_Day21_Protein.LOG)
shapiro.test(residuals(model_Day21_Protein.LOG)) # p-value = 0.08939; non normal via shapiro wilk test
hist((residuals(model_Day21_Protein.LOG))) # histogram of residuals left skew
qqnorm(residuals(model_Day21_Protein.LOG)) # qqplot
leveneTest(model_Day21_Protein.LOG) # p = 0.09729; homogenity of variance 
# plot the data
# D.21.PLOTbox.Protein<- ggboxplot(Day.21.PHYS.2, x = "secondary", y = "mgProtein.mgAFDW",  ylab = "mgProtein.mgAFDW",  fill = "Treat_history",
#                                 palette = c( "#00AFBB", "#FC4E07"),add = "jitter",shape = "tertiary", title = "TP Day 21",
#                                xlab = "Secondary pCO2 Treatment")
# D.21.PLOTbox.Protein <- ggpar(D.21.PLOTbox.Protein, ylim = c(0,35))
# D.21.PLOTbox.Protein

D.21.PLOTbox.Protein <- ggplot(Day.21.PHYS.2, aes(x=tertiary, y=mgProtein.mgAFDW, fill = Treatment, shape=Treatment, colour=Treatment)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0, fill = "white") + 
  theme_classic() +
  labs(y=expression("Total Protein"~(~µg~protein^{-1}*mg~AFDW^{-1})), x=expression("Tertiary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treat_history))), size = 2, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_colour_manual(values=c("skyblue1", "skyblue1", "deepskyblue3", "deepskyblue3", 
                               "blue", "blue","tomato1","tomato1", 
                               "red1", "red1", "firebrick4","firebrick4")) +
  scale_shape_manual(values=c(24,1,2,21,
                              24,1,2,21,
                              24,1,2,21,24,2)) + 
  scale_fill_manual(values=c("white", "red1","white", "red1",
                             "white", "red1", "white", "red1", 
                             "white", "red1", "white", "red1","white", "red1")) +
  ylim(0,35) + 
  scale_x_discrete(labels = c("Ambient","Moderate", "Severe")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
D.21.PLOTbox.Protein

# AFDW ##################################################################################################### #
# run the model on raw data
model_Day21_AFDW <- aov(mgTOTAL_AFDW  ~ Treat_history * secondary * tertiary, data = Day.21.PHYS.2)
summary(model_Day21_AFDW) # significant difference
shapiro.test(residuals(model_Day21_AFDW)) # p-value = 9.291e-06; normal via shapiro wilk test
hist((residuals(model_Day21_AFDW))) # histogram of residuals - negative skewed
qqnorm(residuals(model_Day21_AFDW)) # qqplot
leveneTest(model_Day21_AFDW) # p = 0.4519; homogenity of variance 
TukeyHSD(model_Day21_AFDW, conf.level=0.95) # Treat_history E-A av diff = 3.346883 (E > A); p = 0.0004911

mean_AFDW_D21.tertiary_effect <- Day.21.PHYS.2 %>% 
  dplyr::select(mgTOTAL_AFDW,tertiary)  %>% 
  dplyr::group_by(tertiary) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))
mean_AFDW_D21.tertiary_effect # tertiary
# A = 5.19 ± 0.706
# M = 7.18 ± 1.04 
((7.18-5.19)/7.18)*100 # percent difference 27.71588 %

mean_AFDW_D21.initial_effect <- Day.21.PHYS.2 %>% 
  dplyr::select(mgTOTAL_AFDW,Treat_history)  %>% 
  dplyr::group_by(Treat_history) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,std.error))
mean_AFDW_D21.initial_effect # intital
# A = 4.03 ± 0.712
# E = 8.21 ± 0.920
((8.21-4.03)/8.21)*100 # percent difference 50.91352 %

# transform square root trans
Day.21.PHYS.2$mgTOTAL_AFDW.cbrt <- (Day.21.PHYS.2$mgTOTAL_AFDW)^(1/3)
model_Day21_AFDW.cbrt <- aov(mgTOTAL_AFDW.cbrt  ~ Treat_history * secondary * tertiary, data = Day.21.PHYS.2)
summary(model_Day21_AFDW.cbrt) # significant difference
shapiro.test(residuals(model_Day21_AFDW.cbrt)) # p-value = 0.1901; normal via shapiro wilk test
hist((residuals(model_Day21_AFDW.cbrt))) # histogram of residuals - negative skewed
qqnorm(residuals(model_Day21_AFDW.cbrt)) # qqplot
leveneTest(model_Day21_AFDW.cbrt) # p = 0.8736; homogenity of variance 
# POST HOC
TukeyHSD(model_Day21_AFDW.sqrt, conf.level=0.95) # Treat_history E-A av diff = 3.346883 (E > A); p = 0.0004911
# PLOTS
# D.21.PLOTbox.AFDW<- ggboxplot(Day.21.PHYS.2, x = "secondary", y = "mgTOTAL_AFDW",  ylab = "mgAFDW",  fill = "Treat_history",
#                              palette = c( "#00AFBB", "#FC4E07"),add = "jitter", shape = "tertiary", title = "AFDW Day 21",xlab = "Secondary pCO2 Treatment")
# D.21.PLOTbox.AFDW <- ggpar(D.21.PLOTbox.AFDW, ylim = c(0,25))
# D.21.PLOTbox.AFDW

D.21.PLOTbox.AFDW <- ggplot(Day.21.PHYS.2, aes(x=tertiary, y=mgTOTAL_AFDW, fill = Treatment, shape=Treatment, colour=Treatment)) +
  geom_boxplot(position=position_dodge(0.8), outlier.size = 0, fill = "white") + 
  theme_classic() +
  labs(y=expression("Ash Free Dry Weight"~(mg)), x=expression("Tertiary pCO"[2]~"Exposure")) +
  geom_point((aes(shape = factor(Treat_history))), size = 2, position = position_jitterdodge(jitter.width = 0.2)) +
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
  scale_x_discrete(labels = c("Ambient","Moderate", "Severe")) +
  # theme(legend.title = element_blank()) just ommit the legend title
  theme(legend.position = "none")
D.21.PLOTbox.AFDW






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
                    "white", "red1", "white", "red1", "white", "red1", 
                    "white", "red1", "white", "red1", "white", "red1")) +
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
TAC.TP.D7.PLOTS <- ggarrange(D.7.PLOTboxRESP, D.7.PLOTbox, D.7.PLOTbox.Protein, D.7.PLOTbox.AFDW, nrow = 1, ncol = 4, widths = c(1, 1, 1, 1))
TAC.TP.D7.PLOTS.2 <- ggarrange(FIG.schematic.data, TAC.TP.D7.PLOTS, nrow = 2)
TAC.TP.D7.PLOTS.2# view
TAC.TP.21.PLOTS <- ggarrange( D.21.PLOTboxRESP, plot.DAY21.TAC, D.21.PLOTbox.Protein, D.21.PLOTbox.AFDW, nrow = 1, ncol = 4, widths = c(1, 1, 1, 1))
TAC.TP.21.PLOTS.2 <- ggarrange(FIG.schematic.data, TAC.TP.21.PLOTS, nrow = 2)
TAC.TP.21.PLOTS.2# view
TAC.TP.D7.21_PLOTS <- ggarrange(D.7.PLOTboxRESP, D.7.PLOTbox, D.7.PLOTbox.Protein, D.7.PLOTbox.AFDW,
                                D.21.PLOTboxRESP, plot.DAY21.TAC, D.21.PLOTbox.Protein, D.21.PLOTbox.AFDW, nrow = 2, ncol = 4, widths = c(1, 1, 1, 1, 1, 1, 1, 1))
ggsave(file="Output/PHYSIOLOGY_plots_d7.and.d21.pdf", TAC.TP.D7.21_PLOTS, width = 12, height = 8, units = c("in")) # respiration rate plots
ggsave(file="Output/PHYSIOLOGY_plots_d7.pdf", TAC.TP.D7.PLOTS.2, width = 12, height = 8, units = c("in")) # respiration rate plots
ggsave(file="Output/PHYSIOLOGY_plots_d21.pdf", TAC.TP.21.PLOTS.2, width = 12, height = 8, units = c("in")) # respiration rate plots






























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













