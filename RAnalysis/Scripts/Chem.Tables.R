#Title: Chem tables
#Project: FFAR
#Author: HM Putnam & Sam Gurr
#Edit by: Sam Gurr
#Date Last Modified: 20190909
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


#call cumulative spreadsheet of discrete seawater chemistry
chem<-read.csv("Output/Seawater_chemistry_table_Output_All.csv", header=T, sep=",", na.string="NA", as.is=T) 
chem # view table 

# seperate into tray and tank for pre experiment (~3 month conditioning in heath stack) and 21-day experiment (in tanks or cups)
chem.tray <- chem  %>% 
  dplyr::filter(Treatment %in% c("tray.Elevated", "tray.Ambient")) 
chem.tank <- chem %>% 
  dplyr::filter(Treatment %in% c("tank.Ambient", "tank.Severe", "tank.Moderate"))
#####################   make long table ################################################  #
chem.tray <- chem.tray %>% select(-("Tank")) # drop "Tank" to create long table via meltin next line
chem.tray.LONG <- melt(chem.tray, id.vars=c("Date", "Treatment")) # uses tidyr to make a long table from wide
chem.tank <- chem.tank %>% select(-("Tank")) # drop "Tank" to create long table via meltin next line
chem.tanks.LONG <- melt(chem.tank, id.vars=c("Date", "Treatment")) # uses tidyr to make a long table from wide
####################  divide data into the three experimental periods ################### #
# days 1 - 7
chem.tanks.LONG.1.to.7 <-  chem.tanks.LONG %>% 
  dplyr::filter(Date %in% (20190724:20190731))  
# days 8 - 14
chem.tanks.LONG.8.to.14 <-  chem.tanks.LONG %>% 
  dplyr::filter(Date %in% (20190801:20190807)) 
# days 15 - 21
chem.tanks.LONG.15.to.21 <-  chem.tanks.LONG %>% 
  dplyr::filter(Date %in% (20190808:20190816)) 
################### mean and SEM tables for the conditioning period  ################### # 
TABLE.pre <- ddply(chem.tray.LONG, c("Treatment","variable"), summarise, #Calculate descriptive stats by Treatment and exposure 
                     N = length(na.omit(value)), #count the sample size removing NA
                     mean = mean(value), #calculate average 
                     sem = sd(value)/sqrt(N)) # calcualte the standard error of the mean
TABLE.pre$mean <- signif(TABLE.pre$mean,digits=3) # reduce number of sig digits to three
TABLE.pre$sem <- signif(TABLE.pre$sem,digits=3) # reduce number of sig digits to three
TABLE.pre.WIDE <- reshape(TABLE.pre, idvar = c("Treatment", "N"), timevar = "variable", direction = "wide")
TABLE.pre.WIDE # view table
################### mean and SEM tables for each experimental period ################### # 
# days 1 - 7
TABLE.D.1.7 <- ddply(chem.tanks.LONG.1.to.7, c("Treatment","variable"), summarise, #Calculate descriptive stats by Treatment and exposure 
                         N = length(na.omit(value)), #count the sample size removing NA
                         mean = mean(value), #calculate average 
                         sem = sd(value)/sqrt(N)) # calcualte the standard error of the mean
TABLE.D.1.7$mean <- signif(TABLE.D.1.7$mean,digits=3) # reduce number of sig digits to three
TABLE.D.1.7$sem <- signif(TABLE.D.1.7$sem,digits=3) # reduce number of sig digits to three
TABLE.D.1.7.WIDE <- reshape(TABLE.D.1.7, idvar = c("Treatment", "N"), timevar = "variable", direction = "wide")
TABLE.D.1.7.WIDE # view table
# days 8 - 14
TABLE.D.8.14 <- ddply(chem.tanks.LONG.8.to.14, c("Treatment","variable"), summarise, #Calculate descriptive stats by Treatment and exposure 
                     N = length(na.omit(value)), #count the sample size removing NA
                     mean = mean(value), #calculate average 
                     sem = sd(value)/sqrt(N)) # calcualte the standard error of the mean
TABLE.D.8.14$mean <- signif(TABLE.D.8.14$mean,digits=3) # reduce number of sig digits to three
TABLE.D.8.14$sem <- signif(TABLE.D.8.14$sem,digits=3) # reduce number of sig digits to three
TABLE.D.8.14.WIDE <- reshape(TABLE.D.8.14, idvar = c("Treatment", "N"), timevar = "variable", direction = "wide")
TABLE.D.8.14.WIDE # view table
# days 15 - 21
TABLE.D.15.21 <- ddply(chem.tanks.LONG.15.to.21, c("Treatment","variable"), summarise, #Calculate descriptive stats by Treatment and exposure 
                     N = length(na.omit(value)), #count the sample size removing NA
                     mean = mean(value), #calculate average 
                     sem = sd(value)/sqrt(N)) # calcualte the standard error of the mean
TABLE.D.15.21$mean <- signif(TABLE.D.15.21$mean,digits=3) # reduce number of sig digits to three
TABLE.D.15.21$sem <- signif(TABLE.D.15.21$sem,digits=3) # reduce number of sig digits to three
TABLE.D.15.21.WIDE <- reshape(TABLE.D.15.21, idvar = c("Treatment", "N"), timevar = "variable", direction = "wide")
TABLE.D.15.21.WIDE # view table

############################ FINAL TABLES WITH ST.ERROR ###################################### #
# final table pre - 2 rows 12 columns
FINAL.TABLE.pre <- data.frame(matrix(nrow = 2, ncol = 1))
FINAL.TABLE.pre$N <- TABLE.pre.WIDE$N
FINAL.TABLE.pre$Treatment <- TABLE.pre.WIDE$Treatment
FINAL.TABLE.pre$Salinity <- paste(TABLE.pre.WIDE$mean.Salinity, TABLE.pre.WIDE$sem.Salinity, sep=" ± ")
FINAL.TABLE.pre$Temperature <- paste(TABLE.pre.WIDE$mean.Temperature, TABLE.pre.WIDE$sem.Temperature, sep=" ± ")
FINAL.TABLE.pre$pH <- paste(TABLE.pre.WIDE$mean.pH, TABLE.pre.WIDE$sem.pH, sep=" ± ")
FINAL.TABLE.pre$CO2 <- paste(TABLE.pre.WIDE$mean.CO2, TABLE.pre.WIDE$sem.CO2, sep=" ± ")
FINAL.TABLE.pre$pCO2 <- paste(TABLE.pre.WIDE$mean.pCO2, TABLE.pre.WIDE$sem.pCO2, sep=" ± ")
FINAL.TABLE.pre$HCO3 <- paste(TABLE.pre.WIDE$mean.HCO3, TABLE.pre.WIDE$sem.HCO3, sep=" ± ")
FINAL.TABLE.pre$CO3 <- paste(TABLE.pre.WIDE$mean.CO3, TABLE.pre.WIDE$sem.CO3, sep=" ± ")
FINAL.TABLE.pre$DIC <- paste(TABLE.pre.WIDE$mean.DIC, TABLE.pre.WIDE$sem.DIC, sep=" ± ")
FINAL.TABLE.pre$TA <- paste(TABLE.pre.WIDE$mean.TA, TABLE.pre.WIDE$sem.TA, sep=" ± ")
FINAL.TABLE.pre$Aragonite.Sat <- paste(TABLE.pre.WIDE$mean.Aragonite.Sat, TABLE.pre.WIDE$sem.Aragonite.Sat, sep=" ± ")
FINAL.TABLE.pre <- FINAL.TABLE.pre[,-1] # view table
# final table D.1.7 - 3 rows 12 columns
FINAL.TABLE.D.1.7 <- data.frame(matrix(nrow = 3, ncol = 1))
FINAL.TABLE.D.1.7$N <- TABLE.D.1.7.WIDE$N
FINAL.TABLE.D.1.7$Treatment <- TABLE.D.1.7.WIDE$Treatment
FINAL.TABLE.D.1.7$Salinity <- paste(TABLE.D.1.7.WIDE$mean.Salinity, TABLE.D.1.7.WIDE$sem.Salinity, sep=" ± ")
FINAL.TABLE.D.1.7$Temperature <- paste(TABLE.D.1.7.WIDE$mean.Temperature, TABLE.D.1.7.WIDE$sem.Temperature, sep=" ± ")
FINAL.TABLE.D.1.7$pH <- paste(TABLE.D.1.7.WIDE$mean.pH, TABLE.D.1.7.WIDE$sem.pH, sep=" ± ")
FINAL.TABLE.D.1.7$CO2 <- paste(TABLE.D.1.7.WIDE$mean.CO2, TABLE.D.1.7.WIDE$sem.CO2, sep=" ± ")
FINAL.TABLE.D.1.7$pCO2 <- paste(TABLE.D.1.7.WIDE$mean.pCO2, TABLE.D.1.7.WIDE$sem.pCO2, sep=" ± ")
FINAL.TABLE.D.1.7$HCO3 <- paste(TABLE.D.1.7.WIDE$mean.HCO3, TABLE.D.1.7.WIDE$sem.HCO3, sep=" ± ")
FINAL.TABLE.D.1.7$CO3 <- paste(TABLE.D.1.7.WIDE$mean.CO3, TABLE.D.1.7.WIDE$sem.CO3, sep=" ± ")
FINAL.TABLE.D.1.7$DIC <- paste(TABLE.D.1.7.WIDE$mean.DIC, TABLE.D.1.7.WIDE$sem.DIC, sep=" ± ")
FINAL.TABLE.D.1.7$TA <- paste(TABLE.D.1.7.WIDE$mean.TA, TABLE.D.1.7.WIDE$sem.TA, sep=" ± ")
FINAL.TABLE.D.1.7$Aragonite.Sat <- paste(TABLE.D.1.7.WIDE$mean.Aragonite.Sat, TABLE.D.1.7.WIDE$sem.Aragonite.Sat, sep=" ± ")
FINAL.TABLE.D.1.7 <- FINAL.TABLE.D.1.7[,-1] # view table
# final table D.8.14 - 1 row 12 columns
FINAL.TABLE.D.8.14 <- data.frame(matrix(nrow = 1, ncol = 1))
FINAL.TABLE.D.8.14$N <- TABLE.D.8.14.WIDE$N
FINAL.TABLE.D.8.14$Treatment <- TABLE.D.8.14.WIDE$Treatment
FINAL.TABLE.D.8.14$Salinity <- paste(TABLE.D.8.14.WIDE$mean.Salinity, TABLE.D.8.14.WIDE$sem.Salinity, sep=" ± ")
FINAL.TABLE.D.8.14$Temperature <- paste(TABLE.D.8.14.WIDE$mean.Temperature, TABLE.D.8.14.WIDE$sem.Temperature, sep=" ± ")
FINAL.TABLE.D.8.14$pH <- paste(TABLE.D.8.14.WIDE$mean.pH, TABLE.D.8.14.WIDE$sem.pH, sep=" ± ")
FINAL.TABLE.D.8.14$CO2 <- paste(TABLE.D.8.14.WIDE$mean.CO2, TABLE.D.8.14.WIDE$sem.CO2, sep=" ± ")
FINAL.TABLE.D.8.14$pCO2 <- paste(TABLE.D.8.14.WIDE$mean.pCO2, TABLE.D.8.14.WIDE$sem.pCO2, sep=" ± ")
FINAL.TABLE.D.8.14$HCO3 <- paste(TABLE.D.8.14.WIDE$mean.HCO3, TABLE.D.8.14.WIDE$sem.HCO3, sep=" ± ")
FINAL.TABLE.D.8.14$CO3 <- paste(TABLE.D.8.14.WIDE$mean.CO3, TABLE.D.8.14.WIDE$sem.CO3, sep=" ± ")
FINAL.TABLE.D.8.14$DIC <- paste(TABLE.D.8.14.WIDE$mean.DIC, TABLE.D.8.14.WIDE$sem.DIC, sep=" ± ")
FINAL.TABLE.D.8.14$TA <- paste(TABLE.D.8.14.WIDE$mean.TA, TABLE.D.8.14.WIDE$sem.TA, sep=" ± ")
FINAL.TABLE.D.8.14$Aragonite.Sat <- paste(TABLE.D.8.14.WIDE$mean.Aragonite.Sat, TABLE.D.8.14.WIDE$sem.Aragonite.Sat, sep=" ± ")
FINAL.TABLE.D.8.14 <- FINAL.TABLE.D.8.14[,-1] # view table
# final table D.15.21 - 2 rows 12 columns
FINAL.TABLE.D.15.21 <- data.frame(matrix(nrow = 2, ncol = 1))
FINAL.TABLE.D.15.21$N <- TABLE.D.15.21.WIDE$N
FINAL.TABLE.D.15.21$Treatment <- TABLE.D.15.21.WIDE$Treatment
FINAL.TABLE.D.15.21$Salinity <- paste(TABLE.D.15.21.WIDE$mean.Salinity, TABLE.D.15.21.WIDE$sem.Salinity, sep=" ± ")
FINAL.TABLE.D.15.21$Temperature <- paste(TABLE.D.15.21.WIDE$mean.Temperature, TABLE.D.15.21.WIDE$sem.Temperature, sep=" ± ")
FINAL.TABLE.D.15.21$pH <- paste(TABLE.D.15.21.WIDE$mean.pH, TABLE.D.15.21.WIDE$sem.pH, sep=" ± ")
FINAL.TABLE.D.15.21$CO2 <- paste(TABLE.D.15.21.WIDE$mean.CO2, TABLE.D.15.21.WIDE$sem.CO2, sep=" ± ")
FINAL.TABLE.D.15.21$pCO2 <- paste(TABLE.D.15.21.WIDE$mean.pCO2, TABLE.D.15.21.WIDE$sem.pCO2, sep=" ± ")
FINAL.TABLE.D.15.21$HCO3 <- paste(TABLE.D.15.21.WIDE$mean.HCO3, TABLE.D.15.21.WIDE$sem.HCO3, sep=" ± ")
FINAL.TABLE.D.15.21$CO3 <- paste(TABLE.D.15.21.WIDE$mean.CO3, TABLE.D.15.21.WIDE$sem.CO3, sep=" ± ")
FINAL.TABLE.D.15.21$DIC <- paste(TABLE.D.15.21.WIDE$mean.DIC, TABLE.D.15.21.WIDE$sem.DIC, sep=" ± ")
FINAL.TABLE.D.15.21$TA <- paste(TABLE.D.15.21.WIDE$mean.TA, TABLE.D.15.21.WIDE$sem.TA, sep=" ± ")
FINAL.TABLE.D.15.21$Aragonite.Sat <- paste(TABLE.D.15.21.WIDE$mean.Aragonite.Sat, TABLE.D.15.21.WIDE$sem.Aragonite.Sat, sep=" ± ")
FINAL.TABLE.D.15.21 <- FINAL.TABLE.D.15.21[,-1] # view table

# save output table
write.table(FINAL.TABLE.pre,"C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/Output/Chem.Table.1.pre.csv",sep=",", row.names=FALSE)  # write table to output folder
write.table(FINAL.TABLE.D.1.7,"C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/Output/Chem.Table.2.Days.1-7.csv",sep=",", row.names=FALSE)  # write table to output folder
write.table(FINAL.TABLE.D.8.14,"C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/Output/Chem.Table.3.Days.8-14.csv",sep=",", row.names=FALSE)  # write table to output folder
write.table(FINAL.TABLE.D.15.21,"C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/Output/Chem.Table.4.Days.15-21.csv",sep=",", row.names=FALSE)  # write table to output folder



