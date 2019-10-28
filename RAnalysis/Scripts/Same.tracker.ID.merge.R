#Title: Sample.tracker.ID.merge
#Author: Sam Gurr 
#Edited by: Sam Gurr
#Date Last Modified: 20191008
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 

#Read in required libraries
library(dplyr)          # Version 0.7.6, Packaged: 2018-06-27, Depends: R (>= 3.1.2)Imports: assertthat (>= 0.2.0), bindrcpp (>= 0.2.0.9000), glue (>=1.1.1), magrittr (>= 1.5), methods, pkgconfig (>= 2.0.1), R6(>= 2.2.2), Rcpp (>= 0.12.15), rlang (>= 0.2.0), tibble (>=1.3.1), tidyselect (>= 0.2.3), utils


#Required Data files

# Set Working Directory:
# setwd("~/MyProjects/Geoduck_Conditioning/RAnalysis/") #set working
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/")
#Load Sample Info
samples.TREAT.IDs <- read.csv(file="Sample.Tracker/Geoduck.sample.tracker.csv", header=T) #read sample.info data
samples.BOX.IDs <- read.csv(file="Sample.Tracker/Geoduck.Storage.Box.ID.csv", header=T) #read sample.info data
ID.reference.all <- read.csv(file="Data/Tank.ID.reference.subsequent.csv", header=T) #read sample.info data
ID.reference.all.condensed <- ID.reference.all %>% dplyr::select(Tank.ID, TREATMENT.ID.TOTAL)

# show error (duplicates) in BOX ID -- identified by photos of tubes
n_occur <- data.frame(table(samples.BOX.IDs$Tube.ID)) # list of IDs and number occured
n_occur[n_occur$Freq > 1,]

# NOTE: samples.TREAT.IDs is missing the total IDs for the last three days 
# can merge based on ID.reference.all
samples.TREAT.IDs.DAYS.all.other <-  samples.TREAT.IDs %>%  dplyr::filter(Date < 20190808) # all data BUT the target data for merge (rbind to this after)
unique(samples.TREAT.IDs.DAYS.all.other$Date) # check to see if it grabbed the correct data
samples.TREAT.IDs.DAYS15.21 <- samples.TREAT.IDs %>%  dplyr::filter(Date %in% c("20190808", "20190811", "20190814")) # select days without treatment ID
samples.TREAT.IDs.days15.21MERGE <- merge(samples.TREAT.IDs.DAYS15.21,
                                          ID.reference.all.condensed ,
                                          by = "Tank.ID") # merge with the condensed reference sheet by tank.ID to add treatments
samples.TREAT.IDs.days15.21MERGE$Treatment <- samples.TREAT.IDs.days15.21MERGE$TREATMENT.ID.TOTAL # duplicate to Treatment
samples.TREAT.IDs.days15.21MERGE <- samples.TREAT.IDs.days15.21MERGE %>% dplyr::select(-TREATMENT.ID.TOTAL) # remove to allow an rbind with the other data

samples.TREAT.IDs.FINAL <- rbind(samples.TREAT.IDs.DAYS.all.other,
                                 samples.TREAT.IDs.days15.21MERGE) # merge the two files

# now to bind witht the -80 box ID for analysis ID
samples.TREAT.IDs.FINAL # view the tube.ID column for the raw sample collection
samples.BOX.IDs # view the box reference for each tube

SAMPLE.TRACKER <- merge(samples.TREAT.IDs.FINAL,
                        samples.BOX.IDs,
                        by = "Tube.ID")

# write table output
write.table(SAMPLE.TRACKER,"C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/Sample.Tracker/Sample.Tracker.FINAL.csv",sep=",", row.names=FALSE)  # write table to output folder

