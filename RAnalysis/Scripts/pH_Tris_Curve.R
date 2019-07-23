# modified for personal Lenovo PC at PT Whitney summer 2019
# last modification on 20190625 by SJG
setwd("C:/Users/samjg/Documents/My_Projects/Inragenerational_thresholds_OA/RAnalysis/Data/pH_Calibration_Files/") #set working directory

Calib.Data <-read.table("20190717.csv", header=TRUE, sep=",", na.string="NA", as.is=TRUE) #reads in the data files
model <-lm(mVTris ~ TTris, data=Calib.Data) #runs a linear regression of mV as a function of temperature
coe <- coef(model) #extracts the coeffecients
R2<-summary(model)$r.squared
plot(mVTris ~ TTris, data=Calib.Data)
abline(lm(mVTris ~ TTris, data=Calib.Data))
legend('topleft', legend = bquote(R^2 == .(format(R2, digits = 3))), bty='n')
