###Detecting and treating outliers
rm(list = ls())
data(Prestige, package = "car")
myvar <- Prestige$income
# 3 times IQR
iqr <- IQR(myvar)
myvar[myvar > (3*iqr)]
#Outside 1.5 times the 75th %tile
third.quantile <- quantile(myvar, 0.75)
myvar[myvar > (1.5*third.quantile)]
# asl outlier!
a <- IQR(myvar)
myvar[myvar > (1.5*a + quantile(myvar, 0.75))]

#capping
myvar <- Prestige$income
myvar[myvar > quantile(myvar, 0.95)] <- quantile(myvar, 0.95)
myvar[myvar > quantile(myvar, 0.05)] <- quantile(myvar, 0.05)
