library(car)
data(Prestige)
myVar <- Prestige$income

#different methods to deal with outliers

#3 times IQR
iqr <- IQR(myVar)
myVar[myVar > (iqr*3)]

# outside 1.5 times the .75th quartile
third_quartile <- quantile(myVar, 0.75)
myVar[myVar > third_quartile * 1.5]

#caping the outliers = replace the values above 95th precentile or 99th percentile with respective percentile values itself. Same applies for values below the 5th or 1st percentiles.

myVar[myVar > quantile(myVar,0.95)] <- quantile(myVar, 0.95)
myVar[myVar < quantile(myVar,0.05)] <- quantile(myVar, 0.05)

