### Treating missing values with mice
rm(list = ls())
Prestige_miss <- read.csv("https://raw.githubusercontent.com/selva86/datasets/master/Prestige_miss.csv", sep = ,)
mydata <- (Prestige_miss)
head(mydata)
require(Hmisc)
mydata$education <- impute(mydata$education, mean)
mydata$type <- impute(mydata$type, mode)
install.packages("mice")
require(mice)
mydata <- Prestige_miss
micemod <- mice(mydata) #imputations
mydata2 <- complete(micemod, 1)
