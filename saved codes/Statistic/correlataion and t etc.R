rm(list = ls())
###calculating correlation
swiss <- swiss
cor(swiss)
round(cor(swiss), 2)
cor.test(swiss$Fertility, swiss$Education)
install.packages("Hmisc")
require("Hmisc")
?rcorr
rcorr(as.matrix(swiss))
as.data.frame(A)
detach("package:Hmisc", unload = TRUE)
?detach
###bivariate regression
trees <- trees
trees[1:5, ]
hist(trees$Height)
hist(trees$Girth)
plot(trees$Girth, trees$Height)
abline(lm(trees$Height ~ trees$Girth))
reg1 <- lm(Height ~ Girth, data = trees) #linear regression
summary(reg1)
confint(reg1)
?confint
predict(reg1)#predict height based on girth
predict(reg1, interval = "prediction") #CI for predicted height
lm.influence(reg1)
influence.measures(reg1)
### two sample t test- comparing means
sleep <- sleep
sd <- sleep[,1:2]
hist(sd$extra, col = "#E5E5E5")
boxplot(extra ~ group, data =sd)
boxplot(sd$extra ~ sd$group)
#independant two group t test
t.test(extra ~ group, data =sd)
#t test with options
t.test(extra ~ group,
       data =sd,
       alternative = "less", # one tailed
       conf.level = 0.80)
x <- rnorm(30, mean = 20, sd = 5)
y <- rnorm(30, mean = 22, sd = 5)
t.test(x, y)
### paired t test
t1 <- rnorm(50, mean = 52, sd = 6)
dif <- rnorm(50, mean = 6, sd = 12)
t2 <- t1 + dif
require("MASS")
pairs <- data.frame(t1, t2)
parcoord(pairs, var.label = TRUE)
t.test(t2, t1, paired = TRUE)
t.test(t2, t1,
       paired = TRUE, 
       mu = 6, # null value
       alternative = "greater", # one tailed
       conf.level = 0.99)
### one way ANOVA
x1 <- rnorm(30,40,8)
x2 <- rnorm(30,41,8)
x3 <- rnorm(30,45,8)
x4 <- rnorm(30,45,8)
boxplot(x1, x2, x3, x4)
xdf <- data.frame(cbind(x1, x2, x3, x4)) #xdf1 <- data.frame(x1, x2, x3, x4) is OK too
summary(xdf)
xs <- stack(xdf) #get one column with aoutcome and second column with groups
#conduct ANOVA one way
anova1 <- aov(values ~ ind, data = xs)
anova1
summary(anova1)
# post hoc comparisons
TukeyHSD(anova1)
?pairwise.t.test #other post hocs
?p.adjust #specific methods
### comparing proportions
n5 <- c(rep(100,5))
x5 <- c(65, 60, 60, 50, 45)
prop.test(x5, n5)
n2 <- c(40,40) #n trials
x2 <- c(30, 20) # n successes
prop.test(x2, n2, conf.level = 0.8)
###creating cross tabs for categorial variables
Titanic <- Titanic
ftable(Titanic)
Titanic
#concert table to data frame with one row per obs
?lapply
a <- as.data.frame.table(Titanic)
tdf <- as.data.frame(lapply(as.data.frame.table(Titanic), function(x)rep(x, as.data.frame.table(Titanic)$Freq)))[, -5]
#create contingency table
ttab <- table(tdf$Class, tdf$Survived)
ttab
#call also get cell, row, and column %
#With rounding to get just 2 decimal places
#Multiplied by 100 to make %
round(prop.table(ttab, 1), 2) * 100 #row 
round(prop.table(ttab, 2), 2) * 100 #column
round(prop.table(ttab), 2) * 100 #cell
#chi squared
tchi <- chisq.test(ttab)
tchi
tchi$observed
tchi$expected
tchi$residuals #pearson redidual
tchi$stdres #standardized residual
tchi$expected
### robust statistics for bivariate associations
help(package = "robust")
help(package = "robustbase")
help(package = "MASS") # robust rlm linear model
help(package = "quentreg") # quentile reg
install.packages = "quentreg"
require(quantreg)
?rq
data(engel)
attach(engel)
plot(income, # Create plot frame
     foodexp,
     xlab = "Household Income",
     ylab = "Food expenditure",
     type = "n",
     cex = .5)
points(income, #Points in plot
       foodexp, 
       pch = 16, 
       col = "lightgray",)
taus <- c(.05, .1, .25, .75, .9, .95) #Quantiles
xx <- seq(min(income), max(income), 100)
f <- coef(rq((foodexp) ~ (income), tau = taus)) #coefficients
yy <- cbind(1, xx)%*%f #Y values
for(i in 1:length(taus)){ #For each quantile value...
        lines(xx, yy[, i], col = "darkgray") #draw regression
}
abline(lm(foodexp ~ income), #standard linear reg
       col = "darkred",
       lwd = 2)
abline(rq(foodexp ~ income), #median reg
     col = "blue",
     lwd = 2)
legend(3000, 1000, #Plot legend
       c("mean fit", "medain fit"),
       col = c("darkred", "blue"),
       lty = 1,
       lwd = 2)
detach(engel)
detach("package:robust" , unload = TRUE)
detach("package:quantreg" , unload = TRUE)
detach("package:MASS" , unload = TRUE)
detach("package:rrcov" , unload = TRUE)
rm(list = ls())