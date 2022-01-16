margin.table(HairEyeColor, 2)
eyes <- margin.table(HairEyeColor, 2)
round(prop.table(eyes), 2)
chi1 <- chisq.test(eyes)
chi1
browseURL("http://www.statisticbrain.com/eye-color-distribution-percentages/")
chi2 <- chisq.test(eyes, p=c(.41, .32, .15, .12))
chi2
##########################################################
are <- state.area
are
hist(are)
boxplot(are)
boxplot.stats(are)
summary(are)
mean(are)
mean(are, trim = 0.05)
###########################################################
A <- mtcars
mpg <- A$mpg
hp <- A$hp
qsec <- A$qsec
describe)mpg
describe(mpg)
curve(dnorm(x, mean = mean(mpg), sd = sd(mpg)), col = "Red")
lines(density(mpg))
      
B <- c(mpg,hp,qsec)
describe(B)
View(B)
G <- c(1,4,7)
?c
View(G)    
rm(list = ls())
A <- mtcars
B <- describe(A[, c(1,4,7)])
?describe
subset(mtcars, c("mpg", "hp", "qsec"))   
B[, c(3,4,11,12)]
G <- mtcars[, c(1,4,7)]
D <- describe(mtcars[, c(1,4,7)])
D[, c(3,4,11,12)]
