### Creating bar charts of group means ###
rm(list = ls())
spray <- InsectSprays
means <- aggregate(spray$count ~ spray$spray, FUN = mean)
plot(means)
mean.data <- t(means[-1])
colnames(mean.data) <- means[, 1]
barplot(mean.data)
barplot(means)
#################################################################
### Creating grouped box plots ###
require(MASS)
data(painters)
??MASS
boxplot(painters$Expression ~ painters$School, 
        col = brewer.pal(8, "Pastel2"),
        names = c("Renais.",
                  "Mannerist",
                  "Seicento",
                  "Venetian",
                  "Lombard",
                  "16th C.",
                  "17thC.",
                  "French"),
        boxwex = 0.5, 
        whisklty = 1,
        staplelty = 0,
       outpch = 16,
       outcol = brewer.pal(8, "Pastel2"),
       main = "REZA",
       xlab = "Painter's School",
       ylab = "Expression Ratings",
       border = "Red",
       col.lab = "Red",
       col.axis = "blue")
#clean up!
detach("package:MASS", unload = TRUE)
detach("package:RColorBrewer", unload = TRUE)
rm(list = ls())
###################################################################
###Scatter Plot###
data(cars)
a <- log2(cars+1)
plot(cars,
     pch = 16,
     col = "#A5A5A5",
     main = "REZA",
     xlab = "MPH",
     ylab = "Stopping Distance (feet)",
     col.lab = "Red")
#Linear regression line
abline(lm(cars$dist ~ cars$speed),
       col = "darkred",
       lwd = 2)
#locally weighted scatterplot smoothing
lines(lowess(cars$speed, cars$dist),
      col = "blue",
      lwd = 2)
require("car")
install.packages("car")   
scatterplot(cars$dist ~ cars$speed,
            pch = 16,
            col = "darkblue",
            main = "Reza",
            xlab = "speed",
            ylab =  "distance")

