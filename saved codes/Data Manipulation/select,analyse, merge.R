####Selecting#####
mtcars <- mtcars
rm(list = ls())
mean(mtcars$qsec)
mean(mtcars$qsec[mtcars$cyl == 8])
median(mtcars$hp)
mean(mtcars$mpg[mtcars$hp > median(mtcars$hp)])
cyl.8a <- mtcars[mtcars$cyl == 8,]
mtcars[mtcars$cyl == 8 & mtcars$carb >= 4, ]
############################################################
####Analysing by subgroup###
iris <- iris
mean(iris$Petal.Width)
aggregate(iris$Petal.Width ~ iris$Species, FUN = mean)
?aggregate
aggregate(cbind(iris$Petal.Width, iris$Petal.Length) ~ iris$Species, FUN = mean)
?cbind
###########################################################
###Merging Files###
longley <- longley
data(longley)
a1 <- longley[1:14, 1:6]
a2 <- longley[1:14, 6:7]
b <- longley[15:16, ]
write.table(a1, "C:/Users/RK1994/Desktop/R/longley.a1.txt", sep = "\t")
write.table(a2, "C:/Users/RK1994/Desktop/R/longley.a2.txt", sep = "\t")
write.table(b, "C:/Users/RK1994/Desktop/R/longley.b.txt", sep = "\t")
a1t <- read.table("C:/Users/RK1994/Desktop/R/longley.a1.txt", sep = "\t")
a2t <- read.table("C:/Users/RK1994/Desktop/R/longley.a2.txt", sep = "\t")
a.1.2 <- merge(a1t,a2t, by = "Year") #for combining columns
b <- read.table("C:/Users/RK1994/Desktop/R/longley.b.txt", sep = "\t")
all.data <- rbind(a.1.2, b) #should have the same column names
row.names(all.data) <- NULL
