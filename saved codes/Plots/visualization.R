data("hflights" , package = "hflights")
head(hflights)
hist(hflights$Month, )
a <- table(hflights$Month)
  barplot(a)
plot(hflights$Month)
plot
OD <- table(hflights$Origin , hflights$Month)
OD <- as.data.frame.matrix(OD)
plot(as.numeric(OD[1,]), type = "b", ylim = c(0, 20000), col = "#5566EE", ylab = "N of Dep", xlab = "Month")
text(x= 1:12, y = as.numeric(OD[1, ]+1000), labels = as.numeric(OD[1,]), cex = 0.5)
lines(as.numeric(OD[2, ]), type = "b", ylim = c(0, 20000))
as <- hflights[hflights$DepDelay > 0, ]
d <- boxplot (as$DepDelay ~ as$UniqueCarrier, ylim = c(0,100))
d$out
d$n
aw <- as.numeric(table(hflights$UniqueCarrier))
head(aw)
