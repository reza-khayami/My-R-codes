#zamani ke tozi normal nadarim az wilcox estefade mikonim ke non parametric ast
#vaghti nemidanim che tozi dare aval az qqplot estefade mikonim ta tozi ro befahmim
#example:
x <- c(5, 7, 1, 14, 25, 36, 12)
mean(x)
sd(x)
length(x)
y <- rnorm(7, mean = mean(x), sd= sd(x))
qqplot(x,y) # roye ye khat nist x normal distribution nadarad
#unpaired wilcox should be used
# anthoer example
x <- data.frame(person = 1:7, Pressure = runif(7,11,12), pressure2= runif(7, 10, 11))
wilcox.test(x$Pressure, x$pressure2, paired = T) #also known as manwitny
y <- data.frame(Mean = colMeans(x[,-1]), sd = apply(x[,-1], 2, sd))
y$Group <- rownames(y)
y
ggplot(y, aes(Group, Mean, ymin = Mean-sd, ymax = Mean + sd, fill = Group))+ geom_bar(stat = "identity") +
  geom_errorbar(width = 0.1)
#SE
se <- function(x) sqrt(sd(x))/(length(x)-1)
y <- data.frame(Mean = colMeans(x[,-1]), SE = apply(x[,-1], 2, se))
y
ggplot(y, aes(Group, Mean, ymin = Mean-SE, ymax = Mean + SE, fill = Group))+ geom_bar(stat = "identity") +
  geom_errorbar(width = 0.1)
