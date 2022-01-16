require("datasets")
?lynx
data(lynx)
hist(lynx)
plot(lynx)
?hist
h <- hist(lynx,
          breaks = 11,
          freq = FALSE,
          col = "thistle1",
          main = "Histogram of Annual Canadian Lynx
          Trappings\n1821-1934",
          xlab = "Number of Lynx Trapped")
curve(dnorm(x, mean = mean(lynx), sd = sd(lynx)),
     col = "thistle4",
     lwd = 2,
     add = TRUE)