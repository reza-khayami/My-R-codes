swiss <- (swiss)
data(swiss)
fertility <- swiss$Fertility
h <- hist(fertility,
          prob = TRUE,
          ylim = c(0, 0.04),
          xlim = c(30, 100),
          breaks = 11,
          col = "#E5E5E5",
          border = 0,
          main = "title")
curve(dnorm(x, mean = mean(fertility), sd = sd(fertility)),
      col = "Red",
      lwd = 3,
      add = TRUE)
lines(density(fertility), col = "Blue")
lines(density(fertility, adjust = 3), col = "darkgreen")
rug(fertility, col = "Red")