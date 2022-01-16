rm(list = ls())
### creating clustered bar charts for means
data("warpbreaks")
data <- tapply(warpbreaks$breaks,
               list(warpbreaks$wool,
                    warpbreaks$tension), mean)
barplot(data, 
        beside = TRUE,
        col = c("steelblue3", "thistle3"),
        bor = NA,
        main = "Mean Number of Warp Breaks\nby Tension and Woll",
        xlab = "Tension",
        ylab ="Mean Number of Breaks")
legend(locator(1),
       rownames(data),
       fill = c("steelblue3", "thistle3"))
####scatter plot
data(iris)
install.packages("car")
require(car)
#scatterplot = sp
?sp
sp(Sepal.Width ~ Sepal.Length | Species,
   data = iris,
   xlab = "Sepal Width",
   ylab = "Sepal Length",
   main = "Iris Data",
  )
detach("package:car", unload = TRUE)

### scatterplot matrices
data(iris)
#basic sp matrix
pairs(iris[1:4])
require("RColorBrewer")
#put histograms on the diagonl (from "pairs" help)
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
h <- hist(x, plot = FALSE)
breaks <- h$breaks; nB <- length(breaks)
y <- h$counts; y <- y/max(y)
rect(breaks[-nB], 0, breaks[-1], y, ...)
# Removed "col" = "cyan" frome code block; original below
#rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
pairs(iris[1:4],
      panel = panel.smooth,
      main = "blah blah",
      diag.panel = panel.hist,
      pch = 16,
      col = brewer.pal(3, "Pastel1")[unclass(iris$Species)])
#Similar with "car" package
#Gices Kernel density and rugplot for each variable
require(car)
scatterplotMatrix(~Petal.Length + Petal.Width + Sepal.Length + Sepal.Width | Species,
                  data = iris,
                  col = brewer.pal(3, "Dark2"),
                  main = "blah blah")
#Clean up
palette("default")

###3d sp
data(iris)
