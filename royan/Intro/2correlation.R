rm(list = ls())
options(digits = 4 )
require(ggplot2)
require(pheatmap)
require(reshape2)
require(grid)
theme_complete_bw <- function(base_size = 24, base_family = "") 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      axis.line =         element_blank(),
      axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, colour = "black", vjust = 1),
      axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
      axis.ticks =        element_line(colour = "black"),
      axis.title.x =      element_text(size = base_size, vjust = 0.5),
      axis.title.y =      element_text(size = base_size, angle = 90, vjust = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      axis.ticks.margin = unit(0.1, "cm"),
      
      legend.background = element_rect(colour=NA), 
      legend.key =        element_rect(fill =NA, colour = "black", size = 0.25),
      legend.key.size =   unit(1.5, "lines"),
      legend.text =       element_text(size = base_size * 0.7),
      legend.title =      element_text(size = base_size * 0.8),
      legend.position =   "top",
      
      panel.background = element_rect(fill = "white", colour = NA), 
      panel.border =     element_rect(fill = NA, colour = "black", size=2), 
      panel.grid.major = element_line(colour = NA, size = 0.2), #"grey"
      panel.grid.minor = element_line(colour = NA, size = 0.5), #"grey"
      panel.margin =     unit(0.25, "lines"),
      
      strip.background = element_rect(fill = NA, colour = NA), 
      strip.text.x =     element_text(colour = "black", size = base_size * 0.8),
      strip.text.y =     element_text(colour = "black", size = base_size * 0.8, angle = +90),
      
      plot.background =  element_rect(colour = NA, fill = "white"),
      plot.title =       element_text(size = base_size*.8),
      plot.margin =      unit(c(1, 1, .5, .5), "lines"))
}

x <- read.table("I:\\anatomy\\R\\R introduction\\SharifiZarchi-bio\\Data\\Endoderm.txt", sep = "\t", header = TRUE) 
rownames(x) <- x[,"Type"]

x <- x[,-1]
x <- log2(x+1)
x <- data.frame(Pressure = c( 11.2, 11.5, 11, 11.1, 13, 14),
                Weight = c(70, 75, 67, 73, 100, 120))
cor(x$Pressure, x$Weight)
ggplot(x, aes(x = Pressure, y = Weight)) + geom_point()
x$Beat <-  c(70, 65, 80, 73, 58, 45)
cor(x$Pressure, x$Beat)
a <- ggplot(x, aes(x = Pressure, y = Beat)) + geom_point() +
  geom_smooth(method = "lm")

cor(x)
cor.test(x$Pressure, x$Weight)
cor.test(x$Pressure, x$Weight)
x$Alaki <- c(7, 3, 2, 5, 19, 1)
cor(x)
cor.test(x$Pressure, x$Alaki)
x$Alaki <- runif(6)
x$Alaki
cor.test(x$Pressure, x$Alaki)
###
x <- read.table("I:\\anatomy\\R\\R introduction\\SharifiZarchi-bio\\Data\\Endoderm.txt", sep = "\t", header = TRUE) 
rownames(x) <- x[,"Type"]

x <- x[,-1]
x <- log2(x+1)
x <- na.omit(x) # omits rows U should transpose with t() albeit not in this case
dim(x)
x.cor <- cor(x)
x.cor

####
y <- t(x)
y <- na.omit(y)
dim(y)
#OR y <- t(na,omit(t(y)))
####
head(x.cor)
require(pheatmap)
pheatmap(x.cor)
ggplot(x, aes(NEUROD1, NKX6.1)) + geom_point()
ggplot(x, aes(HLXb9, HHEX)) + geom_point() +
  geom_smooth(method = "lm")
y <- t(x)
y.cor <- cor(y)
y.cor <- na.omit(y.cor) #worst and simplest thing to do !
dim(y.cor)
dim(y[, 2:3])
y[, 2:3] <- y[, 2:3] + runif(18, min = -0.001, max = 0.001)
head(y)
y.cor <- cor(y)
pheatmap(y.cor)

summary(x[, 1])
summary(x)