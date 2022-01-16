get(ws)
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
###PCA
pc <- prcomp(x)
x <- na.omit(x)
any(is.na(x))
summary(pc)
plot(pc)
pcx <- pc$x #matrix
pcx <- data.frame(pcx) #rahat tar
head(pcx)
ggplot(pcx, aes(PC1 , PC2)) + geom_point() 
pcx$Sample <- rownames(pcx)
head(pcx)
pcx$Sample <- substr(pcx$Sample, 1, nchar(pcx$Sample)-2)
ggplot(pcx, aes(PC1 , PC2, color = Sample)) + geom_point(size = 3)+ theme_complete_bw()
pcr <- pc$rotation
head(pcr)
barplot(x$PDX1)
pcr <- data.frame(pcr)
pcr$genes <- rownames(pcr)
ggplot(pcr, aes(PC1 , PC2, label = genes)) + geom_text()+
  theme_complete_bw()
install.packages("rgl")
require(rgl)
plot3d(pcx[, 1:3])