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
y <- t(x) %>% data.frame() 
x <- na.omit(x)
#OR y <- t(x)
#y <- data.frame(y)
ggplot(y, aes(x = DE.1, y = DE.2)) + geom_point()
y$genes <- rownames(y)
#or y <-cbind(x, genes = rownames(x))
p <- ggplot(y, aes(x = DE.1, y = DE.2, label = genes)) + geom_point() 
p
p <- p + geom_text()
#OR
ggplot(y, aes(x = DE.1, y = DE.2, label = genes)) + geom_point() + geom_text() 
#barplot
z <- y[, "DE.1", drop = F]
z <- y[, c("genes", "DE.1")]
ggplot(z, aes(x =genes , y = DE.1, fill = genes)) + geom_bar(stat = "identity") + ylab("Expresion of Genes \n in Defenitive Endoderm")

pdf("R.pdf")
ggplot(z, aes(x =genes , y = DE.1, fill = genes)) + geom_bar(stat = "identity") + ylab("Expresion of Genes \n in Defenitive Endoderm")
dev.off()
?pdf