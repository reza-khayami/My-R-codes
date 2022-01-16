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
###spearman
x <-  c(0, 4, 3, 0, 1)
y <- c(37, 39, 38, 36.8, 37.5)
cor(x, y, method = "spearman")
head(x)
pheatmap(x, border_color = NA, clustering_distance_cols = "correlation", clustering_distance_rows = "correlation" ) # NA $ 0 SD
x <- na.omit(x)
x <- x[-2: -3,]
head(x)
pheatmap(x)
#test khodam!
#pheatmap(x, col = pala)
#pala <- colorRampPalette(c("black", "red","green"))(100)
pheatmap(x, border_color = NA, clustering_distance_cols = "correlation", clustering_distance_rows = "correlation")                 
cor(x[,c("HHEX", "HLXb9")])                         
cor(x[,c("ISL1", "HLXb9")])                         
ggplot(x, aes(HHEX, HLXb9))+ geom_point() + geom_smooth(method = "lm")                    
y <- data.frame(Weight = c(45, 55, 65, 75, 85),
                Height = c(165, 170, 175, 180, 150))
cor(y)
#spearman has more resistance to chert data!
rowMeans(x) 
###apply
apply(x, 1, var)
apply(x, 2, var)
apply(x, 2, mean)
head(x)
apply(x[3,], 1, var) 
#or
var(as.numeric(x[3,]))

apply(x, 1, min)
apply(x, 2, min, na.rm = T)
args(apply)
lapply(1:3, function(x){x ^ 2})
sapply(1:3, function(x){x ^ 2})
t.test(x[4:10, 1], x[11:32, 1]) 
# what if we had 100 genes ?
Mytest <- function(d){ #d = gene mored nazar
  t.test(x[4:10, d], x[11:32, d]) $p.value
}
Mytest(2)
sapply(1:ncol(x), Mytest)
#OR 
Mytest2 <- function(y){ 
  t.test(y[4:10], y[11:32]) $p.value }
apply(x, 2, Mytest2)
###### expression
sd(x[4:6, 1])
se <- function(x){
  sd(x)/sqrt(length(x))
}
#Or
#se <- function(x){
#return(sd(x)/sqrt(length(x)))
#}
a <- se(x[4:6, 1])
length(x[4:6, 1])
options(digits = 4 )
#mean
x$Sample <- rownames(x)
x$Sample <- substr(x$Sample, 1, nchar(x$Sample)-2)
x.m <- aggregate(. ~ Sample, data = x, FUN = mean)
head(x.m)
ggplot(x.m, aes(Sample, HLXb9, fill = Sample))+ geom_bar(stat = "identity")
#### khodam error bar keshidam!
x.d <- aggregate(. ~ Sample, data = x, FUN = sd)
x.d1 <- x.d[,-1]
x.d2 <- x.d1/sqrt(length(x.d1))
ggplot(x.m, aes(Sample, HLXb9, fill = Sample))+ geom_bar(stat = "identity") + geom_errorbar(aes( ymin = x.m[,2]-1.96*x.d2[, 1],                                    ymax = x.m[,2]+1.96*x.d2[, 1]),
                                                                                            width=.2, position=position_dodge(.9))

