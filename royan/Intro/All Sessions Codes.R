library(pheatmap)
library(ggplot2)
library(reshape)
library(grid)

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

setwd("~/Desktop/R Class/")
x = read.delim("~/Desktop/R Class/Endoderm.txt")
rownames(x) = x[,1]
x = x[,-1]
x = log2(x+1)
x <- na.omit(x)

pdf("Desktop/R Class/Heatmap.pdf")
pheatmap(x,fontsize_row=5,border_color=NA)
dev.off()


x$Gene = rownames(x)
googooli=x
p = ggplot(googooli, aes(x=DE.1, y=DE.2, label=Gene)) + geom_point()
p = p + geom_text()
p


y = x[,c("Gene","DE.1")]
p1 = ggplot(y, aes(x=Gene,y=DE.1,fill=Gene))+geom_bar(stat="identity")
p2 = p1 + ylab("Expression of genes in Defenitive Endoderm (log2)")

pdf("~/Desktop/R Class/Barplot.pdf")
p1
p2
dev.off()

x.cor <- cor(x)

y <- t(x)
y[,2:3]=y[,2:3]+runif(18,min=-0.001,max=+0.001)
y.cor <- cor(y)
pheatmap(x.cor)
pheatmap(y.cor)

ggplot(x,aes(NEUROD1, NKX6.1))+geom_point()
ggplot(x,aes(HLXb9, HHEX))+geom_point()


MyTest <- function(i) {
    t.test(x[4:10,i], x[11:32,i])$p.value
}

MyTest2 <- function(y) {
    t.test(y[4:10], y[11:32])$p.value
}

se <- function(x) sd(x)/sqrt(length(x))



pc <- prcomp(x)

pcx <- data.frame(pc$x)

pcx$Sample <- rownames(pcx)
pcx$Sample <- substr(pcx$Sample, 1, nchar(pcx$Sample)-2)
ggplot(pcx, aes(PC1,PC2,color=Sample))+geom_point(size=3)+theme_complete_bw()

pcr <- data.frame(pc$rotation)
pcr$Gene <- rownames(pcr)
ggplot(pcr, aes(PC1,PC2,label=Gene))+geom_text()+theme_complete_bw()

x$Sample <- substr(rownames(x), 1, nchar(rownames(x))-2)

y <- x[,c("HLXb9","Sample")]


AOV <- function(Gene) {
    y <- x[,c(Gene,"Sample")]
    colnames(y)[1]="Gene"
    anova(aov(Gene ~ Sample, y))[1,5]
}

sapply(genes, AOV)


x.t <- t(x)
x.t <- x.t[,-2:-3]
xc <- cor(x.t, method="spearman")
pdf("Heatmap.pdf",width=14,height=14)
pheatmap(xc)
dev.off()

pc <- prcomp(x)
pcx <- data.frame(pc$x)
pcx$Sample <- rownames(pcx)
pcx$Sample <- substr(pcx$Sample, 1, nchar(pcx$Sample)-2)
pdf("PCA.pdf",width=10,height=10)
ggplot(pcx, aes(PC1, PC2, color=Sample))+geom_point(size=3)+theme_complete_bw()
dev.off()


x.m <- melt(as.matrix(x[3:10,]))
colnames(x.m) <- c("Sample","Gene","Exp")
anova(aov(Exp~Gene+Sample,x.m))
