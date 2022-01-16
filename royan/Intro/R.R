# Load packages -----------------------------------------------------------
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
#Clean up
#cleanup <- theme(panel.grid.major = element_blank(),
#                panel.grid.minor = element_blank(),
#               panel.background = element_blank(),
#              axis.line = element_line(color = "black"))


# Load Data ---------------------------------------------------------------
x <- read.table("G:\\tutorial\\genetics and bioinformatics\\introduction to bioinformatics\\SharifiZarchi-bio\\Data\\Endoderm.txt", sep = "\t", header = TRUE) 
rownames(x) <- x[,"Type"]
head(x)
x <- x[,-1]
x <- log2(x+1)
x <- na.omit(x)
y <- t(x)
y <- data.frame(y)
# OR y <- t(x) %>% data.frame() 

####Intro to Graph
ggplot(y, aes(x = DE.1, y = DE.2)) + geom_point()
y$genes <- rownames(y)
head(y)
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


# Correlation -------------------------------------------------------------
rm(list= ls())
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

x <- read.table("G:\\tutorial\\genetics and bioinformatics\\introduction to bioinformatics\\SharifiZarchi-bio\\Data\\Endoderm.txt", sep = "\t", header = TRUE) 
rownames(x) <- x[,"Type"]

x <- x[,-1]
x <- log2(x+1)
x <- na.omit(x) # omits rows U should transpose with t() albeit not in this case
dim(x)
x.cor <- cor(x)
x.cor

y <- t(x)
y <- na.omit(y)
dim(y)
#OR y <- t(na,omit(t(y)))

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


# Gene expression Analysis ------------------------------------------------

##boxplot
require(reshape2)
x.m <- melt(x)
ggplot(x.m, aes(x = variable, y = value)) +
  geom_boxplot()
#gene expression distribiution in differet samples
x <- na.omit(x)
class(x) #data.frame
y <- t(x)
class(y) #matrix
head(y)
y.m <- melt(y) # = x1 <- as.matrix(x)
               #x2 <- melt(x1)
               #head(x2)
head <- melt(y)
head(y.m)
ggplot(y.m, aes(x = Var2, y = value)) +
  geom_boxplot()
#make first plot better
ggplot(x.m, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() + theme_classic()
#why y.m and x.m are different?
class(x)
class(y)
# omiting outliers
ggplot(x.m, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) 
#violin plot
ggplot(x.m, aes(x = variable, y = value, fill = variable)) +
  geom_violin() 
?Cairo

# Aplly functions ---------------------------------------------------------

###spearman cor
x <-  c(0, 4, 3, 0, 1)
y <- c(37, 39, 38, 36.8, 37.5)
cor(x, y, method = "spearman")
###
x <- read.table("G:\\tutorial\\genetics and bioinformatics\\introduction to bioinformatics\\SharifiZarchi-bio\\Data\\Endoderm.txt", sep = "\t", header = TRUE) 
rownames(x) <- x[,"Type"]
head(x)
x <- x[,-1]
x <- log2(x+1)
head(x)
pheatmap(x)
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


# Expression and Error bar ------------------------------------------------

#SD DE HLXb9
sd(x[4:6, 1])
se <- function(x){
  sd(x)/sqrt(length(x))
}
# because the function is just one line you can write it as follows:
# se <- function(x) sd(x)/sqrt(length(x))
#Another way:
#se <- function(x){
#e = (sd(x)/sqrt(length(x)))
#return(e)
#}
a <- se(x[4:6, 1])
length(x[4:6, 1])
options(digits = 4 )
#mean
head(x)
x$Sample <- rownames(x)
x$Sample <- substr(x$Sample, 1, nchar(x$Sample)-2)
#nchar returns lenghth of the charachter
#nchar("DE") = 2
#substr("Hello", 4, 5)
x.m <- aggregate(. ~ Sample, data = x, FUN = mean)
# . means everything
head(x.m)
ggplot(x.m, aes(Sample, HLXb9, fill = Sample))+ geom_bar(stat = "identity")

### khodam error bar keshidam!
x.d <- aggregate(. ~ Sample, data = x, FUN = sd)
head(x.d)
x.d1 <- x.d[,-1]
length(x.d1)
x.d2 <- x.d1/sqrt(length(x.d1)) # ehtemalan bayad makhraj 3 bashad jaye length(x.d1)
head(x.d2)
ggplot(x.m, aes(Sample, HLXb9, fill = Sample))+ geom_bar(stat = "identity") + geom_errorbar(aes( ymin = x.m[,2]-1.96*x.d2[, 1],                                    ymax = x.m[,2]+1.96*x.d2[, 1]),
                            width=.2, position=position_dodge(.9))



# PCA ---------------------------------------------------------------------


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
head(pcr)
ggplot(pcr, aes(PC1 , PC2, label = genes)) + geom_text()+ 
  theme_complete_bw()
install.packages("rgl")
require(rgl)
plot3d(pcx[, 1:3])

#### Anova
x$Sample <- substr(rownames(x), 1, nchar(rownames(x))-2)
y <- x[,c("HLXb9", "Sample")]
t.test(y[y$Sample == "SC",1], y[y$Sample == "DE",1])
#OR
#subset(y, Sample=="SC")
#subset(y, Sample=="SC")[,1]
#t.test(subset(y, Sample=="SC")[,1], subset(y, Sample=="DE")[,1])
#other examples
subset(x, HLXb9 > 2 * PDX1)
alaki <- data.frame(Gene = paste("Gene", 1:100), FC = runif(100,-10,10)) # paste data ra be ham michasbanad paste("hello" , "world") = "hello world"
subset(alaki, FC > 2)
subset(alaki, abs(FC) > 5)
###
var(subset(y, Sample=="SC")[,1])  #=[1] 0.1676977
var(subset(y, Sample=="DE")[,1]) # [1] 0.5520587
# variances are not equal and the greater variance is more than 2x of the lesser variance  # Brown-Forsythe & Welch F tests SHOULD BE CONSIDERED
# here anova is conducted just as example 
a <- aov(HLXb9 ~ Sample, data = y)
summary(a)
a
b <- anova(a)
b$`Pr(>F)`
b["Pr(>F)"]
#for every gene
x$Sample <- substr(rownames(x), 1, nchar(rownames(x))-2)
AOV <- function(Gene){
  y <- x[,c(Gene, "Sample")]
  colnames(y)[1] = "Gene"
  anova(aov(Gene ~ Sample, data = y))[1,5]} 
AOV("HLXb9")
genes <- colnames(x)
head(x)
genes <- genes[genes!="Sample"]
genes
#OR genes <- genes[-length(genes)]
sapply(genes, AOV)
ggplot(x, aes (Sample, HLXb9, fill = Sample))+ geom_boxplot() + +theme_complete_bw()
# sort graph
x$Sample <- factor(x$Sample, levels = unique(x$Sample))
x$Sample
ggplot(x, aes (Sample, HLXb9, fill = Sample))+ geom_boxplot()


# Analysis of Results -----------------------------------------------------



# if we want to know that how repeatable our experiments were done : our repeats were similar
# 1.correlation: cor = calculates data on columns and pressure etc should be on rows
# if high correlation was seen between replicates ( cor between replicates were higher than cor between other conditions) we've done our job well.
# highrarcial clustering:  if the replicates were next to each other we've done Ok
# or heatmap of cor(x)
# 2. PCA: replicates should be near each other ##the best##
# 3. two way anova could be used: because we have two variables: genes and replicates

#data analysis
#1. with correlation
x <- read.table("G:\\tutorial\\genetics and bioinformatics\\introduction to bioinformatics\\SharifiZarchi-bio\\Data\\Endoderm.txt", sep = "\t", header = TRUE) 
rownames(x) <- x[,"Type"]
head(x)
x <- x[,-1]
x <- log2(x+1)
x <- na.omit(x)
x.t <- t(x) # correlation is done on columns
cor(x.t) # two columns are 1
x.t <- x.t[,-2:-3]
xc <- cor(x.t) #Pearson
head(x.t)
pheatmap(xc) # DE3!
xc <- cor(x.t, method = "spearman")
pheatmap(xc) # SC1! chon varaiton kam , PP va PE !

# 2. with PCA
pc <- prcomp(x)
dim(pc$rotation) # genes
head(pc$x) # samples
pcx <- data.frame(pc$x)
pcx$samples <- rownames(pcx)
head(pcx)
pcx$sample <- substr(pcx$sample, 1, nchar(pcx$samples)-2)
ggplot(pcx, aes(PC1, PC2, color = sample)) + geom_point(size = 3)+
  theme_complete_bw()
#with ANOVA
head(x)
x.t <- t(x)
head(x.t)
#anova(aov(x~S1+S2+S3+..., x)) two way anova S = variable
#inja var1= genes, var2 = samples var3 = expression 3 columns
x.m <- melt(as.matrix(x)) # on data frame returns two columns on matrix reutrns 3 columns
head(x.m)
colnames(x.m) <- c("Sample", "Gene", "Exp") 
dim(x.m)
anova(aov(Exp ~ Gene + Sample, data = x.m))
x.m1 <- melt(as.matrix(x[4:6,])) # different repeats between one cell
colnames(x.m1) <- c("Sample", "Gene", "Exp") 
anova(aov(Exp ~ Gene + Sample, data = x.m1)) #  conclusion : genes are different while samples are not= OK!
x.m2 <- melt(as.matrix(x[3:10,])) # just one SC!
colnames(x.m2) <- c("Sample", "Gene", "Exp") 
anova(aov(Exp ~ Gene + Sample, data = x.m2)) # both pvalues are sig


# Non parametric ----------------------------------------------------------


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

