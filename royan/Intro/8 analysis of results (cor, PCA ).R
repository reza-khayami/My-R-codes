# if we want to know that how repeatable our experiments were done : our repeats were similar
# 1.correlation: cor = calculates data on columns and pressure etc should be on rows
# if high correlation was seen between replicates ( cor between replicates were higher than cor between other conditions) we've done our job well.
# highrarcial clustering:  if the replicates were next to each other we've done Ok
# or heatmap of cor(x)
# 2. PCA: replicates should be near each other ##the best##
# 3. two way anova could be used: because we have two variables: genes and replicates

#analysis of data 
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