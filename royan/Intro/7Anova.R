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
####
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