
# -||| Microarray Analysis |||- -------------------------------------------

rm(list = ls())
getwd()
setwd("~/R.working/royan/Advanced1/Marray/")
res <- read.delim("results/results.txt", header = T)
aml.up <- subset(res, logFC >1 & adj.P.Val < 0.05 ) 
head(aml.up)
aml.up.gene <- unique(aml.up$Gene.symbol) # get the list of genes
length(aml.up.gene) # N of genes
dim(aml.up) # 528 probes



# -||| Load Packages |||- -------------------------------------------------

require(limma) # DE analysis
require(GEOquery) # dl date form geo
require(pheatmap)
require(ggplot2)
require(gplots)
require(dplyr)
require(reshape2)
require(Biobase)
cleanup <- theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(color = "black"))
#or 
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




# -||| Rscript in geo2r |||- ----------------------------------------------

series <- "GSE9476"
#if a series has multiple platforms you should note what platform you want to analysis

platform <- "GPL96"


# -||| load Series and Platform Data from GEO |||- ------------------------

gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/") # destdir = R saves the data in this address and wont download it every time

length(gset)
class(gset) #list
head(gset)

if (length(gset) > 1) idx <- grep(platform, attr(gset, "names"))   else idx <- 1
gset <- gset[[idx]] 
# or gset <- gset[[1]] 
#why [[]] and not []? [] returns list as a list [[]] returns as its true structure


# -||| Group Names for All Samples |||- -----------------------------------

gr <- c("CD34", rep("BM", 10), rep("CD34", 7), rep("AML", 26), rep("PB", 10), rep("CD34", 10))
length(gr)

ex <- exprs(gset)
# extracts expression matrix
dim(ex)

### log2scale if required

max(ex) 
# = 15.73261 no need for log2 here

min(ex) 
# 0.7730247

#  ex <- log2(ex+1)
#  exprs(gset) <- ex



# -||| Quality Control |||- -----------------------------------------------

pdf("results/boxplot.pdf", width = 64)
boxplot(ex)
dev.off()

#better with ggplot2 (selfmade!)
# pdf("results/boxplotggplot.pdf", width = 64)
# ex1 <- as.data.frame(t(ex))
# ex1 <- melt(ex)
# ex1$GR <- rep(gr, each = nrow(ex))
# ggplot(ex1, aes(Var2, value, fill = GR)) + geom_boxplot()+ cleanup
# dev.off()

 ### if Normalizing required  :   ex <- normalizeQuantiles(ex)
                                # exprs(gset) <- ex
                               #boxplot(ex)

#correlation heatmap  

pdf("Results/CorHeatmap.pdf", width= 15, height = 15)
pheatmap(cor(ex), labels_row = gr, labels_col = gr)
dev.off()
#pheatmap(cor(ex), labels_row = gr, labels_col = gr, color = redgreen(256), border_color = NA)

###PCA
#for genes

pc <- prcomp(ex)
pdf("Results/pc.pdf")
plot(pc)
plot(pc$x) # genes
dev.off()
#PCA on only expression didn't differentiate between genes with different expression very effective

#PCA based on difference from mean
# if a gene has high expression on both case and normal or has no epression on case and normal will be put aside

?scale
ex.scale <- t(scale(t(ex), scale = F))
#scale: if T x-mean/sd which we dont want if F x-mean  
# why transpose? beacuase scale can only work on columns and we want to scale genes not samples

pc2 <- prcomp(ex.scale)
plot(pc2)
plot(pc2$x[, 1:2])

#PCA for samples

pcr <- data.frame(pc2$rotation[, 1:3], Group = gr) 
# first three columns

pc2$rotation
head(pcr)
# if you want to add sth to matrix for example strings it will convert numerics to strings but data frame is not like that and each column can have different structure

pdf("results/PCA_samples.pdf") 
ggplot(pcr, aes(PC1, PC2, color = Group)) + geom_point()+theme_bw()
dev.off()

#PC1 has succesfully differntiated cd34 form AML and PB
#BM is near AML (AML should be hetrogenous cause its cancer!)
#data is good becuase different groups are differentiated well
#cd34 probably has two sub populations



# -||| DE Analysis |||- ---------------------------------------------------

# text vs factor? factor is string in appearance but actually is a number forexample men =1 women = 2 :

x <- c("MEN", "WOMEN")
y <- factor(x)
?factor
as.numeric(y)
levels(y)

gr <- factor(gr)
gset$description <- gr 
# gset description is the same as the group we defined (gr) but its text we need to make it as factor

design1 <- model.matrix(~ description + 0, gset) 
#creats a design matrix with entities in description and nothing else from gset
# AML+BM+CD34+PB 0= off 1= on rows are samples

colnames(design) <- levels(gr)
fit <- lmFit(gset, design) 
# fits linear model to data 
#R^2 = var(mean)- var(line)/ var(mean)

cont.matrix <- makeContrasts(AML-CD34, levels=design) 
# what are we comparing?

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01) 
# basyian model

tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
#fdr is Benjamini & Hochberg (False discovery rate)
#could be sorted by logFC or ...
#number = 250 first genes could make it Inf default is 250
#colnames came from annotation = AnnotGPL=TRUE in getGEO

head(tT$Gene.ID)
#gene ID is always constant but gene symbol could change 
#accession number could have multiple ids for one gene but gene ID is unique
#tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

head(tT)
tT <- subset(tT, select=c("Gene.symbol","Gene.ID","adj.P.Val","logFC"))
write.table(tT, "Results/AML_CD34.txt", row.names=F, sep="\t", quote = F)
head(tT)


# -||| Pathway Analysis |||- ----------------------------------------------

aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
dim(aml.up)


# in the final file there are genes with multiple names such as CALM3///CALM2///CALM1 which we don't want them to be like that so: 

# 1. we can delete those excessive names:
#aml.up.genes <- sub("///".*, "", aml.up.genes) .* means every letter  after that, double quote after that is the things we want to change the names to
# so in this example we change all /// and every letters after that to nothing

# 2. better way is to keep other names without /// :

aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))
write.table(aml.up.genes, file = "Results/AML_CD34_Up.txt", quote = F, row.names = F, col.names = F)


aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))
write.table(aml.down.genes, file = "Results/AML_CD34_Down.txt", quote = F, row.names = F, col.names = F)






