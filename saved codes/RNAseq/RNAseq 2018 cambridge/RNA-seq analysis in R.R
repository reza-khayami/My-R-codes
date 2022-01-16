
# Pre-processsing RNA-seq data --------------------------------------------
setwd("C:/Users/RK1994/Documents/R.working/saved codes/RNAseq 2018 cambridge/")
# Data Import
library(DESeq2)
library(tidyverse)
require(dplyr)
require(data.table)
library(ggfortify)

# Read the sample information into a data frame
sampleinfo <- read.delim("data/SampleInfo.txt", stringsAsFactors=F)

# Read the data into R
seqdata <- read.delim("data/GSE60450_Lactation.featureCounts", 
                      comment = "#",
                      stringsAsFactors=F)
#Format Data

countdata <- seqdata %>%
  column_to_rownames("Geneid") %>% 
  # turn the geneid column into rownames
  rename_all(str_remove, ".bam") %>% 
  # remove the ".bam" from the column names
  select(sampleinfo$Sample) %>% 
  # keep sample columns using sampleinfo$Sample
  as.matrix()
head(countdata)

# countdata1 <- seqdata
# rownames(countdata1) <- countdata1$Geneid
# countdata1 <- countdata1[,-1]
# countdata1 <- countdata1[,6:17]
# colnames(countdata1) <- sampleinfo$Sample
# countdata1 <- as.matrix(countdata1)

#Filtering the genes

dim(countdata)
keep <- rowSums(countdata) > 5
countdata <- countdata[keep,]
dim(countdata)

#Quality assessment

##Library sizes bar plot

librarySizes <- colSums(countdata)
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, 
        main="Barplot of library sizes")
abline(h=20e6, lty=2)

#Count distribution boxplots
# Get log2 counts per million

logcounts <- log2(countdata + 1) 
# Count data is not normally distributed

# make a colour vector
statusCol <- as.numeric(factor(sampleinfo$Status)) +1
statusCol
# Check distributions of samples using boxplots
boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=statusCol)

# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(as.matrix(logcounts)), col="blue")

# Principle Component Analysis
rlogcounts <- rlog(countdata)
# run PCA
pcDat <- prcomp(t(rlogcounts))
# plot PCA
autoplot(pcDat)
# Lets add colour to look at the clustering for Status
autoplot(pcDat,
         data = sampleinfo, 
         colour="Status", 
         size=5)
# and now status
# Lets add colour to look at the clustering for Cell Type
autoplot(pcDat,
         data = sampleinfo, 
         colour="CellType", 
         size=5)
# We could use shape for one of the factors
autoplot(pcDat,
         data = sampleinfo, 
         colour="Status", 
         shape="CellType",
         size=5)
# Specify some clearer shapes to use that have a black outline and use fill
autoplot(pcDat,
         data = sampleinfo, 
         fill="Status", 
         shape="CellType",
         size=5) +
  scale_shape_manual(values=c(21, 24)) +
  guides(fill = guide_legend(override.aes=list(shape=22)))
# setting shape to FALSE causes the plot to default to using the labels
autoplot(pcDat,
         data = sampleinfo, 
         colour="CellType", 
         shape=FALSE,
         label.size=6)
?ifelse
sampleinfo <- sampleinfo %>% 
  mutate(CellType=ifelse(Sample=="MCL1.DG", "basal", CellType)) %>% 
  mutate(CellType=ifelse(Sample=="MCL1.LA", "luminal", CellType)) %>% 
  mutate(Group=str_c(CellType, ".", Group))

write_csv(sampleinfo, "results/SampleInfo_Corrected.txt")

autoplot(pcDat,
         data = sampleinfo, 
         fill="Status", 
         shape="CellType",
         size=5) +
  scale_shape_manual(values=c(21, 24)) +
  guides(fill = guide_legend(override.aes=list(shape=22)))
# PCA beyond the first two dimensions

autoplot(pcDat,
         data = sampleinfo, 
         fill = "Status", 
         shape = "CellType",
         size = 5,
         x = 2,
         y = 3) +
  scale_shape_manual(values=c(21, 24)) +
  guides(fill = guide_legend(override.aes=list(shape=22)))

# Another alternative:
# Interactive MDS Plot with Glimma
library(Glimma)
glMDSPlot(rlogcounts, 
          labels = sampleinfo$Sample, 
          groups = sampleinfo[,c("CellType", "Status")], 
          folder = "mds")

# Hierarchical clustering with heatmaps
#We don't want to plot a heatmap of all 22013 genes, so let's select data for the 500 most variable genes and plot the heatmap.

# We estimate the variance for each row in the logcounts matrix
countVar <- apply(rlogcounts, 1, var)
# Get the row numbers for the top 500 most variable genes
highVar <- order(countVar, decreasing=TRUE)[1:500]
# Subset logcounts matrix
hmDat <- rlogcounts[highVar,]

library(gplots)
library(RColorBrewer)
# Get some nicer colours
mypalette <- brewer.pal(11, "PiYG")
# http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$CellType]
# Plot the heatmap
heatmap.2(hmDat, 
          col=rev(morecols(50)),
          trace="column", 
          main="Top 500 most variable genes across samples",
          ColSideColors=col.cell,scale="row")
dev.off()

#Convert counts to DESeqDataSet object

# first lets check that our rows and columns match
all(sampleinfo$Sample == colnames(countdata))
# create the design formula
design <- as.formula(~ CellType)
# create the DESeqDataSet object
ddsObj <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = sampleinfo,
                                 design = design)
# Apply normalisation to DDS object
ddsObj <- estimateSizeFactors(ddsObj)
ddsObj@colData$sizeFactor
# The MCL1.LE and MCL1.LF have much smaller normalisation factors, and MCL1.LA and MCL1.LB have the largest

library(limma)
logcounts <- log2(countdata + 1)
par(mfrow=c(1,2))
plotMA(logcounts, array = 7)
abline(h=0,col="grey")
plotMA(logcounts, array = 11)
abline(h=0,col="grey")

normalizedCounts <- counts(ddsObj, normalized=TRUE) 
logNormalizedCounts <- log2(normalizedCounts + 1)
par(mfrow=c(1,2))
plotMA(logNormalizedCounts, array = 7)
abline(h=0,col="grey")
plotMA(logNormalizedCounts, array = 11)
abline(h=0,col="grey")

save(countdata, sampleinfo, file="results/preprocessing.RData")


# Differential Expression of RNA-seq data ---------------------------------



# load the RData object we created in the previous session
load("results/preprocessing.RData")
ls()

#Create a DESeqDataSet object with the raw data
#Creating the design model formula

# Use the standard R 'formula' syntax for an additive model
design <- as.formula(~ CellType + Status)

modelMatrix <- model.matrix(design, data = sampleinfo)
modelMatrix
#get virgin in intercept its nicer!
sampleinfo$Status <- factor(sampleinfo$Status, 
                            levels = c("virgin", "pregnant", "lactate"))
modelMatrix <- model.matrix(design, data = sampleinfo)
modelMatrix

# Build a DESeq2DataSet

# create the DESeqDataSet object
ddsObj <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = sampleinfo,
                                 design = design)

#Data exploration
# Let's plot a PCA from vst transformed data
vstcounts <- vst(ddsObj, blind=TRUE) 
plotPCA(vstcounts, intgroup=c("Status", "CellType"))

#Differential expression analysis with DESeq2
#The DESeq2 work flow

ddsObj <- estimateSizeFactors(ddsObj)
ddsObj <- estimateDispersions(ddsObj)
ddsObj <- nbinomWaldTest(ddsObj)

#The DESeq command
#In practice the 3 steps above can be performed in a single step using the DESeq wrapper function

# rebuild a clean DDS object
ddsObj <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = sampleinfo,
                                 design = design)
# Run DESeq
ddsObj <- DESeq(ddsObj)
res <- results(ddsObj, alpha=0.05)
head(res)
# By default, results has returned the contrast encoded by the final column in the model matrix. DESeq2 has the command resultsNames that allows us to view the contrasts that are available directly from the model matrix.

resultsNames(ddsObj)
#Let's just rename res so that we know which contrast results it contains.
resLvV <- res
rm(res)
resLvV
resPvV <- results(ddsObj, 
                  name="Status_pregnant_vs_virgin", 
                  alpha = 0.05)
resPvV
#Top 100
topGenesPvV <- as.data.frame(resPvV) %>%
  rownames_to_column("GeneID") %>% 
  arrange(padj) %>% 
  head(100)
topGenesPvV

resBvL <- results(ddsObj, 
                  name="CellType_luminal_vs_basal", 
                  alpha = 0.05)

topGenesBvL <- as.data.frame(resBvL) %>%
  rownames_to_column("GeneID") %>% 
  arrange(padj) %>% 
  head(200)

#Contrasts
#Suppose we want to find differentially expressed genes between pregnant and lactate. We don't have a parameter that explicitly will allow us to test that hypothesis. We need to provide a contrast.
resultsNames(ddsObj)

resPvL <- results(ddsObj,
                  contrast=c("Status", "pregnant", "lactate"), 
                  alpha = 0.05)
resPvL

# Comparing two design models
#Suppose we thought that maybe status were irrelevant and really the only differences might be between cell types. We could fit a simpler model, this would give us more degrees of freedom and therefore more power, but how would we know if it was a better model of not? We can compare the two models using the "log ratio test" (LRT).

designC <- as.formula(~ CellType )
# Compare the designs
ddsObjC <- DESeq(ddsObj, test="LRT", reduced=designC)
resCvCS <- results(ddsObjC)
resCvCS
#The null hypothesis is that there is no significant difference between the two models, i.e. the simpler model is sufficient to explain the variation in gene expression between the samples. If the the adjusted p-value for a gene passes a significance threshold (e.g. padj < 0.05) then we should consider using the more complex model for this gene. In practice we would usually choose one model or the other and apply it to all genes.

#Testing log2 fold change versus a threshold
#A common practice when considering the results of a differential expression analysis is to filter out genes that are statistically significant but have a low fold change. Perhaps you are only interested in very strong response to the experimental conditions and decide to eliminate from consideration all genes with an absolute fold change lower than 2x.

#Rather than do this by filtering after running the differential analysis, DESeq2 provides a means for testing the hypothesis that fold change is greater (or actually lesser than if you want) a given threshold.
resPvL2 <- results(ddsObj,
                   contrast=c("Status", "pregnant", "lactate"), 
                   alpha = 0.05,
                   lfcThreshold=0.5, 
                   altHypothesis="greaterAbs")
sum(resPvL2$padj<0.05, na.rm = T)
sum(resPvL$padj<0.05 & abs(resPvL$log2FoldChange)>=2^0.5, na.rm = T)
save(resLvV, ddsObj, sampleinfo, file="results/DE.RData")


# Clean Code --------------------------------------------------------------
# Read the sample information into a data frame
sampleinfo <- read.delim("data/SampleInfo.txt", stringsAsFactors=F)

# Read the data into R
seqdata <- read.delim("data/GSE60450_Lactation.featureCounts", 
                      comment = "#",
                      stringsAsFactors=F)
#Format Data

countdata <- seqdata %>%
  column_to_rownames("Geneid") %>% 
  # turn the geneid column into rownames
  rename_all(str_remove, ".bam") %>% 
  # remove the ".bam" from the column names
  select(sampleinfo$Sample) %>% 
  # keep sample columns using sampleinfo$Sample
  as.matrix()
head(countdata)
keep <- rowSums(countdata) > 5
countdata <- countdata[keep,]
# Principle Component Analysis
rlogcounts <- rlog(countdata)
# run PCA
pcDat <- prcomp(t(rlogcounts))

autoplot(pcDat,
         data = sampleinfo, 
         fill="Status", 
         shape="CellType",
         size=5) +
  scale_shape_manual(values=c(21, 24)) +
  guides(fill = guide_legend(override.aes=list(shape=22)))
# setting shape to FALSE causes the plot to default to using the labels
autoplot(pcDat,
         data = sampleinfo, 
         colour="CellType", 
         shape=FALSE,
         label.size=6)
design <- as.formula(~ CellType + Status)

modelMatrix <- model.matrix(design, data = sampleinfo)

ddsObj <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = sampleinfo,
                                 design = design)
# Run DESeq
ddsObj <- DESeq(ddsObj)
res <- results(ddsObj, alpha=0.05)
resLvV <- res
rm(res)
resLvV
resPvV <- results(ddsObj, 
                  name="Status_pregnant_vs_virgin", 
                  alpha = 0.05)
resPvV
#Top 100
topGenesPvV <- as.data.frame(resPvV) %>%
  rownames_to_column("GeneID") %>% 
  arrange(padj) %>% 
  head(100)
topGenesPvV

resBvL <- results(ddsObj, 
                  name="CellType_luminal_vs_basal", 
                  alpha = 0.05)

topGenesBvL <- as.data.frame(resBvL) %>%
  rownames_to_column("GeneID") %>% 
  arrange(padj) %>% 
  head(200)
save(resLvV, ddsObj, sampleinfo, file="results/DE.RData")
