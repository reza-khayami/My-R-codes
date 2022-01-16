
# # Install and Load Packages and Data ---------------------------------------------


# If you haven't already: Install a few necessary packages
# These are some of the packages you need to install. If you have already installed them, no need to install again. Say Yes to questions about using personal libraries and about updates. Say No to installing from source.

# Install packages from CRAN:
install.packages("pheatmap")
install.packages("reshape2")
install.packages("gplots")
install.packages("RColorBrewer")
install.packages("ggplot2")
require(reshape2)
require(pheatmap)
require(gplots)
require(RColorBrewer)
require(ggplot2)

# Install bioconductor packages:
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2") 
# Remember to load packages using library(package_name)
# 
# Get and set the working directory
# For this course I recommend that you do not change the working directory and that you store any files for input in this directory.

getwd() # shows the directory where R is currently looking for files and saving files to
setwd("~/R.working/saved codes/RNAseq oslo university/data/") # You can change the working directory
# Reading in data
# Reading in the count data from the exercises. 
# Downloaded from here:http://folk.uio.no/jonbra/MBV-INF4410_2017/R/Mnemiopsis_count_data.txt

Mnemiopsis_count_data = read.table(file = "Mnemiopsis_count_data.txt", header = T, sep = "\t")

# And some metadata needed for later, downloaded from here: 
# http://folk.uio.no/jonbra/MBV-INF4410_2017/R/Mnemiopsis_col_data.txt
Mnemiopsis_col_data = read.table(file = "Mnemiopsis_col_data.txt", header = T, sep = "\t")

# Inspect the first 8 lines of the files
head(Mnemiopsis_count_data, 8)
head(Mnemiopsis_col_data, 8)


# # Quality Assessment ----------------------------------------------------


# Transformation and Visualization-------------------------------------
# Visualizing the between-sample distribution of counts
# If we plot the counts directly, the few extreme outliers (high counts) will dominate the plot completely.

boxplot(Mnemiopsis_count_data)
# And by plotting a histogram we see that the majority of genes have zero or very low counts

hist(Mnemiopsis_count_data[,1]) # Plotting only the first sample (column 1)
# It is therefore common to log2-transform counts before visualization

pseudoCount <-  log2(Mnemiopsis_count_data + 1) # log-transform to make numbers on scale (+1 to avoid zeroes)

boxplot(pseudoCount)
# And the histogram again

hist(pseudoCount[,1])
# We can use the ggplot2 package to make nicer plots :-) But it is also much more complicated and probably requires some googling around.


pseudoCount <-  as.data.frame(pseudoCount)
df  <-  melt(pseudoCount, variable.name = "Samples", value.name = "count") # reshape the matrix 

df <-  data.frame(df, Condition = substr(df$Samples, 1, 4))
ggplot(df, aes(x = Samples, y = count, fill = Condition)) +
  geom_boxplot() + 
  xlab("") +
  ylab(expression(log[2](count + 1)))

# We can also visualise the count distributions in a density plot.

ggplot(df, aes(x = count, colour = Samples, fill = Samples)) + ylim(c(0, 0.17)) +
  geom_density(alpha = 0.2, size = 1.25) + 
  facet_wrap(~ Condition) +
  theme(legend.position = "top") +
  xlab(expression(log[2](count + 1)))


# DESeq2 offers transformations for count data that stabilize the variance across the mean: the regularized logarithm (rlog) and the variance stabilizing transformation (VST). These have slightly different implementations, discussed a bit in the DESeq2 paper and in the very extensive web tutorial, but a similar goal of stablizing the variance across the range of values. Both produce log2-like values for high counts. Here we will use the regularized log transformation implemented with the rlog function.

# But first we need to create a DESeqDataSet which is the core object of DESeq2. There are several ways to create this depending on the type of input data.

library(DESeq2) # load the DESeq2 package
dds = DESeqDataSetFromMatrix(countData = Mnemiopsis_count_data,
                             colData = Mnemiopsis_col_data,
                             design = ~ condition) 
# we're testing for the different condidtions

dds

# Do the rlog transformation:
  rld <- rlogTransformation(dds)
  rld2 <- rlog(dds)
# See the effect of transformation
par( mfrow = c( 1, 2 ) )
plot(log2( 1 + counts(dds)[ , 1:2] ),
     pch=16, cex=0.3, main = "log2")
plot(assay(rld)[ , 1:2], # The assay function returns the count values of rld in a matrix
     pch=16, cex=0.3, main = "rlog")

# Clustering the sample-to-sample distances----------------------------
library("RColorBrewer") # Load a package giving more colors
library("pheatmap") # load a package for making heatmaps

distsRL <- dist(t(assay(rld))) # Calculate distances using transformed (and normalized) counts
?dist
mat <- as.matrix(distsRL) # convert to matrix
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition)) # set rownames in the matrix
colnames(mat) = NULL # remove column names

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) # set colors
pheatmap(mat,
         clustering_distance_rows=distsRL,
         clustering_distance_cols=distsRL,
         col=colors)
# PCA plot-------------------------------------------------------------
plotPCA(rld, intgroup=c("condition"))

# Systematic bias------------------------------------------------------
# We should also investigate any systematic bias in the sequencing data, such as whether one sample has been sequenced more deeply than others.

plot(pseudoCount[,4], pseudoCount[,7], pch = 20, xlab = "Aboral 4", ylab = "Oral 3") + # pch specifies the type of symbols. Filled dots in this case.
  abline(0,1, col = "red") # make a diagonal line

# We can also do this using ggplot2:
  
  x <-  pseudoCount[, 4] # extract fourth column
y <- pseudoCount[, 7] # extract seventh column
M <- x - y # M-values (differences between samples)
A <- (x + y)/2 # A-values (averages)
df <- data.frame(A, M)
ggplot(df, aes(x = A, y = M)) +
  geom_point(size = 1.5, alpha = 1/5) +
  geom_hline(aes(yintercept = 0), color = "blue3") +
  stat_smooth(se =                                                                 FALSE, method = "loess", color = "red3") + 
  xlab("Average gene count") +
  ylab("Count difference") +
  ggtitle("Aboral 4 vs. Oral 3")

# Normalization--------------------------------------------------------
# Normalization is done by DESeq2. It will estimate a size factor (scaling factor) which all the genes in a sample will be multiplied with. Ideally the size factor should be 1, which means that no normalization will take place.

dds <-  estimateSizeFactors(dds)
sizeFactors(dds)

# We can plot the normalized counts
# (compare with the previous box plot):
  
  norm_counts <-  counts(dds, normalized = TRUE) # Extract the normalized counts
pseudoCount <-  log2(norm_counts + 1) # convert to log-scale for visualization
df <- melt(pseudoCount) # transpose the matrix
df <- data.frame(df, Condition = substr(df$Var2, 1, 4))

ggplot(df, aes(x = Var2, y = value, fill = Condition)) + geom_boxplot() + xlab("") +
  ylab(expression(log[2](count + 1)))

# Differential expression analysis-------------------------------------

# As we have already specified an experimental design when we created the DESeqDataSet, we can run the differential expression pipeline on the raw counts with a single call to the function DESeq:
  dds<- DESeq(dds)
# Building the results table
res <- results(dds)
a <- mcols(res, use.names=TRUE)
summary(res)
# MA plot
plotMA(res, ylim=c(-7,7))
# Red points have adjusted p value < 0.1.

# We also see that there is a lot of noise associated with log2 fold changes from low count genes. DESeq2 can produce shrunken log2 fold changes to reduce the noice.

resShrink <- lfcShrink(dds, coef=2)
plotMA(resShrink, ylim=c(-5,5))
# Filtering the results
# Extract only the significant genes (padj < 0.1) from res

resSig <- subset(res, padj < 0.1)
resSig

# Show the 10 most strongest up-regulated in aboral based on fold change (i.e. most negative fold change because test vas oral vs. aboral)

head(resShrink[ order(resShrink$log2FoldChange), ], 10)

# Plot the top gene
plotCounts(dds, "ML327424a", "condition")

# And the 10 most down-regulated
head(resShrink[ order(resShrink$log2FoldChange, decreasing=TRUE), ], 10)

# And plot the top gene
plotCounts(dds, "ML34341a", "condition")

# Show the most highly significantly expressed genes (ordered by adjusted p-value)
head(res[order(res$padj),], 5) # order by padjusted and print the top 5

# And plot the top gene
plotCounts(dds, "ML087114a", "condition")

# There were a few thousand DE genes. This is perhaps too much? And from the MA-plot we see that when genes are highly expressed they are called DE even when the fold change is really small. It is possible to filter the results based on fold change. In this case, genes with a fold change of 2 (or -2) - i.e.doubling/halving.

resLFC1 <- results(dds, lfcThreshold = 1)
summary(resLFC1)
plotMA(resLFC1, ylim=c(-7,7)) +
  abline(h = 1, col = "blue") +
  abline(h = -1, col = "blue")

resLFC1Srhunk <- lfcShrink(dds, coef=2, res=resLFC1)
plotMA(resLFC1Srhunk, ylim=c(-5,5))+
  abline(h = 1, col = "blue") +
  abline(h = -1, col = "blue")

# We can also change the p.adjusted cutoff to more strict than the default (0.1) by changing the alpha (equivalent to FDR)

res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)
summary(res.05)
res.05Shrunk = lfcShrink(dds, coef = 2, res=res.05)
plotMA(res.05Shrunk, alpha = 0.05, ylim=c(-5,5)) +
  abline(h = 1, col = "blue") +
  abline(h = -1, col = "blue")

# Make a list of the 5 most significantly up-regulated genes in the aboral organ

resUpReg <- res[which(res$log2FoldChange < 0), ] # get the upregulated genes
head(resUpReg[order(resUpReg$padj),], 5) # order by padjusted and print the top 5

# Plot the top gene
plotCounts(dds, "ML327424a", "condition")

# Heatmap of the most significantly differentially expressed genes----
library("pheatmap")
?assay
mat <- assay(rld)[ head(order(res$padj),30), ] 
# select the top 30 genes with the lowest padj
mat <- mat - rowMeans(mat) # Subtract the row means from each value
# Optional, but to make the plot nicer:
df <- as.data.frame(colData(rld)[,c("condition")]) # Create a dataframe with a column of the conditions
colnames(df)  <-  "condition" # Rename the column header
rownames(df)  <-  colnames(mat) # add rownames
# and plot the actual heatmap
pheatmap(mat, annotation_col=df)
#Access the values of a specific gene----------------------------------
assay(rld)["ML327424a",]
