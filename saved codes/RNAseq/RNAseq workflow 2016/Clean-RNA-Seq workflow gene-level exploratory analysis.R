#doi: 10.12688/f1000research.7035.1     read if you have problems!


# Packages ----------------------------------------------------------------

library(airway)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments) # for single core alingment 
library(BiocParallel)      #multiple core alignment windows can't do it!
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu) # for poisson dist
library(ggplot2)
library(genefilter) # for rowVars
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ReportingTools) #creates html file
library(Gviz)
library(sva)
library(fission)
#Notes----------------------------------------------------------------------------------------------------------
citation("pkgName") # to see the information about how to cite the package

#to install bioconductor packages:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pkgName")

sessionInfo()
# Locating alignment files ---------------------------------------------------------------
dir <- system.file("extdata", package="airway", mustWork=TRUE)
# In this directory, we find the eight BAM files (and some other files):
list.files(dir)
csvfile <- file.path(dir,"sample_table.csv") #load datail table
(sampleTable <- read.csv(csvfile,row.names=1))
filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000) 
seqinfo(bamfiles[1])

# Defining gene models ----------------------------------------------------

gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
(ebg <- exonsBy(txdb, by="gene"))
# Read counting step ------------------------------------------------------
register(SerialParam())
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )

# SummarizedExperiment ----------------------------------------------------
se
assayNames(se)
head(assay(se), 3)
colSums(assay(se))
rowRanges(se)
str(metadata(rowRanges(se)))
colData(se)
(colData(se) <- DataFrame(sampleTable))

# Starting from SummarizedExperiment --------------------------------------
data("airway")
se <- airway
se$dex <- relevel(se$dex, "untrt")
se$dex
round( colSums(assay(se)) / 1e6, 1 )
dds <- DESeqDataSet(se, design = ~ cell + dex)
# Starting from count matrices --------------------------------------------
countdata <- assay(se)
head(countdata, 3)
coldata <- colData(se)
(ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ cell + dex))
# We will continue with the object generated from the SummarizedExperiment section
# Pre-filtering the dataset -----------------------------------------------
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

# The rlog transformation -------------------------------------------------
?rlog
rld <- rlog(dds, blind=FALSE) 
head(assay(rld), 3)
#log2 vs rlog
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)
#*Sequencing depth correction is done automatically for the rlog method (and for varianceStabilizingTransformation).


# Sample distances --------------------------------------------------------

#assess overall similarity between samples
#the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns
sampleDists <- dist( t( assay(rld) ) ) #euclidian dist
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues") ) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)


# PCA plot ----------------------------------------------------------------

plotPCA(rld, intgroup = c("dex", "cell"))

#manual plot
(data <- plotPCA(rld, intgroup = c( "dex", "cell"), returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data$data, aes(PC1, PC2, color=dex, shape=cell)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
# probably should be ylab = data$labels[1] and x lab = data$labels[2]


# MDS plot ----------------------------------------------------------------
# multidimensional scaling (MDS)
# This is useful when we don't have a matrix of data, but only a matrix of distances
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=dex,shape=cell)) + geom_point(size=3)

# Creating the same plot for the PoissonDistance
mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
mdsPois <- cbind(mdsPoisData, as.data.frame(colData(dds)))
ggplot(mdsPois, aes(X1,X2,color=dex,shape=cell)) + geom_point(size=3)


# Differential expression analysis ----------------------------------------
dds <- DESeq(dds)


# Building the results table ----------------------------------------------

# Calling results without any arguments will extract the estimated log2 fold changes and p values for the last 
(res <- results(dds))
mcols(res, use.names=TRUE)
summary(res)
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)

# If we want to raise the log2 fold change threshold, so that we test for genes that show more substantial changes due to treatment, we simply supply a value on the log2 scale.

resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)



# Other comparisons -------------------------------------------------------

# In general, the results for a comparison of any two levels of a variable can be extracted using the contrast argument to results.
results(dds, contrast=c("cell", "N061011", "N61311"))


# Multiple testing --------------------------------------------------------
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
# #DESeq2 uses the Benjamini-Hochberg (BH) adjustment
# We subset the results table to these genes and then sort it by the log2 fold change estimate to get the significant genes with the strongest down-regulation
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
# ... and with the strongest up-regulation:
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])


# Plotting results --------------------------------------------------------
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("dex"))

# We can also make custom plots using the ggplot function from the ggplot2 package

data <- plotCounts(dds, gene=topGene, intgroup=c("dex","cell"), returnData=TRUE)
data
ggplot(data, aes(x=dex, y=count, color=cell)) +
  scale_y_log10() +
  geom_point(position=position_jitter(width=.1,height=0), size=3)

ggplot(data, aes(x=dex, y=count, fill=dex)) +
  scale_y_log10() +
  geom_dotplot(binaxis="y", stackdir="center")

ggplot(data, aes(x=dex, y=count, color=cell, group=cell)) +
  scale_y_log10() + geom_point(size=3) + geom_line()

# An MA-plot provides a useful overview for an experiment with a two-group comparison
plotMA(res, ylim=c(-5,5))
# We can also make an MA-plot for the results table in which we raised the log2 fold change threshold
plotMA(resLFC1, ylim=c(-5,5))
topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

# Another useful diagnostic plot is the histogram of the p values. This plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram.
hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50", border="white")
?hist
as <- as.data.frame(res)


# Gene clustering ---------------------------------------------------------
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("cell","dex")])
pheatmap(mat, annotation_col=df)


# Independent filtering ---------------------------------------------------
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~",round(signif(.5*qs[-1] + .5*qs[-length(qs)],2)))
ratios <- tapply(resLFC1$pvalue, bins, function(p) mean(p < .05, na.rm=TRUE))
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")



# Annotating and exporting results ----------------------------------------

columns(org.Hs.eg.db)

                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

# Now the results have the desired external gene IDs:



# Exporting results -------------------------------------------------------


resOrderedDF <- as.data.frame(resOrdered)[1:100,]
write.csv(resOrderedDF, file="results.csv")

# Another more sophisticated package for exporting results from various Bioconductor analysis packages is the ReportingTools package. 

htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)


# Plotting fold changes in genomic space ----------------------------------

(resGR <- results(dds, lfcThreshold=1, format="GRanges"))
# We need to add the symbol again for labeling the genes on the plot:
resGR$symbol <- mapIds(org.Hs.eg.db, names(resGR), "SYMBOL", "ENSEMBL")


# We create a subset of our full results, for genes within the window We add the gene symbol as a name, if the symbol exists or is not duplicated in our subset.
window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)

# We create a vector specifying if the genes in this subset had a low false discovery rate.
sig <- factor(ifelse(resGRsub$padj < .1 & !is.na(resGRsub$padj),"sig","notsig"))
sig
options(ucscChromosomeNames=FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name="gene ranges", feature=sig)
d <- DataTrack(resGRsub, data="log2FoldChange", baseline=0,
               type="h", name="log2 fold change", strand="+")
plotTracks(list(g,d,a), groupAnnotation="group", notsig="grey", sig="hotpink")


# Removing hidden batch effects -------------------------------------------

# Finally we specify that we want to estimate 2 surrogate variables.
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ dex, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv=2)
svseq$sv
# Because we actually do know the cell lines, we can see how well the SVA method did at recovering these variables
par(mfrow=c(2,1),mar=c(3,5,3,1))
stripchart(svseq$sv[,1] ~ dds$cell,vertical=TRUE,main="SV1")
abline(h=0)
stripchart(svseq$sv[,2] ~ dds$cell,vertical=TRUE,main="SV2")
abline(h=0)

# Finally, in order to use SVA to remove any effect on the counts from our surrogate variables, we simply add these two
# surrogate variables as columns to the DESeqDataSet and then add them to the design:
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + dex
# We could then produce results controlling for surrogate variables by running DESeq with the new design:
ddssva <- DESeq(ddssva)



# Time course experiments -------------------------------------------------

data("fission")
ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)

ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + minute)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),],4)
data <- plotCounts(ddsTC, which.min(resTC$padj),
                   intgroup=c("minute","strain"), returnData=TRUE)
ggplot(data, aes(x=minute, y=count, color=strain, group=strain)) +
  geom_point() + stat_smooth(se=FALSE,method="loess") + scale_y_log10()
resultsNames(ddsTC)
res30 <- results(ddsTC, name="strainmut.minute30", test="Wald")
res30[which.min(resTC$padj),]    

# We can furthermore cluster significant genes by their profiles. We extract a matrix of the shrunken log2 fold changes using the coef function:
betas <- coef(ddsTC)
colnames(betas)

library("pheatmap")
topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)


