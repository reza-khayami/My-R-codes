# read if you have problems!
# doi: 10.12688/f1000research.7035.1 

# Locating alignment files ---------------------------------------------------------------

library(airway)

# The R function system.file can be used to find out
# where on your computer the files from a package have 
# been installed.
dir <- system.file("extdata", package="airway",
                   mustWork=TRUE)


# In this directory, we find the eight BAM files 
# (and some other files):
list.files(dir)


# Typically, we have a table with detailed information
#for each of our samples that links samples to the
#associated FASTQ and BAM files. For your own project,
#you might create such a comma-separated value (CSV) file
#using a text editor or spreadsheet software such as Excel


#load detail table
csvfile <- file.path(dir,"sample_table.csv") 
(sampleTable <- read.csv(csvfile,row.names=1))
# the parentheses () around the entire code of the
# last line above is an R trick to print the output of 
# the function in addition to saving it to sampleTable

# The following tools can be used generate count matrices:
#summarizeOverlaps, featureCounts, or htseq-count


# We now proceed using the summarizeOverlaps method of
# counting.

# construct the full paths to the files

filenames <- file.path(dir, paste0(sampleTable$Run,
                                   "_subset.bam"))
file.exists(filenames)

library(Rsamtools)

bamfiles <- BamFileList(filenames, yieldSize=2000000) 
# only process 2 million reads at a time

#We need to check if the chromosome names are matched
#between bam and gtf (here called "seqnames")
#in the alignment files
seqinfo(bamfiles[1])


# Defining gene models ----------------------------------------------------
# will read the gene model from an Ensembl GTF file
# TxDb object : a database that can be used to generate
# a variety of range-based objects, such as exons,
# transcripts, and genes

# There are other options for constructing a TxDb.
# For the known genes track from the UCSC Genome Browser

# one can use the pre-built Transcript DataBase:
# TxDb.Hsapiens.UCSC.hg19.knownGene. 

# If the annotation file is accessible from AnnotationHub
# (as is the case for the Ensembl genes),
# a pre-scanned GTF file can be imported using
# makeTxDbFromGRanges. 

# Finally, the makeTxDbFromBiomart function can be used 
#to automatically pull a gene

# model from Biomart using biomaRt15.

library("GenomicFeatures")

# We indicate that none of our sequences (chromosomes) 
# are circular using a 0-length character vector

gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))

# The following line produces a GRangesList of all the
# exons grouped by gene. Each element of the list is 
# a GRanges object of the exons for a gene.
(ebg <- exonsBy(txdb, by="gene"))


# Read counting step ------------------------------------------------------

library(GenomicAlignments) 
library(BiocParallel) #multiple core alignment windows
#can't do it!
register(SerialParam()) #single

options(MulticoreParam=quote(MulticoreParam(workers=4)))
register(MulticoreParam())
param <- MulticoreParam(workers = 4)
# The following call creates the SummarizedExperiment 
# object with counts:
?summarizeOverlaps
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE,
                        BPPARAM = param)

# The mode argument describes what kind of read overlaps
# will be counted

# Setting singleEnd to FALSE indicates that the experiment
# produced paired-end reads, and we want to count a pair
# of reads (a fragment) only once toward the count for
# a gene.

# ignore.strand to TRUE
# The fragments argument can be used when singleEnd=FALSE to specify if unpaired reads should be counted (yes if fragments=TRUE).

??GenomicAlignments # see Counting reads with
#                     summarizeOverlaps
# vignette


# SummarizedExperiment ----------------------------------------------------

se
assayNames(se)
head(assay(se), 3)
colSums(assay(se))

# The rowRanges, when printed, only shows the first
# GRanges, and tells us there are 19 more elements:
rowRanges(se)

# The rowRanges also contains metadata about the
#construction of the gene model in the metadata slot.
#Here we use a helpful R function, str, to display the
#metadata compactly:

  str(metadata(rowRanges(se)))
  colData(se)  
  # The colData slot, so far empty, should contain all
  # the metadata. Because we used a column of sampleTable
  # to produce the bamfiles vector, we know the columns of
  # se are in the same order as the rows of sampleTable.
  # We can assign the sampleTable as the colData of the summarized experiment, by converting it into a DataFrame and using the assignment function:
    (colData(se) <- DataFrame(sampleTable))

# The DESeqDataSet, sample information, and the design formula ------------

# The simplest design formula for differential expression
# would be ~ condition, where condition is a column in 
# colData(dds) that specifies which of two (or more
# groups) the samples belong to. 
# For the airway experiment, we will specify ~ cell + dex
# meaning that we want to test for the effect of 
# dexamethasone (dex) controlling for the effect of 
# different cell line (cell). 
# We can see each of the columns just using the
# $ directly on the SummarizedExperiment or DESeqDataSet:

  
se$cell
se$dex
se$
# it is prefered in R that the first level of a factor
# be the reference level (e.g. control, or untreated
# samples), so we can relevel the dex factor like so:

se$dex <- relevel(se$dex, "untrt")
se$dex
#For running DESeq2 models, you can use R's formula notation to express any fixed-effects experimental design
# Note that DESeq2 uses the same formula notation as, for instance, the lm function of base R.
??results
# If the research aim is to determine for which genes the effect of treatment is different across groups, then interaction terms can be included and tested using a design such as ~ group + treatment +
#   group:treatment. 
# See the manual page for ?results for more examples

# For a full example of using the HTSeq Python package for read counting, please see the pasilla vignette.
# For an example of generating the DESeqDataSet from files produced by htseq-count, please see the DESeq2 vignette.

# In the following sections, we will demonstrate the construction of the DESeqDataSet from two starting points:
# . from a SummarizedExperiment object
# . from a count matrix and a sample information table


# Starting from SummarizedExperiment --------------------------------------

# We now use R's data command to load a prepared SummarizedExperiment that was generated from the publicly available sequencing data files associated with the Himes et al.1 paper, described above. 
# The steps we used to produce thisobject were equivalent to those you worked through in the previous sections, except that we used all the reads and all the genes. 
# For more details on the exact steps used to create this object, type vignette("airway") into your R session.

data("airway")
se <- airway

# Again, we want to specify that untrt is the reference level for the dex variable:
se$dex <- relevel(se$dex, "untrt")
se$dex

# We can quickly check the millions of fragments that uniquely aligned to the genes (the second argument of round tells how many decimal points to keep).
round( colSums(assay(se)) / 1e6, 1 )

# Supposing we have constructed a SummarizedExperiment using one of the methods described in the previous section, we now need to make sure that the object contains all the necessary information about the samples, i.e., a table with metadata on the count matrix's columns stored in the colData slot

colData(se)

# Once we have our fully annotated SummarizedExperiment object, we can construct a DESeqDataSet object from it that will then form the starting point of the analysis.
# We add an appropriate design for the analysis:
  library("DESeq2")
dds <- DESeqDataSet(se, design = ~ cell + dex)

# Starting from count matrices --------------------------------------------
# 
# In this section, we will show how to build an DESeqDataSet supposing we only have a count matrix and a table of sample information

# While the previous section would be used to construct a DESeqDataSet from a SummarizedExperiment, here we first extract the individual object (count matrix and sample info) from the SummarizedExperiment in order to build it back up into a new object - only for demonstration purposes.

# In practice, the count matrix would either be read in from a file or perhaps generated by an R function like featureCounts from the Rsubread package12.
# to see theactual data, i.e., here, the fragment counts, we use the assay function.

countdata <- assay(se)
head(countdata, 3)
coldata <- colData(se)

# We now have all the ingredients to prepare our data object in a form that is suitable for analysis, namely:
#   . countdata: a table with the fragment counts
# . coldata: a table with information about the samples
(ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ cell + dex))

# We will continue with the object generated from the SummarizedExperiment section --------






# Exploratory analysis and visualization ----------------------------------

# There are two separate paths in this workflow; the one we will see first involves transformations of the counts in order to visually explore sample relationships. 
# In the second part, we will go back to the original raw counts for statisticaltesting. This is critical because the statistical testing methods rely on original count data (not scaled or transformed) for calculating the precision of measurements


# Pre-filtering the dataset -----------------------------------------------

# Our count matrix with our DESeqDataSet contains many rows with only zeros, and additionally many rows with only a few fragments total. 
# In order to reduce the size of the object, and to increase the speed of our functions, we can remove the rows that have no or nearly no information about the amount of gene expression. 
# Here we remove rows of the DESeqDataSet that have no counts, or only a single count across all samples:
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

# The rlog transformation -------------------------------------------------

# the rlog transformation is provided for applications other than differential testing.
# For differential testing we recommend the DESeq function applied to raw counts, as described later in this workflow, which also takes into account the dependence of the variance of counts on the mean value during the dispersion estimation step.
# The function rlog returns a SummarizedExperiment object that contains the rlog-transformed values in its assay slot.
?rlog
rld <- rlog(dds, blind=FALSE) 
# blind T for  QA (quality assurance)
# blind F for downstream analysis
head(assay(rld), 3)

# We specify blind=FALSE, which means that differences between cell lines and treatment should not add to thevariance-mean profile of the experiment. 
# However, the experimental design is not used directly in the transformation,only in estimating the global amount of variability in the counts. For a fully unsupervised transformation, one can set
# blind=TRUE (which is the default).

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
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues") ) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Another option is poisson dist
#takes the inherent variance structure of counts into consideration
# The PoissonDistance function takes the original count matrix
# (not normalized) with samples as rows instead of columns, so we need to transpose the counts in dds
install.packages("PoiClaClu")
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)


# PCA plot ----------------------------------------------------------------

# The percent of the total variance that is contained in the direction is printed in the axis label.
# Note that these percentages do not add to 100%, because there are more dimensions that contain the remaining variance (although each of these remaining dimensions will explain less than the two that we see).

plotPCA(rld, intgroup = c("dex", "cell"))

#manual plot
(data <- plotPCA(rld, intgroup = c( "dex", "cell"), returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))

library("ggplot2")
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

# Running the differential expression pipeline
# we can run the differential expression pipeline on the raw counts with a single call to the function DESeq
dds <- DESeq(dds)


# Building the results table ----------------------------------------------

# Calling results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula. If there are more than 2 levels for this variable, results will extract the results table for a comparison of the last level over the first level. 
# This comparison is printed at the top of the output: dex trt vs untrt

(res <- results(dds))

# As res is a DataFrame object, it carries metadata with information on the meaning of the columns:
mcols(res, use.names=TRUE)
# The first column, baseMean, is a just the average of the normalized count values, dividing by size factors, taken over all samples in the DESeqDataSet
# The column log2FoldChange is the effect size estimate. It tells us how much the gene's expression seems to have changed due to treatment with dexamethasone in comparison to untreated samples

# DESeq2 performs for each gene a hypothesis test to see whether evidence is sufficient to decide against the null hypothesis that there is zero effect of the treatment on the gene and that the observed difference between treatment and control was merely caused by experimental variability
summary(res)

# there are two ways to be more strict about which set of genes are considered significant:
#   . lower the false discovery rate threshold (the threshold on padj in the results table)
# . raise the log2 fold change threshold from 0 using the lfcThreshold argument of results

# If we lower the false discovery rate threshold, we should also tell this value to results(), so that the function will use an alternative threshold for the optimal independent filtering step:

res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)

# If we want to raise the log2 fold change threshold, so that we test for genes that show more substantial changes due to treatment, we simply supply a value on the log2 scale.

resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

# Sometimes a subset of the p values in res will be NA ("not available"). This is DESeq's way of reporting that all counts for this gene were zero, and hence no test was applied. In addition, p values can be assigned NA if the gene was excluded from analysis because it contained an extreme count outlier

citation("pkgName") # to see the information about how to cite the package


# Other comparisons -------------------------------------------------------

# In general, the results for a comparison of any two levels of a variable can be extracted using the contrast argument to results.
# The user should specify three values:
# the name of the variable,
# the name of the level for the numerator
# the name of the level for the denominator
# Here we extract results for the log2 of the fold change of one cell line over another:
results(dds, contrast=c("cell", "N061011", "N61311"))


# Multiple testing --------------------------------------------------------
# In high-throughput biology, we are careful to not use the p values directly as evidence against the null, but to correct for multiple testing.
# What would happen if we were to simply threshold the p values at a low value, say 0.05? There are 5722 genes with a p value below 0.05 among the 29391 genes, for which the test succeeded in reporting a p value:

sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
# #DESeq2 uses the Benjamini-Hochberg (BH) adjustment
# We subset the results table to these genes and then sort it by the log2 fold change estimate to get the significant genes with the strongest down-regulation
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
# ... and with the strongest up-regulation:
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])


# Plotting results --------------------------------------------------------
# 
# A quick way to visualize the counts for a particular gene is to use the plotCounts function that takes as arguments the DESeqDataSet, a gene name, and the group over which to plot the counts

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

#Since the clustering is only relevant for genes that actually carry a signal, one usually would only cluster a subset of the most highly variable genes
# for demonstration, let us select the 20 genes with the highest variance across samples. We will work with the rlog transformed counts:
library("genefilter") # for rowVars
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)

# The heatmap becomes more interesting if we do not look at absolute expression strength but rather at the amount by which each gene deviates in a specific sample from the gene's average across all samples.
# Hence, we center each genes'values across samples, and plot a heatmap
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("cell","dex")])
pheatmap(mat, annotation_col=df)


# Independent filtering ---------------------------------------------------
# 
# we have no chance of seeing differential expression, because the low read counts suffer from such high Poisson noise that any biological
# effect is drowned in the uncertainties from the sampling at a low rate. We can also show this by examining the ratio of small p values (say, less than, 0.05) for genes binned by mean normalized count

# In the following code chunk, we create bins using the quantile function, bin the genes by base mean using cut, rename the levels of the bins using the middle point, calculate the ratio of p values less than 0.05 for each bin, and finally plot these ratios
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~",round(signif(.5*qs[-1] + .5*qs[-length(qs)],2)))
ratios <- tapply(resLFC1$pvalue, bins, function(p) mean(p < .05, na.rm=TRUE))
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")

#At first sight, there may seem to be little benefit in filtering out these genes. After all, the test found them to be nonsignificant
# anyway. However, these genes have an influence on the multiple testing adjustment, whose performance improves if such genes are removed.
# By removing the low count genes from the input to the FDR procedure, we can find more genes to be significant among those that we keep, and so improved the power of our test. 
# This approach is known as independent filtering.
# The DESeq2 software automatically performs independent filtering that maximizes the number of genes with adjusted p value less than a critical value (by default, alpha is set to 0.1). 
# This automatic independent filtering is performed by,and can be controlled by, the results function.


# Annotating and exporting results ----------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
# We can use the mapIds function to add individual columns to our results table. We provide the row names of our results table as a key, and specify that keytype=ENSEMBL.
# The column argument tells the mapIds function which information
# we want, and the multiVals argument tells the function what to do if there are multiple possible values for a single input value. Here we ask to just give us back the first one that occurs in the database.
# To add the gene symboland Entrez ID, we call mapIds twice.
res$symbol <- mapIds(org.Hs.eg.db,
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

# The call to as.data.frame is necessary to convert the DataFrame object (IRanges package) to a data.frame object that can be processed by write.csv. Here, we take just the top 100 genes for demonstration.
resOrderedDF <- as.data.frame(resOrdered)[1:100,]
write.csv(resOrderedDF, file="results.csv")

# Another more sophisticated package for exporting results from various Bioconductor analysis packages is the ReportingTools package. 
# ReportingTools will automatically generate dynamic HTML documents, including links to external databases using gene identifiers and boxplots summarizing the normalized counts across groups.
# The simplest version of creating a dynamic ReportingTools report is performed with the following code

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ReportingTools")

library("ReportingTools")
htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)


# Plotting fold changes in genomic space ----------------------------------

# If we have used the summarizeOverlaps function to count the reads, then our DESeqDataSet object is built on top of ready-to-use Bioconductor objects specifying the genomic ranges of the genes. 
# We can therefore easily plot our differential expression results in genomic space. 
# While the results function by default returns a DataFrame, using the
# format argument, we can ask for GRanges or GRangesList output.

(resGR <- results(dds, lfcThreshold=1, format="GRanges"))
# We need to add the symbol again for labeling the genes on the plot:
  resGR$symbol <- mapIds(org.Hs.eg.db, names(resGR), "SYMBOL", "ENSEMBL")

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("Gviz")
  
  # We will use the Gviz package for plotting the GRanges and associated metadata: the log fold changes due to dexamethasone treatment.
  library("Gviz")
  
  # The following code chunk specifies a window of 1 million base pairs upstream and downstream from the gene with the smallest p value.
  # We create a subset of our full results, for genes within the window We add the gene symbol as a name, if the symbol exists or is not duplicated in our subset.
?strand
    window <- resGR[topGene] + 1e6
  strand(window) <- "*"
  resGRsub <- resGR[resGR %over% window]
  naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
  resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)
  
  # We create a vector specifying if the genes in this subset had a low false discovery rate.
  sig <- factor(ifelse(resGRsub$padj < .1 & !is.na(resGRsub$padj),"sig","notsig"))
  sig
  # We can then plot the results using Gviz functions.
  # We create an axis track specifying our location in the  genome, a track that will show the genes and their names, colored by significance, and a data track that will draw vertical bars showing the moderated log fold change produced by DESeq2, which we know are only large when the effect is well supported by the information in the counts.
  options(ucscChromosomeNames=FALSE)
  g <- GenomeAxisTrack()
  a <- AnnotationTrack(resGRsub, name="gene ranges", feature=sig)
  d <- DataTrack(resGRsub, data="log2FoldChange", baseline=0,
                 type="h", name="log2 fold change", strand="+")
  plotTracks(list(g,d,a), groupAnnotation="group", notsig="grey", sig="hotpink")
  

# Removing hidden batch effects -------------------------------------------

  # Suppose we did not know that there were different cell lines involved in the experiment, only that there was treatmentwith dexamethasone.  
  # The cell line effect on the counts then would represent some hidden and unwanted variation that   might be affecting many or all of the genes in the dataset.
  # We can use statistical methods designed for RNA-seq from  the sva package to detect such groupings of the samples, and then we can add these to the DESeqDataSet design,  in order to account for them.
  # The SVA package uses the term surrogate variables for the estimated variables that we  want to account for in our analysis
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("sva")
  library("sva")
  
  # Below we obtain a matrix of normalized counts for which the average count across samples is larger than 1. 
  # # As we described above, we are trying to recover any hidden batch effects, supposing that we do not know the cell line information.
  # So we use a full model matrix with the dex variable, and a reduced, or null, model matrix with only an intercept  term. 
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

    # DESeq2 can be used to analyze time course experiments, for example to find those genes that react in a conditionspecific manner over time, compared to a set of baseline samples. Here we demonstrate a basic time course analysis    with the fission data package, that contains gene counts for an RNA-seq time course of fission yeast    
    # The yeast were  exposed to oxidative stress, and half of the samples contain a deletion of the gene atf21. We use a design formula
    # that models the strain difference at time 0, the difference over time, and any strain-specific differences over time (the interaction term strain:minute).
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install("fission")
    library("fission")
    data("fission")
    ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)
    # The following chunk of code performs a likelihood ratio test, where we remove the strain-specific differences over  time.
    # Genes with small p values from this test are those which at one or more time points after time 0 showed a strainspecific  effect.
    # Note therefore that this will not give small p values to genes that moved up or down over time in the  same way in both strains.
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
    
    sessionInfo()
    