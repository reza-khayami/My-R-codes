
# #load pakages -----------------------------------------------------------
getwd()

setRepositories() # 1 2
install.packages("DESeq2")
browseVignettes("DESeq2")
require(data.table)
require(DESeq2)
require(ggplot2)
require(limma)
require(gplots)
require(reshape2)
require(dplyr)



# ###A brief example ------------------------------------------------------
#making samples
setwd("C:/Users/RK1994/Documents/R.working/royan/Advanced1/RNAseq/data")
mypath <- "C:/Users/RK1994/Documents/R.working/royan/Advanced1/RNAseq/data"
A <- function(x){
  z = (matrix(round(runif(100, 0, 10^5)), nrow = 100))
rownames(z) <- paste0("Gene", 1:100)
z}

for (i in 1:6) {
  write.table(A(x), file = paste0(mypath,
        paste0("Sample", i, ".txt")),
        col.names = F,
        quote = F,
        sep ="\t")
}

#makes 6 txt files with 100 numbers
# its better to use file.path instead of paste

##load multple data sets and tabulation

files <- list.files("." , "*.txt") 
# "." means here and "*.txt" mean all txt file in the location

cn <- lapply(files, read.delim, header = F,
             stringsAsFactors = F)
cnDT <- lapply(files, fread)
#cnDT is the same as cn with only one difference: data.table
head(cn[[1]])

#what if we wanted to put all these files togheter?

cn <- do.call(cbind, cn)
cnDT <- do.call(cbind, cnDT)
#puts the files in a table together

head(cn)
# head(cnDT)

#multiple gene1:100 in rows
rownames(cn) <- cn[,1]

cn <- cn[, -seq(1, ncol(cn), 2)]
cnDT <- cnDT[, -seq(3, ncol(cnDT), 2), with = F]

colnames(cn) <- sub(".txt", "", files) 
# replace .txt with nothing in files
colnames(cnDT) <- c("Genes", sub(".txt", "", files))
                        



####Analysis with DEseq2 -------------------------------------------------
# dds <- DESeqDataSetFromMatrix(
#countData = cts,
#colData = coldata,
#design= ~ batch + condition)


# ###load files -----------------------------------------------------------


setwd("C:/Users/RK1994/Documents/R.working/royan/Advanced1/RNAseq/data")

files <- list.files("." , "*.count")

cn <- lapply(files, read.delim,
                 header = F, stringsAsFactors = F,
                 comment.char = "_")

#htseqwrites summary of data at the end of it and we dont't
#want it 
#the comments start with _ so we use comment.char= "_" to
#note them as comments 

cn <- do.call(cbind, cn)
colnames(cn) <- sub("count", "", files) 
rownames(cn) <- cn[,1]
cn <- cn[, -seq(1, ncol(cn), 2)]
colnames(cn) <- sub(".count", "", files) 
colSums(cn)
#Samples are too different so we need to normalize and 
#DEseq2 does that


# ###Desing Matrix --------------------------------------------------------


gr <- factor(c(rep("YRI", 3), rep("GBR", 2), "YRI",
               rep("GBR", 4), rep("YRI",2)))
# its better to be factor
gender <- factor(c(rep("M", 2), rep("F", 2),
                   "M", "F", "F", "M", "M", "F", "M", "F"   ))
# need to be the real sample names but we couldn't find them
coldata <- data.frame(group = gr,
                      type = "paired-end", gender = gender )

cds <- DESeqDataSetFromMatrix(cn, coldata, design = ~group)


# ### Data Normalization --------------------------------------------------

cds <- DESeq(cds)
cnt <- log2(1+counts(cds, normalized = T))

# returns number of counts
min(cnt) # 0 

boxplot(cnt)
#PCA could also be used here
#if rows are samples and genes are columns :
#PCA get the difference between the samples the other way around
#is for gene difference
# x: contains PCs for drwaing a graph
a <- prcomp(t(cnt))

# ### finding DE genes ----------------------------------------------------

dif <- results(cds, c("group", "YRI", "GBR"))
head(dif)
dif <- data.frame(results(cds, c("group", "YRI", "GBR")))
head(dif)

dif$padj <- p.adjust(dif$pvalue, method = "BH")
# we can adjust pvalue our selves to remove NAs

dif <- dif[order(dif$padj), ]
# sorts the table based on adj.pv
head(dif)

#Volcano plot
ggplot(dif, aes(log2FoldChange, -log10(padj), color = padj)) +
  geom_point() + theme_bw()
#instead of padj pvalue could be used
?rpois
