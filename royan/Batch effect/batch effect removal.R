##########################Batch effect#################################
#LOAD PACKAGES----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(GEOquery)
library(limma)
library(ggplot2)
library(sva)
library(VennDiagram)
library(parallel)
library(pheatmap)
library(DESeq2)
#########you can download full table or use annotation table downloaded with GEOquery

#note: if you have lot of data to read: -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ds.plat <- c("GPL570" , "GPL6244")
# ds.id <- c("GSE21222", "GSE46872")
#qc = F
# data.files <- paste0(ds.id, "_series_matrix.txt")
# raw <- lapply(data.files, fread, skip = "!platform_table_begin")
# raw <- lapply(raw, function(x) {rownames(x) <- x[,1];x[,-1]})
# data <- lapply(row, function(x) if(max(x) > 100) log2(x+1) else x)
#if (qc) {  #change qc to T if you want to draw boxplot
# pdf("Results/boxplot.pdf", width = 20)
# invisible(lapply(data, boxplot))
# dev.off()
#}
# # plat.files <- paste0(ds.plat, ".annot")
# ann <- lapply(paste0("", plat.files), fread, skip = "!platform_table_begin" )
# ann[[1]]
# ann[[1]] <- ann[[1]][, c("ID", "Gene symbol", "Gene ID")]
# ann[[1]] <- data.frame(sapply(ann[[1]], function(x) sub("///.*", ",", x)))
# annNames <- c("ID", "symbol", "Gene_ID")
# colnames(ann[[1]]) <- annNames
# do the same for other data sets you need to do the dirtywork of cleaning the datasets!
#data <- lapply(data, function(x) data.frame(ID = rownames(x), x))
#data <- lapply(seq_along(data), function(i) merge( ann[[i]], data[[i]], by = "ID"))
#you can delete bad data sets!
#.....     data <- lapply(data, function(x) aggregate(. ~ symbol, x, mean)  
# mclapply if multi core data <- mclapply(data, function(x) aggregate(. ~ symbol, x, mean), mc.cores =   detectCores())
#mdata <- Reduce(function(...) merge(..., by = "symbol"), data)
#rownames(mdata) <- mdata$symbol
#delet symbol column
#batch <- rep(seq_along(data), times = sapply(data, ncol)-1 )
#mdata.c <- Combat(mdata, batch)
#save(mdata.c, mdata, hgnc, file = "Datasets.Rdata")
#setwd()
#load("Datasets.Rdata")
##if (qc) {  #change qc to T if you want to draw boxplot
# pdf("Results/BoxplotAll.pdf", width = 40)
# boxplot(mdata)
#boxplot(mdata.q)
#boxplot(mdata.c)
# dev.off()}
#make two datasets with factor for samples e.g. (naive , primed) and authors and year of publish
#pca
#ANNOTATION WITH GEOquery-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#_____________________________LOAD DATA_______________________________#

##first read the data and save all of them in one folder expression matrixes

# datasets used :GSE21222, GSE46872

getwd()
setwd("C:/Users/RK1994/Documents/R.working/royan/Advanced1/Batch effect/DATA")
data1 <- read.delim("GSE21222_series_matrix.txt",comment.char = "!") 
data2 <- read.delim("GSE46872_series_matrix.txt",comment.char = "!")

max(data1[,-1])
max(data2[,-1])
min(data1[,-1])
#data1 is not log scale
data1[,-1] <- log2(1+data1[,-1])


#or annotation get annotation with GEOquery
setwd("C:/Users/RK1994/Documents/R.working/royan/Advanced1/Batch effect/Ann/")
annot1 <- fread("GPL570.annot", skip = "!platform_table_begin", data.table = F)
annot2<- fread("GPL6244.annot", skip = "!platform_table_begin", data.table = F)
#check tables
head(annot1)
head(annot2)
dim(data1)
dim(data2)
dim(annot2)
colnames(annot1)
colnames(annot2)




#_________________________PREPARING THE TABLES_________________________#

head(annot2$`Gene ID`)
#rename symbols so that only one entry remains for each symbols
annot1$`Gene symbol` <- sub("///.*","", annot1$`Gene symbol`) 
annot1$`Gene ID` <- sub("///.*","",annot1$`Gene ID`)
annot2$`Gene symbol` <- sub("///.*","", annot2$`Gene symbol`) 
annot2$`Gene ID` <- sub("///.*","",annot2$`Gene ID`)

#rename column names
annot1 <- annot1[, c("ID", "Gene symbol", "Gene ID")]
colnames(annot1) <- c("ID", "Symbol", "Gene_ID")
head(annot1)

annot2 <- annot2[, c("ID", "Gene symbol", "Gene ID")]
colnames(annot2) <- c("ID", "Symbol", "Gene_ID")
head(annot2)

# we want to use one of the columns which has  more common entries 
length(intersect(annot1$Symbol, annot2$Symbol))
length(intersect(annot1$Gene_ID, annot2$Gene_ID))

dim(data2)
dim(annot2)
#data2 and annot2 have differents numbers of rows

#setting IDs as rownames
rownames(annot1) <- annot1$ID
rownames(annot2) <- annot2$ID
head(annot1$ID)
annot1["117_at",]

#assigning Gene IDs in a column named Entrez in each data file from annotations which have the same entry as the data
data1$Entrez <- annot1[as.character(data1$ID_REF), "Gene_ID"]
data2$Entrez <- annot2[as.character(data2$ID_REF), "Gene_ID"]
#checking validity
sum(is.na(data1$Entrez))
sum(is.na(data2$Entrez))
dim(data2)
annot2[data2$ID_REF,]
dim(annot2)
data1[85, "Entrez"]
annot1["1552368_at",]
# removing ID_REF
data1 <- data1[,-1]
data2 <- data2[,-1]

#___________________________MERGING THE TABLES_________________________#

#dealing with duplicates
data1 <- aggregate(. ~ Entrez, data1, FUN = mean)
data2 <- aggregate(. ~ Entrez, data2, FUN = mean)

all <- merge(data1, data2, "Entrez")
dim(all)
#Removing NAs
all <- subset(all, Entrez != "")
boxplot(all[,-1])



#D
#ANNOTATION WITH DOWLOADED FULL TABLE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#_____________________________LOAD DATA_______________________________#

##first read the data and save all of them in one folder expression matrixes

# datasets used :GSE21222, GSE46872

getwd()
setwd("C:/Users/RK1994/Documents/R.working/royan/Advanced1/Batch effect/DATA/")
data1 <- read.delim("GSE21222_series_matrix.txt",comment.char = "!") 
data2 <- read.delim("GSE46872_series_matrix.txt",comment.char = "!")

max(data1[,-1])
max(data2[,-1])
min(data1[,-1])
#data1 is not log scale
data1[,-1] <- log2(1+data1[,-1])

#Download full table for platform
setwd("C:/Users/RK1994/Documents/R.working/royan/Advanced1/Batch effect/Ann/")
annot1 <- fread("GPL570-55999.txt",header = F, data.table = F)
#if header = T, encounters error dont know why!
colnames(annot1) <- annot1[1,]
annot1 <- annot1[-1,]
annot2 <- fread("GPL6244-17930.txt", data.table = F)


#check tables
head(annot1)
head(annot2)
dim(data1)
dim(data2)
dim(annot2)
colnames(annot1)
colnames(annot2)


# if you have lot of data to read: 
# ds.plat <- c("GPL570" , "GPL6244")
# ds.id <- c("GSE21222", "GSE46872")
# data.files <- paste0(ds.id, "_series_matrix.txt")
# raw <- lapply(data.files, fread, skip = "!platform_table_begin")
# raw <- lapply(raw, function(x) {rownames(x) <- x[,1];x[,-1]})
# data <- lapply(row, function(x) if(max(x) > 100) log2(x+1) else x)
# # plat.files <- paste0(ds.plat, ".annot")
# ann <- lapply(paste0("", plat.files), fread, skip = "!platform_table_begin" )
# ann[[1]]
# ann[[1]] <- ann[[1]][, c("ID", "Gene symbol", "Gene ID")]
# ann[[1]] <- data.frame(sapply(ann[[1]], function(x) sub("///.*", ",", x)))
# annNames <- c("ID", "symbol", "Gene_ID")
# colnames(ann[[1]]) <- annNames


#_________________________PREPARING THE TABLES_________________________#
head(annot2$gene_assignment)
head(annot1$`Gene Symbol`)

#rename symbols so that only one entry remains for each symbols
annot1$`Gene Symbol` <- sub(" .*","", annot1$`Gene Symbol`) 
annot1$ENTREZ_GENE_ID <- sub(" .*","",annot1$ENTREZ_GENE_ID)
#. means every letter next to it and * means every repeat of it
colnames(annot1)

strsplit2("adad // adsdad // dadadda // ", " // ")
s <- strsplit2(annot2$gene_assignment, " /* ") 

# /* means / and every duplicate of /
head(s[4,])

#rename column names
annot1 <- annot1[, c("ID", "Gene Symbol", "ENTREZ_GENE_ID")]
colnames(annot1) <- c("ID", "Symbol", "ENTREZ")
head(annot1)
head(s[2, ])
dim(s)
annot2$Symbol <- s[,2]
annot2$ENTREZ <- s[,5]
annot2 <- annot2[, c("ID", "Symbol", "ENTREZ")]
head(annot2)
# we want to use one of the columns which has  more common entries 
length(intersect(annot1$Symbol, annot2$Symbol)) #18037
length(intersect(annot1$ENTREZ, annot2$ENTREZ)) #18072
length(unique(annot1$ENTREZ)) 
length(unique(annot2$ENTREZ)) 
# we'll lose approximatly 10 precent of data ~=2000 not very bad!


dim(data2)
dim(annot2)
#data2 and annot2 have differents numbers of rows

#setting IDs as rownames
rownames(annot1) <- annot1$ID
rownames(annot2) <- annot2$ID
head(annot1$ID)
annot1["117_at",]

#assigning Gene IDs in a column named Entrez in each data file from annotations which have the same entry as the data
data1$ENTREZ <- annot1[as.character(data1$ID_REF), "ENTREZ"]
data2$ENTREZ <- annot2[as.character(data2$ID_REF), "ENTREZ"]
length(unique(data2$ENTREZ))
#checking validity
sum(is.na(data1$Entrez))
sum(is.na(data2$Entrez))
dim(data2)
annot2[data2$ID_REF,]
dim(annot2)
data1[85, "ENTREZ"]
annot1["1552368_at",]
# removing ID_REF
data1 <- data1[,-1]
data2 <- data2[,-1]

#___________________________MERGING THE TABLES_________________________#

#dealing with duplicates
data1 <- aggregate(. ~ ENTREZ, data1, FUN = mean)
data2 <- aggregate(. ~ ENTREZ, data2, FUN = mean)

all <- merge(data1, data2, "ENTREZ")
dim(all)
#tred lightly with merge!
#Removing NAs
all <- subset(all, ENTREZ != "")
rownames(all) <- all$ENTREZ
all <- subset(all, select = -ENTREZ)
max(all) #make sure all data are log scaled
boxplot(all)

#QUALITY CONTROL--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#make a file to define which sample is from which data

batch <- factor(c(rep(1, ncol(data1)-1), rep(2, ncol(data2)-1)))
gr <- c(rep("Primed", 12), rep("Naive", 10), "Primed", "Naive",
        "Primed", rep("Naive", 4), "Primed")

#_______________________________PCA____________________________________# 
#if rows are samples and genes are columns : PCA get the difference between the samples the other way around is for gene difference
# x: contains PCs for drwaing a graph

pc <- prcomp(t(all), scale = TRUE)
pcx <- data.frame(pc$x[, 1:3], batch, gr)

pca.var <- pc$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

ggplot(data=pcx, aes(x=PC1, y=PC2, color = batch, shape = gr)) + geom_point() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() + 
  ggtitle("My PCA Graph")


#or
all.m <- all - rowMeans(all)
pcm <- prcomp(all.m)
pcmr <- data.frame(pcm$r[, 1:3], batch, gr)

pca.var <- pcm$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

ggplot(data=pcmr, aes(x=PC1, y=PC2, color = gr, shape = batch)) + geom_point() + 
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() + 
  ggtitle("My PCA Graph")

#NORMALIZATION----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#_____________________QUANTILE NORMALIZATION___________________________#

all.q <- normalizeQuantiles(all)
boxplot(all.q)

all.m <- all.q - rowMeans(all.q)
pcm <- prcomp(all.m)
pcmr <- data.frame(pcm$r[, 1:3], batch, gr)

pca.var <- pcm$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

ggplot(data=pcmr, aes(x=PC1, y=PC2, color = gr, shape = batch)) + geom_point() + 
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() + 
  ggtitle("My PCA Graph")
#it's not good
#_____________________SURROGATE VARIABLE ANALYSIS______________________#
library(sva)
allc <- ComBat(as.matrix(all), batch = batch)
boxplot(allc)
all.m <- allc - rowMeans(allc)
pcm <- prcomp(all.m)
pcmr <- data.frame(pcm$r[, 1:3], batch, gr)

pca.var <- pcm$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

ggplot(data=pcmr, aes(x=PC1, y=PC2, color = gr, shape = batch)) + geom_point() + 
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() + 
  ggtitle("My PCA Graph")