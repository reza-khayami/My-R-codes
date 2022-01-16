setwd("~/testsurv/")
getwd()
library(Biobase)
library(GEOquery)
BiocManager::install("GEOquery")
# load series and platform data from GEO
gset <- getGEO('GSE2990', GSEMatrix =TRUE, getGPL=FALSE)
x <- exprs(gset[[1]])
# remove Affymetrix control probes
x <- x[-grep('^AFFX', rownames(x)),]
boxplot(x)
# transform the expression data to Z scores
x <- t(scale(t(x)))
# extract information of interest from the phenotype data (pdata)
idx <- which(colnames(pData(gset[[1]])) %in%
               c('age:ch1', 'distant rfs:ch1', 'er:ch1',
                 'ggi:ch1', 'grade:ch1', 'node:ch1',
                 'size:ch1', 'time rfs:ch1'))
metadata <- data.frame(pData(gset[[1]])[,idx],
                       row.names = rownames(pData(gset[[1]])))
# remove samples from the pdata that have any NA value
discard <- apply(metadata, 1, function(x) any( is.na(x) ))
metadata <- metadata[!discard,]
# filter the Z-scores expression data to match the samples in our pdata
x <- x[,which(colnames(x) %in% rownames(metadata))]
# check that sample names match exactly between pdata and Z-scores 
all((colnames(x) == rownames(metadata)) == TRUE)
## [1] TRUE
# create a merged pdata and Z-scores object
coxdata <- data.frame(metadata, t(x))
# tidy column names
colnames(coxdata)[1:8] <- c('Age', 'Distant.RFS', 'ER',
                            'GGI', 'Grade', 'Node',
                            'Size', 'Time.RFS')
# prepare phenotypes
coxdata$Distant.RFS <- as.numeric(coxdata$Distant.RFS)
coxdata$Time.RFS <- as.numeric(gsub('^KJX|^KJ', '', coxdata$Time.RFS))
coxdata$ER <- factor(coxdata$ER, levels = c(0, 1))
coxdata$Grade <- factor(coxdata$Grade, levels = c(1, 2, 3))

# With the data prepared, we can now apply a Cox survival model independently for each gene (probe) in the dataset against RFS.
# Here we will use RegParallel to fit the Cox model independently for each gene.

library(survival)
BiocManager::install("RegParallel")
library(RegParallel)
res <- RegParallel(
  data = coxdata,
  formula = 'Surv(Time.RFS, Distant.RFS) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(coxdata)[9:ncol(coxdata)],
  blocksize = 2000,
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95)
res <- res[!is.na(res$P),]
res
# 3, annotate top hits with biomaRt
# 
# Filter by Log Rank p<0.01
res <- res[order(res$LogRank, decreasing = FALSE),]
final <- subset(res, LogRank < 0.01)
probes <- gsub('^X', '', final$Variable)
BiocManager::install("biomaRt")
library(biomaRt)
mart <- useMart('ENSEMBL_MART_ENSEMBL', host='useast.ensembl.org')
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(mart = mart,
                     attributes = c('affy_hg_u133a',
                                    'ensembl_gene_id',
                                    'gene_biotype',
                                    'external_gene_name'),
                     filter = 'affy_hg_u133a',
                     values = probes,
                     uniqueRows = TRUE)

# 4, encode statistically significant genes as Low | Mid | High and plot survival curves
# extract RFS and probe data for downstream analysis
survplotdata <- coxdata[,c('Time.RFS', 'Distant.RFS',
                           'X203666_at', 'X205680_at')]
colnames(survplotdata) <- c('Time.RFS', 'Distant.RFS',
                            'CXCL12', 'MMP10')
# set Z-scale cut-offs for high and low expression
highExpr <- 1.0
lowExpr <- -1.0
survplotdata$CXCL12 <- ifelse(survplotdata$CXCL12 >= highExpr, 'High',
                              ifelse(survplotdata$CXCL12 <= lowExpr, 'Low', 'Mid'))
survplotdata$MMP10 <- ifelse(survplotdata$MMP10 >= highExpr, 'High',
                             ifelse(survplotdata$MMP10 <= lowExpr, 'Low', 'Mid'))

# relevel the factors to have mid as the ref level
survplotdata$CXCL12 <- factor(survplotdata$CXCL12,
                              levels = c('Mid', 'Low', 'High'))
survplotdata$MMP10 <- factor(survplotdata$MMP10,
                             levels = c('Mid', 'Low', 'High'))
install.packages("survminer")
library(survminer)

ggsurvplot(survfit(Surv(Time.RFS, Distant.RFS) ~ CXCL12,
                   data = survplotdata),
           data = survplotdata,
           risk.table = TRUE,
           pval = TRUE,
           break.time.by = 500,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE)
