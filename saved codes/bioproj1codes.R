

####### 1. Install Packages-----------
setwd("/media/rk94/01CD77B76BF0B4F0/R")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)


install.packages("pheatmap")

BiocManager::install("Rgraphviz")
BiocManager::install("graph")
install.packages("CePa")

install.packages("survminer")




###### 2. Load Packages------------
package_list <- c("TCGAbiolinks", "ggplot2", "dplyr", "DESeq2",
                  "estimate", "CePa", "BiocParallel",
                  "survminer", "pheatmap", "RColorBrewer", "gplots",
                  'RegParallel')

lapply(package_list, library, character.only = TRUE)


####### Multi Core-----------
options(MulticoreParam=quote(MulticoreParam(workers=4)))
register(MulticoreParam(4))
param <- MulticoreParam(workers = 4)
bplapply(1:4, mean, BPPARAM = param)




###### 3. Search Data------------
patient_list <-  c("TCGA-AB-2803", "TCGA-AB-2805",
                  "TCGA-AB-2806", "TCGA-AB-2807",
                  "TCGA-AB-2808", "TCGA-AB-2810",
                  "TCGA-AB-2811", "TCGA-AB-2812",
                  "TCGA-AB-2813", "TCGA-AB-2814",
                  "TCGA-AB-2815", "TCGA-AB-2816",
                  "TCGA-AB-2817", "TCGA-AB-2818",
                  "TCGA-AB-2819", "TCGA-AB-2820",
                  "TCGA-AB-2821", "TCGA-AB-2822",
                  "TCGA-AB-2823", "TCGA-AB-2824",
                  "TCGA-AB-2825", "TCGA-AB-2826",
                  "TCGA-AB-2828", "TCGA-AB-2830",
                  "TCGA-AB-2832", "TCGA-AB-2833",
                  "TCGA-AB-2834", "TCGA-AB-2835",
                  "TCGA-AB-2836", "TCGA-AB-2837",
                  "TCGA-AB-2838", "TCGA-AB-2839",
                  "TCGA-AB-2840", "TCGA-AB-2841",
                  "TCGA-AB-2842", "TCGA-AB-2843",
                  "TCGA-AB-2844", "TCGA-AB-2845",
                  "TCGA-AB-2846", "TCGA-AB-2847",
                  "TCGA-AB-2848", "TCGA-AB-2849",
                  "TCGA-AB-2851", "TCGA-AB-2853",
                  "TCGA-AB-2854", "TCGA-AB-2855",
                  "TCGA-AB-2856", "TCGA-AB-2857",
                  "TCGA-AB-2858", "TCGA-AB-2859",
                  "TCGA-AB-2860", "TCGA-AB-2861",
                  "TCGA-AB-2862", "TCGA-AB-2863",
                  "TCGA-AB-2865", "TCGA-AB-2866",
                  "TCGA-AB-2867", "TCGA-AB-2868",
                  "TCGA-AB-2869", "TCGA-AB-2870",
                  "TCGA-AB-2871", "TCGA-AB-2872",
                  "TCGA-AB-2873", "TCGA-AB-2874",
                  "TCGA-AB-2875", "TCGA-AB-2877",
                  "TCGA-AB-2879", "TCGA-AB-2880",
                  "TCGA-AB-2881", "TCGA-AB-2882",
                  "TCGA-AB-2884", "TCGA-AB-2885",
                  "TCGA-AB-2886", "TCGA-AB-2887",
                  "TCGA-AB-2888", "TCGA-AB-2889",
                  "TCGA-AB-2890", "TCGA-AB-2891",
                  "TCGA-AB-2895", "TCGA-AB-2896",
                  "TCGA-AB-2897", "TCGA-AB-2898",
                  "TCGA-AB-2899", "TCGA-AB-2900",
                  "TCGA-AB-2901", "TCGA-AB-2903",
                  "TCGA-AB-2904", "TCGA-AB-2908",
                  "TCGA-AB-2909", "TCGA-AB-2910",
                  "TCGA-AB-2911", "TCGA-AB-2912",
                  "TCGA-AB-2913", "TCGA-AB-2914",
                  "TCGA-AB-2915", "TCGA-AB-2916",
                  "TCGA-AB-2917", "TCGA-AB-2918",
                  "TCGA-AB-2919", "TCGA-AB-2920",
                  "TCGA-AB-2921", "TCGA-AB-2924",
                  "TCGA-AB-2925", "TCGA-AB-2927",
                  "TCGA-AB-2928", "TCGA-AB-2929",
                  "TCGA-AB-2930", "TCGA-AB-2931",
                  "TCGA-AB-2932", "TCGA-AB-2933",
                  "TCGA-AB-2934", "TCGA-AB-2935",
                  "TCGA-AB-2936", "TCGA-AB-2937",
                  "TCGA-AB-2938", "TCGA-AB-2939",
                  "TCGA-AB-2940", "TCGA-AB-2941",
                  "TCGA-AB-2942", "TCGA-AB-2943",
                  "TCGA-AB-2944", "TCGA-AB-2946",
                  "TCGA-AB-2948", "TCGA-AB-2949",
                  "TCGA-AB-2950", "TCGA-AB-2952",
                  "TCGA-AB-2954", "TCGA-AB-2955",
                  "TCGA-AB-2956", "TCGA-AB-2959",
                  "TCGA-AB-2963", "TCGA-AB-2964",
                  "TCGA-AB-2965", "TCGA-AB-2966",
                  "TCGA-AB-2967", "TCGA-AB-2969",
                  "TCGA-AB-2970", "TCGA-AB-2971",
                  "TCGA-AB-2972", "TCGA-AB-2973",
                  "TCGA-AB-2975", "TCGA-AB-2976",
                  "TCGA-AB-2977", "TCGA-AB-2978",
                  "TCGA-AB-2979", "TCGA-AB-2980",
                  "TCGA-AB-2981", "TCGA-AB-2982",
                  "TCGA-AB-2983", "TCGA-AB-2984",
                  "TCGA-AB-2985", "TCGA-AB-2986",
                  "TCGA-AB-2987", "TCGA-AB-2988",
                  "TCGA-AB-2990", "TCGA-AB-2991",
                  "TCGA-AB-2992", "TCGA-AB-2993",
                  "TCGA-AB-2994", "TCGA-AB-2995",
                  "TCGA-AB-2996", "TCGA-AB-2998",
                  "TCGA-AB-2999", "TCGA-AB-3000",
                  "TCGA-AB-3001", "TCGA-AB-3002",
                  "TCGA-AB-3005", "TCGA-AB-3006",
                  "TCGA-AB-3007", "TCGA-AB-3008",
                  "TCGA-AB-3009", "TCGA-AB-3011",
                  "TCGA-AB-3012")

query <- GDCquery(project = "TCGA-LAML",
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  legacy = TRUE,
                  platform = "Illumina HiSeq",
                  experimental.strategy	= "RNA-Seq",
                  barcode = patient_list
                  )
# DESeq2 uses unnormalized counts so we download
# unnormalized counts
df <- query[[1]][[1]]
rows <- rownames(df[lapply(df$tags, "[", 1) == "unnormalized", ])
query[[1]][[1]] <- query[[1]][[1]][rows, ]

#clinicl data query
query_clin <- GDCquery(project = "TCGA-LAML",
                       data.category = "Clinical",
                       legacy = TRUE,
                       barcode = patient_list
                       )
#clean up
rm(df, rows)

###### 4. Download Data------------
GDCdownload(query, method = "api")
GDCdownload(query_clin, method = "api")







###### 5. Prepare------------
Data <- GDCprepare(query, summarizedExperiment = TRUE)
clinical <- GDCprepare_clinic(query_clin, "patient")
                    
                    
?





###### 6. Prepare for ESTIMATE------------
#normalize
genetable <- assay(Data)
genetable <- apply(genetable, 1:2, round)
AMLexpr <- vst(genetable+1, blind = FALSE)

# change (gene | number) to (gene) for example:
#A1BG | 1 to A1BG
row <- rownames(AMLexpr)
row <- sub("\\|.*", "", row)
rownames(AMLexpr) <- sub("\\|.*", "", rownames(AMLexpr))

#shorten col names
colnames(AMLexpr) <- noquote(substr(colnames(AMLexpr), 9, 12))


#find out if the rows are duplicated
dups <- duplicated(rownames(AMLexpr))
rownames(AMLexpr)[dups]

# a table of rownames with their assigned ensembl ID
# just to check out

#Get ensembl IDs
aws <- Data@rowRanges@elementMetadata@listData[["ensembl_gene_id"]]
aws[which(rownames(AMLexpr) == "?")]

#make the table
ungenes <- data.frame(
  ensembl = aws[c((which(rownames(AMLexpr) == "?")),
                  (which(rownames(AMLexpr) == "SLC35E2")))],
  rowname = c((which(rownames(AMLexpr) == "?")),
              (which(rownames(AMLexpr) == "SLC35E2"))),
  genes = Data@rowRanges@elementMetadata@listData[["gene_id"]][c((which(rownames(AMLexpr) == "?")),
                                                                 (which(rownames(AMLexpr) == "SLC35E2")))])


View(ungenes)


# change ? and other duplicated names in rows 
# with real gene names
duplist <-c((which(rownames(AMLexpr) == "?")),
            (which(rownames(AMLexpr) == "SLC35E2")))
geneNames <- Data@rowRanges@elementMetadata@listData[["gene_id"]]
row[duplist] <- geneNames[duplist]
rownames(AMLexpr) <- row

#check if there are any remained dups
rownames(AMLexpr)[duplicated(rownames(AMLexpr))]
(which(rownames(AMLexpr) == "UBE2Q2P2"))
row[19879] <- "UBE2Q2P2-2"
rownames(AMLexpr) <- noquote(row)
rownames(AMLexpr)[duplicated(rownames(AMLexpr))]
#save the table
write.table(AMLexpr, "AMLexpr.txt", quote = FALSE,
            col.names = TRUE, row.names = TRUE,
            sep = "\t")

#clean up
rm( genetable, dups, aws, ungenes,
          duplist, geneNames, row )

###### 7. ESTIMATE------------
# get the address of the file

AMLcancerExpr <- file.path("AMLexpr.txt")
filterCommonGenes(input.f=AMLcancerExpr,
                  output.f="AMLgenes.gct",
                  id="GeneSymbol")
print("OH YES!")

estimateScore("AMLgenes.gct",
              "AML_estimate_score.gct",
              platform="illumina")
print("OMG!")

####### 8. Load gct file and classify patients------------
AMLgct <- file.path("AML_estimate_score.gct")
gcttable <- read.gct(AMLgct)
Scores <- as.data.frame(t(gcttable))





Scores <- arrange(Scores, rownames(Scores))

Scores <- arrange(Scores, desc(ImmuneScore))
Scores$IM <- c(rep(1, 87), rep(2, 86))

Scores$IM <- factor(Scores$IM, 
                       levels = c(1,2),
                       labels = c("high", "low"))

Scores <- arrange(Scores, desc(StromalScore))
Scores$ST <- c(rep(1, 87), rep(2, 86))
Scores$ST <- factor(Scores$ST, 
                    levels = c(1,2),
                    labels = c("high", "low"))
#clean up
rm( gcttable, AMLgct, AMLexpr)







####### 9. Clinical Data--------
clinical <- GDCprepare_clinic(query_clin,"patient")


a <- clinical[,c("leukemia_french_american_british_morphology_code",
                      "gender",
                      "days_to_death",
                      "bcr_patient_barcode",
                      "vital_status",
                      "days_to_last_followup",
                      "acute_myeloid_leukemia_calgb_cytogenetics_risk_category") ]


rownames(a) <- a[,"bcr_patient_barcode"]
a <- subset (a, select = -bcr_patient_barcode)
names(a)[names(a) == "leukemia_french_american_british_morphology_code"] <- "FAB"
names(a)[names(a) == "days_to_death"] <- "OS"
names(a)[names(a) == "acute_myeloid_leukemia_calgb_cytogenetics_risk_category"] <- "calgb"
head(a)

head(Scores)
rownames(Scores) <- sub("X", "TCGA-AB-",
                        rownames(Scores))

#merge tables
Scores <- arrange(Scores,rownames(Scores))
Data_Table <- cbind(Scores, a)

#check if cbind is OK in an idiatic way!

if(identical(Data_Table$FAB, a$FAB) &
   identical(Data_Table$gender, a$gender) &
   identical(Data_Table$OS, a$OS) &
   identical(Scores$IM, Data_Table$IM) &
   identical(Scores$ST, Data_Table$ST) &
   identical(Scores$ESTIMATEScore,
             Data_Table$ESTIMATEScore) &
   identical(Scores$ImmuneScore,
             Data_Table$ImmuneScore) &
   identical(Scores$StromalScore,
             Data_Table$StromalScore)
   ){
  print("OK!")
}

#another way!
aaa <- function(x,y,z){
  tab <- cbind(y,z)
  for (i in 1:dim(tab)[2]){
    if(identical(x[,i], tab[,i])){
      print("OK")
    }
  }
}
aaa <- function(x,y,z){
  tab <- cbind(y,z)
  for (i in 1:dim(tab)[2]){
    if(identical(x[,i], tab[,i])){
      print("OK")
    }
  }
}

aaa(Data_Table, Scores, a)

#clean up
rm(a, aa, aaa, AMLcancerExpr, fabname, MLgct,
   Scores)

#10. Survival Analysis ------------------
?TCGAanalyze_survival
Data_Table$vital_status <- clinical$vital_status
Data_Table$days_to_death <- clinical$days_to_death
Data_Table$days_to_last_follow_up <- clinical$days_to_last_followup
Data_Table <- Data_Table[, -8]
Data_Table$VT <- ifelse(Data_Table$vital_status == "Dead",
                        T, F)
Data_Table$OS <- ifelse(Data_Table$VT,
                        Data_Table$days_to_death,
                        Data_Table$days_to_last_follow_up)
TCGAanalyze_survival(
  Data_Table,
  clusterCol = "ST",
  legend = "Legend",
  labels = NULL,
  risk.table = FALSE,
  xlim = NULL,
  main = "Kaplan-Meier Overall Survival Curves",
  ylab = "Probability of survival",
  xlab = "Time since diagnosis (days)",
  filename = "survival.pdf",
  color = c("Red", "Blue"),
  height = 8,width = 12,
  dpi = 300,
  pvalue = TRUE,
  conf.int = FALSE,
)

######11. plot for FAB ------

#theme
cleanup <- theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(color = "black"),
                 axis.text.x = element_text(angle = 45,
                                            vjust = 0.5,
                                            hjust=1),
                 axis.title.x = element_text(vjust=-0.2),
                 axis.text.x.bottom = element_text(vjust=0.7),
                 legend.position="none")
########################################################
FAB_Immuneplot <- ggplot(Data_Table,
                         aes(FAB, ImmuneScore, fill = FAB,
                             color = FAB)) + cleanup

FAB_Immuneplot + geom_dotplot(binaxis='y',
                              stackdir='center',
                              stackratio=1.5,
                              dotsize=0.5,
                              position = "dodge"
                              ) +
stat_summary(fun.data=mean_sdl,
             fun.args = list(mult=1), 
                 geom="pointrange",
             color="#888888")

### 12. Anova for FAB-------------
fit <- aov(ImmuneScore ~ FAB, data = Data_Table)
a <- anova(fit)
fit2 <- aov(StromalScore ~ FAB, data = Data_Table)
b <- anova(fit)
Data_Table$cyto <- clinical$acute_myeloid_leukemia_calgb_cytogenetics_risk_category
bb <- c("TCGA-AB-2810","TCGA-AB-2895")
for (i in bb){
print(which(rownames(Data_Table) == i))}
cc <- c(6,79)
fit3 <- aov(ImmuneScore ~ cyto, data = Data_Table[-cc,])
c <- anova(fit3)
TukeyHSD(fit3)
###13. differential gene expression analysis -------

##check the data



# get a backup
Data1 <- Data

# add IM and ST to Data
Data@colData@listData$IM <- Data_Table$IM
Data@colData@listData$ST <- Data_Table$ST



#round counts in case they are not integers
Data@assays@data@listData[["raw_count"]] <- round(Data@assays@data@listData[["raw_count"]])
dds <- DESeqDataSet(Data, design = ~ IM + ST)




# omit rows with count-Sum < 1
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)




# normalize the reads for ploting
rld <- vst(dds, blind=FALSE)
head(assay(rld), 3)




# get sample dist
sampleDists <- dist( t( assay(rld) ) )
sampleDists




# heatmap sample distance
sampleDistMatrix <- as.matrix( sampleDists, rownames = T )
rows <- rownames(sampleDistMatrix)
rownames(sampleDistMatrix) <- sub("TCGA-AB-", "", rownames(sampleDistMatrix))
rownames(sampleDistMatrix) <- sub("-.*", "", rownames(sampleDistMatrix))

rownames(sampleDistMatrix) <- paste(rownames(sampleDistMatrix),
                                    rld$IM, rld$ST, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
               clustering_distance_rows=sampleDists,
               clustering_distance_cols=sampleDists,
               col=colors,
               fontsize_row = 3)


poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)






# PCA
plotPCA(rld, intgroup = c("IM", "ST"))
(data <- plotPCA(rld, intgroup = c( "IM", "ST"),
                 returnData=TRUE))

percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2,shape = ST, color=IM)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  stat_ellipse()






## Differential expression analysis
dds <- DESeq(dds)


# check results with FDR 5%
res <- results(dds, alpha = 0.05)



# get IM and ST genes without any filter for volcano plot


resIM <- results(dds, contrast=c("IM", "low", "high"))
resST <- results(dds, contrast=c("ST", "low", "high"))

summary(resIM)
summary(resST)

##plot volcano IM
resIMT <- data.frame(resIM)
resIMT$FC <- ifelse(resIMT$log2FoldChange > log2(1.5) & resIMT$pvalue < 0.05,
                  "red",
                  ifelse(resIMT$log2FoldChange < -log2(1.5) & resIMT$pvalue < 0.05,
                         "blue", "grey"))


ggplot(data=resIMT, aes(x=log2FoldChange,
                       y=-log10(pvalue))) +
  geom_point( size=0.5 ,aes(color=as.factor(FC))) +
  xlim(c(-10, 10)) + ylim(c(0, 33)) +
  xlab("log2 fold change") + 
  ylab("-log10 p-value")  + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  geom_vline(xintercept=c(-log2(1.5), log2(1.5)), col="black", 
             linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="black",
             linetype = "dashed") +
  scale_color_manual(values=c(blue = "royalblue", grey ="grey", red="red2"),
                     labels=c(grey="NS", 
                              blue=paste("FDR Q<", 0.05, sep=""),
                              red=paste("FDR Q<", 0.05, " & LogFC>|",
                                        1.5, "|", sep=""))) + 
  geom_text(aes(8,30,label = "p = 0.05", vjust = -0.5)) + 
  scale_x_continuous(breaks = sort(c(seq(-10, 10, by=5), round(-log2(1.5), 2)
                                     ,round(log2(1.5),2))))



## Volcano ST
resSTT <- data.frame(resST)
resSTT$FC <- ifelse(resSTT$log2FoldChange > log2(1.5) & resSTT$pvalue < 0.05,
                    "red",
                    ifelse(resSTT$log2FoldChange < -log2(1.5) & resSTT$pvalue < 0.05,
                           "blue", "grey"))


ggplot(data=resSTT, aes(x=log2FoldChange,
                        y=-log10(pvalue))) +
  geom_point( size=0.5 ,aes(color=as.factor(FC))) +
  xlim(c(-5.5, 5.5)) + ylim(c(0, 16)) +
  xlab("log2 fold change") + 
  ylab("-log10 p-value")  + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  geom_vline(xintercept=c(-log2(1.5), log2(1.5)), col="black", 
             linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="black",
             linetype = "dashed") +
  scale_color_manual(values=c(blue = "royalblue", grey ="grey", red="red2"),
                     labels=c(grey="NS", 
                              blue=paste("FDR Q<", 0.05, sep=""),
                              red=paste("FDR Q<", 0.05, " & LogFC>|",
                                        1.5, "|", sep=""))) + 
  geom_text(aes(8,14,label = "p = 0.05", vjust = -0.5)) + 
  scale_x_continuous(limits = c(-4.5, 5.5), breaks = sort(c(seq(-5.5, 5.5, by=2.5), round(-log2(1.5), 2)
                                     ,round(log2(1.5),2))))







# MA plot
plotMA(res, ylim=c(-5,5))
plotMA(resIM, ylim=c(-5,5))
plotMA(resST,ylim=c(-5,5))

# It's better to use lfcShrink with apeglm method from apeglm package

##for IM

#normal
resNorm <- lfcShrink(dds, coef=2, type="normal")

plotMA(resNorm, main="normal")

#apeglm
resLFC <- lfcShrink(dds,
                    coef="IM_low_vs_high",
                    type="apeglm")
plotMA(resLFC, main="normal")






#Dispersion Plot
plotDispEsts(dds)

# meanSD plot

meanSdPlot(assay(rld))


#get DEGs for ST & IM

resIMsig <- results(dds, contrast=c("IM", "low", "high"), alpha = 0.05,
                    lfcThreshold = log2(1.5))
resSTsig <- results(dds, contrast=c("ST", "low", "high"), alpha = 0.05,
                    lfcThreshold = log2(1.5))


# convert to data.frame

DEGIM <- data.frame(resIMTSig <- resIMsig[ which(resIMsig$padj < 0.05 & abs(resIMsig$log2FoldChange) > log2(1.5) ), ])

DEGST <- data.frame(resSTTSig <- resSTsig[ which(resSTsig$padj < 0.05 & abs(resSTsig$log2FoldChange) > log2(1.5) ), ])


# Prepare data for a heatmap IM
DEGIM$Genes <- rownames(DEGIM)
DEGIM$Genes <- sub("\\|.*", "", DEGIM$Genes)

DFrld <- data.frame(assay(rld))

DEGT <- DFrld[rownames(DEGIM),]

#change rownames
colnames(DEGT) <- sub("TCGA.AB.", "", colnames(DEGT))
colnames(DEGT) <- sub("\\..*", "", colnames(DEGT))
rownames(DEGT) <- sub("\\|.*", "", rownames(DEGT))



dfIM <- data.frame(colData(dds)["IM"])
rownames(dfIM) <- sub("TCGA-AB-", "", rownames(dfIM))
rownames(dfIM) <- sub("-.*", "", rownames(dfIM))
dfIM$ID <- paste0(rownames(dfIM),"-", dfIM$IM)
rownames(dfIM) <- dfIM$ID
a <- rownames(dfIM)
dfIM <- as.data.frame(dfIM[,1], row.names = a)
colnames(dfIM) <- "IM"

# put high and low together

DEGT_arranged_IM <- t(DEGT)
DEGT_arranged_IM <- cbind(DEGT_arranged_IM, dfIM)
DEGT_arranged_IM <- arrange(DEGT_arranged_IM, desc(DEGT_arranged_IM$IM))
DEGT_arranged_IM <- DEGT_arranged_IM[, 1:561]
DEGT_arranged_IM <- t(DEGT_arranged_IM)
dfIM <- arrange(dfIM, desc(dfIM$IM))

pdf("IMDEGheatmap.pdf")
pheatmap(DEGT_arranged_IM,
         cluster_rows=T,
         show_rownames=T,
         cluster_cols = T,
         clustering_distance_cols = "binary",
         annotation_col = dfIM,
         scale = "row",
         fontsize = 3,
         fontsize_row = 1,
         treeheight_row = 0,
         treeheight_col = 0,
         color = bluered(255),
            border_color= NA)


#manual scale
DEGTscale3 <- t(DEGTscale)
DDD <- cbind(DEGTscale3, dfIM)
DEGTscale3 <- arrange(DDD, DDD$IM)
DEGTscale3 <- t(DEGTscale3[, 1:561])




# Prepare data for a heatmap IM
DEGST$Genes <- rownames(DEGST)
DEGST$Genes <- sub("\\|.*", "", DEGST$Genes)

DFrld <- data.frame(assay(rld))

DEGSTT <- DFrld[rownames(DEGST),]

#change rownames
colnames(DEGSTT) <- sub("TCGA.AB.", "", colnames(DEGSTT))
colnames(DEGSTT) <- sub("\\..*", "", colnames(DEGSTT))
rownames(DEGSTT) <- sub("\\|.*", "", rownames(DEGSTT))



dfST <- data.frame(colData(dds)["ST"])
rownames(dfST) <- sub("TCGA-AB-", "", rownames(dfST))
rownames(dfST) <- sub("-.*", "", rownames(dfST))
dfST$ID <- paste0(rownames(dfST),"-", dfST$ST)
rownames(dfST) <- dfST$ID
b <- rownames(dfST)
dfST <- as.data.frame(dfST[,1], row.names = b)
colnames(dfST) <- "ST"

# put high and low together

DEGSTT_arranged_ST <- t(DEGSTT)
DEGSTT_arranged_ST <- cbind(DEGSTT_arranged_ST, dfST)
DEGSTT_arranged_ST <- arrange(DEGSTT_arranged_ST,
                              desc(DEGSTT_arranged_ST$ST))
DEGSTT_arranged_ST <- DEGSTT_arranged_ST[, 1:219]
DEGSTT_arranged_ST <- t(DEGSTT_arranged_ST)
dfST <- arrange(dfST, desc(dfST$ST))

pdf("STDEGheatmap.pdf")
pheatmap(DEGSTT_arranged_ST,
         cluster_rows=T,
         show_rownames=T,
         cluster_cols = T,
         clustering_distance_cols = "binary",
         annotation_col = dfST,
         scale = "row",
         fontsize = 3,
         fontsize_row = 1,
         treeheight_row = 0,
         treeheight_col = 0,
         color = bluered(255),
         border_color= NA)
################################################################
ddsST <- DESeqDataSet(Data, design = ~ ST)
nrow(ddsST)
ddsST <- ddsST[ rowSums(counts(ddsST)) > 1, ]
nrow(ddsST)
ddsST <- DESeq(ddsST)


STres_sig <- results(ddsST, contrast=c("ST", "low", "high"),
                     alpha = 0.05,
                     lfcThreshold = log2(1.5))
resSTsig <- results(dds, contrast=c("ST", "low", "high"),
                    alpha = 0.05,
                    lfcThreshold = log2(1.5))
summary(STres_sig)

# convert to data.frame

DEGSTH <- data.frame(STres_sig[ which(STres_sig$padj < 0.05 & STres_sig$log2FoldChange > log2(1.5) ), ])
DEGSTL <- data.frame(STres_sig[ which(STres_sig$padj < 0.05 & STres_sig$log2FoldChange < -log2(1.5) ), ])
####################################################
ddsIM <- DESeqDataSet(Data, design = ~ IM)
nrow(ddsIM)
ddsIM <- ddsIM[ rowSums(counts(ddsIM)) > 1, ]
nrow(ddsIM)
ddsIM <- DESeq(ddsIM)


IMres_sig <- results(ddsIM, contrast=c("IM", "low", "high"),
                     alpha = 0.05,
                     lfcThreshold = log2(1.5))

DEGIMH <- data.frame(IMres_sig[ which(IMres_sig$padj < 0.05 & IMres_sig$log2FoldChange > log2(1.5) ), ])
DEGIML <- data.frame(IMres_sig[ which(IMres_sig$padj < 0.05 & IMres_sig$log2FoldChange < -log2(1.5) ), ])
#################################################################
up <- intersect(rownames(DEGSTH), rownames(DEGIMH))
down <- intersect(rownames(DEGSTL), rownames(DEGIML))

Tcount <- assay(Data)
TcountH <- Tcount[up,]
TcountL <- Tcount[down,]


dataAMLhigh <- log2(TcountH)
tokenStop<- 1

HighKM <- NULL

for( i in 1: round(nrow(dataAMLhigh))){
  message( paste( i, "of ", round(nrow(dataAMLhigh))))
  tokenStart <- tokenStop
  tokenStop <-i
  tabSurvKM<-TCGAanalyze_SurvivalKM(clinical,
                                    dataAMLhigh,
                                    Genelist = rownames(dataAMLhigh)[tokenStart:tokenStop],
                                    Survresult = F,
                                    ThreshTop=0.67,
                                    ThreshDown=0.33)
  
  HighKM <- rbind(HighKM,tabSurvKM)
}

HighKM <- HighKM[HighKM$pvalue < 0.05,]
HighKM <- HighKM[order(HighKM$pvalue, decreasing=F),]

dataAMLlow <- log2(TcountL)
tokenStop<- 1

LowKM <- NULL

for( i in 1: round(nrow(dataAMLlow)/100)){
  message( paste( i, "of ", round(nrow(dataAMLlow)/100)))
  tokenStart <- tokenStop
  tokenStop <-i*100
  tabSurvKM<-TCGAanalyze_SurvivalKM(clinical,
                                    dataAMLlow,
                                    Genelist = rownames(dataAMLlow)[tokenStart:tokenStop],
                                    Survresult = F,
                                    ThreshTop=0.67,
                                    ThreshDown=0.33)
  
  LowKM <- rbind(LowKM,tabSurvKM)
}

LowKM <- LowKM[LowKM$pvalue < 0.01,]
LowKM <- LowKM[order(LowKM$pvalue, decreasing=F),]

listsq <- c(rownames(HighKM), rownames(LowKM))
GN <- data.frame(rowData(Data))

GNfinal <- GN[listsq, 2]

GNfinala <- GN[listsq, ]
write.table(GNfinal, "GNfinal.txt", sep = "\t", row.names = F,
            col.names = F, quote = F)
######################################################################

rownames(DEGIM) <- sub("\\|.*", "", rownames(DEGIM))
rownames(DEGST) <- sub("\\|.*", "", rownames(DEGST))

length(intersect(rownames(DEGIM), rownames(DEGST)))
aaa <- Data@rowRanges@elementMetadata@listData[["ensembl_gene_id"]]
aaa$genes <- sub("\\|.*", "", aaa$genes)

aw <- rownames(DEGIM)

DEGIM$EN <- aaa[which(intersect(rownames(DEGIM), aaa$genes)), 2]
s <- intersect(rownames(DEGIM), aaa$genes)
ENIDs <- aaa[which(s %in% aaa$genes),]
  
dd <- intersect(rownames(DEGIM), rownames(DEGST))
commonDEGs <- DEGST[dd,]
write.table(ENIDs$ENSEMBL, "commonDEGs.txt", quote = F, row.names = F,
            sep = ",")
ENIDs <- aaa[which(commonDEGs$Genes %in% aaa$genes),]

#survival each gene
rldexpr <- assay(vst(ddsIM, blind = FALSE))
rldexpr <- rldexpr[common,]
rlddata <- data.frame(Data_Table[,11:12], t(rldexpr))

res5 <- RegParallel(
  data = rlddata,
  formula = 'Surv(OS, VT) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(rlddata)[3:ncol(rlddata)],
  blocksize = 481,
  p.adjust = "BH")
res5 <- res5[!is.na(res5$P),]
res5 <- res5[order(res5$LogRank, decreasing = FALSE),]
final <- subset(res5, LogRank < 0.05)
final$Variable <- gsub("\\..*", "",final$Variable)

genes <- final$Variable
colnames(rlddata2) <- gsub("\\..*", "",colnames(rlddata2))

rlddata2 <- rlddata
rldgenes <- rlddata2[,c("OS","VT",genes)]
library(survminer)
apply(rldgenes, 2, median)
apply(rldgenes, 2, mean)

rldgenes$GNGT2 <- ifelse(rldgenes$GNGT2.2793 >= median(rldgenes$GNGT2.2793),
                         "High", "Low")
rldgenes <- rldgenes[,-83]
rldgenes[,3:ncol(rldgenes)] <- apply(rldgenes[,3:ncol(rldgenes)], 2,
                                     function(x){ifelse(x >= mean(x),
                                                        "High", "Low")})

rldgenes2 <- rlddata2[,c("OS","VT",genes)]
rldgenes2[,3:ncol(rldgenes2)] <- apply(rldgenes2[,3:ncol(rldgenes2)], 2,
                                     function(x){ifelse(x >= median(x),
                                                        "High", "Low")})

r
KM1 <- final$Variable[1:8]
KM1 <-gsub("\\..*", "",KM1)
a <- which("IL1R2.7850" == colnames(rldgenes))
for (i in KM1){
  pdf(paste0(i, ".pdf"), onefile = FALSE)
  print(ggsurvplot(survfit(as.formula(paste0("Surv(OS, VT) ~", i)),
                     data = rldgenes),
             data = rldgenes,
             risk.table = FALSE,
             pval = TRUE,
             break.time.by = 500,
             ggtheme = theme_minimal(),
             risk.table.y.text.col = FALSE,
             risk.table.y.text = FALSE))
  while (!is.null(dev.list()))  dev.off()}


ggsurvplot(survfit(Surv(OS, VT) ~ IL10,
                   data = rldgenes2),
           data = rldgenes2,
           risk.table = TRUE,
           pval = TRUE,
           break.time.by = 500,
           ggtheme = theme_classic(),
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE)

setwd("/media/rk1994/01CD77B76BF0B4F0/R/Results/KEGG/")
dev.off()
