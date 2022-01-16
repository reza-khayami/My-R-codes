setwd("C:/Users/RK1994/Desktop/cyto")
#load data
cc <- read.table(file = 'GO_CC.txt', sep = '\t', header = T,comment.char = "",quote="")
BP <- read.table("GO_BP.txt", sep = '\t', header = T,comment.char = "",quote="")
BP$cat <- rep("Biological Processes",74)


cc$cat <- rep("Cellular Components",18)

MF <- read.table("GO_MF.txt", sep = '\t', header = T,comment.char = "", quote="")
MF$cat <- rep("Molecular Functions",19)

GO <- rbind(cc,BP,MF)


#get top30
library("dplyr")
GO <- arrange(GO, GO$FDR)
GO30 <- GO[1:30,]
colnames(GO30)[c(2,3, 6,9)] <- c("GO Terms", "Gene Count", "False Discovery Rate", "Category")

GO30$Term <- sub(".*~","",GO30$Term)
GO30 <- arrange(GO30, GO30$X.)
GO30 <- arrange(GO30, GO30$Category)
GO30$Term <- factor(GO30$Term,
                                levels = unique(GO30$Term))



#plot
pdf("GO.pdf", width = 10, height = 10)
ggplot(GO30, aes(X., Term, color=-log10(FDR), shape= Category, size = Count))  +
  geom_jitter() +
  theme_classic() +
  labs(x="Percentage (Involved Genes/Total Genes)",
      y="",
      size="Gene Count",
       col="-log10(False Discovery Rate)",
      shape="Category") +scale_color_distiller(palette = "PRGn", direction = -1) + theme(text = element_text(size=10)) + ggtitle("GOTERMS")
dev.off()







