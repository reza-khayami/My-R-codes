# Download the data from github (click the "raw" button, save as a text file called "results.txt").
# https://gist.github.com/stephenturner/806e31fce55a8b7175af
res <- read.table("R.working/results.txt", header=TRUE)
head(res)
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))



###ggplot
fold_changes <- c(rnorm(20000, 0, 2))
pvalues <- runif(n=20000, min=1e-50, max=.1)
dif <- data.frame(fc =fold_changes,pv =pvalues)
dif$thershold <- ifelse(dif$fc > 1 & dif$pv < 0.01, "red", 
                        ifelse(dif$fc < -1 & dif$pv < 0.01, -1, "blue"))
ggplot(data=dif, aes(x=fc, y=-log10(pv))) +
  geom_point( size=1 ,aes(color=as.factor(thershold))) +
  xlim(c(-10, 10)) + ylim(c(0, 6)) +
  xlab("log2 fold change") + 
  ylab("-log10 p-value")  + 
  annotate("label", x =c(-8,5), y = 4.75, label = c("400","120"), col=c("red","steelblue"))+
  annotate("text", x =c(-8,5), y = 5, label = c("50 FC>4","8FC <-4"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")+  # for manual coloring


scale_color_manual(values=c(NS="grey30", FC="forestgreen", FDR="royalblue", FC_FDR="red2"), labels=c(NS="NS", FC=paste("LogFC>|", FCCutoff, "|", sep=""), FDR=paste("FDR Q<", AdjustedCutoff, sep=""), FC_FDR=paste("FDR Q<", AdjustedCutoff, " & LogFC>|", FCCutoff, "|", sep=""))) 
  
geom_text(aes(-8,0,label = "p = 0.05", vjust = -0.5)) + scale_x_continuous(breaks = sort(c(seq(-10, 10, by=5), -1,1)))