setwd("E:/university files/MSc/Articles/Dr Kerachian/RS96/results/RPI/RPIseq/new/")
files <- list.files('.')
cn <- lapply(files, read.delim, header = T,  #
             stringsAsFactors = F)           #
lnc <- lapply(files, fread)   
names(lnc) <- files

cn <- do.call(cbind, cn)


for (i in 1:6){
  lnc[[i]][[1]] <- gsub(".*\\|", "", lnc[[i]][[1]])
}
for (i in 1:6){
  lnc[[i]][[1]] <- sub("_HUMAN.*", "", lnc[[i]][[1]])
}

#ggplot

lncmelt <- melt(lnc)
colnames(lncmelt)[2:4] <- c("Method", "Score", "Transcript")
lncmelt$Transcript <- sub( ".txt", "", lncmelt$Transcript)


mine.heatmap <- ggplot(data = lncmelt,
                       mapping = aes(x = Transcript,                                                                   y = Header,                                                                   fill = Score)) +
  geom_tile() +
  facet_grid(~ Method)+
  scale_fill_gradient(low = "#FDE725", high = "#440154")  +
  theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=10)) 

lncmelt$Tra <- paste0(lncmelt$Method, lncmelt$Transcript)

pdf("11.pdf")
dev.off()
