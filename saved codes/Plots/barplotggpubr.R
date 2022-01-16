library(readxl)
temp1 <- read_excel("temp1.xlsx")
temp1 <- read_excel("tem2.xlsx")

View(temp1)
temp1$variable <- factor(temp1$variable)
temp1$classification <- factor(temp1$classification)

temp2 <- reshape2::melt(temp1)
colnames(temp2) <- c("variable", "Classification", "Group", "Frequency")

myplots <- list()
for (i in 1:length(levels(temp2$variable))){
myplots[[i]] <- ggpar(ggbarplot(
  temp2[temp2$variable== levels(temp2$variable)[i],],
  x = "Group",
  y = "Frequency", 
  fill = "Classification",
  position = position_dodge(0.9),
  palette = get_palette(c("npg"), 10),
  title = LETTERS[i],
  ylab = "Frequency",
  color = "black",
  legend = "right")  +
  labs_pubr()+
  font("xy", size = 7)+
  font("x.text", size = 7)+
  font("y.text", size = 7)+ 
  font("legend.title", size = 7) +
  font("legend.text", size = 7) +
  scale_y_continuous(breaks = seq(0,70, 10)),legend.title = levels(temp2$variable)[i] )
  }



w <- ggarrange(plotlist =myplots, ncol = 3, nrow = 3) 
ggsave("table2.png", plot = w, device = "png", width = 3840, height = 2160, units = "px", dpi = 300)



