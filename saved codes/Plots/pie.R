pie <- data.frame(Labs = c("Missense",
                           "Synonymous",
                           "Non-coding",
                           "Frameshift",
                           "Splice junction loss",
                           "nonsense",
                           "Others"),
                  Count = c(525,170,98,79,60,53,15),
                  Percent = c(52, 17, 9.8, 7.9, 6, 5.3, 1.5))
library(ggplot2)
library(dplyr)

#Calculate positions
count.data <- pie %>%
  arrange(desc(Labs)) %>%
  mutate(lab.ypos = cumsum(Percent) - 0.5*Percent)

#Prepare labels
count.data$Percent2 <- paste0(count.data$Count," (", count.data$Percent, "%)")
mycols <- c("#f29e08", "#4cda41", "#fe4614","#343434", "#224483", "#c037a4", "#d21a0c" )

#----------------------------------------------------
ggplot(count.data, aes(x = "", y = Percent, fill = Labs)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = mycols)+
  geom_text(aes(y = lab.ypos, label = Percent2), color = "white")+
   theme_void()+ labs(fill = "Variant Type")+
  scale_fill_brewer(palette="Set1")

#---------------------------------------------
library(ggrepel)
ggplot(count.data, aes(x = "", y = Percent, fill = Labs)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0)+
  geom_label_repel(aes(y = lab.ypos,label = Percent2), size=5, show.legend = F, nudge_x = 1)+
  scale_fill_brewer(palette="Set3")+
  theme_void() + labs(fill = "Variant Type")

library(ggpubr)
library(ggsci)
ggpie(count.data, "Count", label = "Percent2",
      lab.pos = "in", lab.font = "white",
      fill = "Labs", color = "white",
      pallete = "jco")
