library(ggpubr)
library(ggsci)
a <- data.frame(Status = rep(c("Lnc", "Control"), each = 2),
                Necrosis = c(71.4,58.9, 15.4, 9.44),
                `Late apoptosis` = c(20.5, 16.3, 55.9, 48.5),
                `Early apoptosis` = c(2.04, 3.81, 5.43,7.83))
a <- reshape2::melt(a)
a$variable[which(a$variable == "Late.apoptosis")] <- "Late apoptosis"
compare_means(value ~ Status, data = a, 
          group.by = "variable", method = "t.test")


ggbarplot(a,
          x = "variable",
          y = "value",
          fill = "Status",
          add = "mean_se",
          palette = pal_jco()(6)[c(5,6)],
          legend = "none",
          color = "black",
          ylab = "Percentage",
          xlab = "",
          position = position_dodge(0.8))  + labs_pubr()+  font("x.text", size = 9)+
  stat_compare_means(aes(group = Status), label = "p.signif", label.y = 80, data = a, method = "t.test") +
  font("y.text", size = 9)
ggsave("Aza_apoptos.pdf", )