library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(rstatix)

bxp <- ggboxplot(total, x = "Status", y = "Folic", fill = "Status", outlier.shape = NA) + ylab("Folic Acid (ng/ml)") 

stat.test <- total %>%
  pairwise_wilcox_test(Folic ~ Status) %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "Status")

bxpf <- bxp + stat_pvalue_manual(
  stat.test, label = "p = {p}",
  vjust = -1, bracket.nudge.y = 1
) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
a <- ggpar(bxpf, legend = "none")


bxp2 <- ggboxplot(total[-32,], x = "Status", y = "V12", fill = "Status", outlier.shape = NA) + ylab("Vit B12 (pg/ml)") 

stat.test <- total[-32] %>%
  pairwise_wilcox_test(V12 ~ Status) %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "Status")

bxp2f <- bxp2 + stat_pvalue_manual(
  stat.test, label = "p = {p}",
  vjust = -1, bracket.nudge.y = 1
) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
b <- ggpar(bxp2f, legend = "none")

grid.arrange(a, b, nrow = 1)
