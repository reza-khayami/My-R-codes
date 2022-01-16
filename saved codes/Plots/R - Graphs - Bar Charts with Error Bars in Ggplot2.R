rm(list= ls())
####clean up
cleanup <- theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(color = "black"),
                #legend.position = "none"
                )

chick <- read.csv("C://Users//RK1994//Downloads//c4 ChickFlick.csv", header = TRUE)
head(chick)
colnames(chick) <- c("gender", "film", "arousal")

chick$gender <- factor(chick$gender, 
                       levels = c(1,2),
                       labels = c("Men", "Women"))
chick$film <- factor(chick$film, 
                     levels = c(1,2),
                     labels = c("Bridget Jones", "Memento"))
# Bar chart one indipendent variable
chickbar <- ggplot(chick, aes(film, arousal))
chickbar + 
  stat_summary(fun.y = mean, geom = "bar", fill = "white", color = "black" ) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width = 0.90), width = 0.02)+ cleanup + xlab("Movie Watched by Participants") + ylab("Arousal Levels")
# Bar chart : Two independent variables

chickbar2 <- ggplot(chick, aes(film, arousal, fill = gender))
chickbar2 + stat_summary(fun.y= mean, geom = "bar", position= "dodge") + stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width = 0.90), width = 0.2) + xlab ("Film Watched")+ ylab("Arousal Level") + cleanup + scale_fill_manual(name = "Gender of Participants",values = c("black", "grey"))

                               