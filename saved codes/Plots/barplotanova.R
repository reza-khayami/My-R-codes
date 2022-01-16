library(readxl)
library(ggpubr)
library(ggsci)
library(rstatix)
library(tidyverse)

#Load data----------------------------------------------
file <- "D:/university files/MSc/Thesis/SLC30A10-3/results/mtt table.xlsx"
mtt_table <- read_excel(file, 
                        sheet = "Sheet2",
                        col_types = c("numeric",
                                      "numeric",
                                      "numeric",
                                      "numeric",
                                      "text"))
#melt
mtt_table2 <- as.tbl(reshape2::melt(mtt_table))

#We have to change colnames after mel because later the function  add_xy_position brain farts if a column is named variable
mtt_table2$Time <- as.factor(mtt_table2$Time)
mtt_table2$Time <-  factor(mtt_table2$Time,levels(mtt_table2$Time)[c(2:5,1)])

colnames(mtt_table2)[2] <- "Groups"


#First check data normality------------------------------
ggdensity(mtt_table2[mtt_table2$Time == "120h",]$value, 
          main = "Density plot of value",
          xlab = "")
ggqqplot(mtt_table2$value)
shapiro.test(mtt_table2$value)

#If the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality. Otherwise its not normal.


#If data is normal:---------------------------------------
#Anova-----------
stat.testb = list()

for (i in levels(as.factor(mtt_table2$Time))){
stat.testb[[i]] <-  aov(value ~ Groups, data = mtt_table2[mtt_table2$Time == i,]) %>% tukey_hsd()

}
stat.test1 = do.call(rbind, stat.testb)
stat.test2 <- stat.test1 %>% filter(group1 == "control")



stat.test2$Time <- rep(names(stat.testb), each = 3)
stat.test2$Time <- as.factor(stat.test2$Time)
stat.test3 <- stat.test2[,c("Time","term","group1", "group2","null.value","estimate","conf.low" ,"conf.high","p.adj" , "p.adj.signif"  )]
stat.test3$term <- ifelse(stat.test3$term == "Groups", "value","")
colnames(stat.test3)[2] <- ".y."


#I don't know why the code below doesn't work so I do it another way
stat.test4 <- stat.test3 %>%
  add_xy_position(x = "Time", fun = "mean_se", dodge = 0.8)


#here is the workaround
stat.test <- mtt_table2 %>%
  group_by(Time) %>%
  t_test(value ~ Groups, ref.group = "control") 

stat.test <- stat.test %>%
  add_xy_position(x = "Time", fun = "mean_se", dodge = 0.8)

stat.test4 <- cbind(stat.test3, stat.test[,-c(1:11)])


bp <- ggbarplot(mtt_table2,
       x = "Time",
       y = "value",
       fill = "Groups",
       color = "black",
       add = "mean_se",
       legend = "right",
       ylab = "Cell viability (%)",
       position = position_dodge(0.8)
       )  +
  labs_pubr()+
  font("xy", size = 9)+
  font("x.text", size = 9)+
  font("y.text", size = 9)+ 
  font("legend.title", size = 7) +
  font("legend.text", size = 7) +
  scale_y_continuous(breaks = seq(0,100, 20))+
  scale_fill_manual(values = c("#E5E5E5","#F1C0E8", "#A3C4F3", "#CFBAF0"))+
  scale_color_manual(values = c("#E5E5E5","#F1C0E8", "#A3C4F3", "#CFBAF0")) 
  

bp + 
  stat_pvalue_manual(
    stat.test4[stat.test4$p.adj.signif != "ns",], label = "p.adj.signif", tip.length = 0.01,
    bracket.nudge.y = -2,
    y.position = c(120,125)
  )
#t-test-----------------------------------

stat.test <- mtt_table2 %>%
  group_by(Time) %>%
  t_test(value ~ Groups, ref.group = "Control", p.adjust.method = "BH") 

stat.test <- stat.test %>%
  add_xy_position(x = "Time", fun = "mean_se", dodge = 0.8)


bp <- ggbarplot(mtt_table2,
                x = "Time",
                y = "value",
                fill = "Groups",
                color = "Groups",
                add = "mean_se",
                legend = "right",
                ylab = "Cell viability (%)",
                position = position_dodge(0.8)
)  +
  labs_pubr()+
  font("xy", size = 9)+
  font("x.text", size = 9)+
  font("y.text", size = 9)+ 
  font("legend.title", size = 7) +
  font("legend.text", size = 7) +
  scale_y_continuous(breaks = seq(0,100, 20))+
  scale_fill_manual(values = c("#E5E5E5","#F1C0E8", "#A3C4F3", "#CFBAF0"))+
  scale_color_manual(values = c("#E5E5E5","#F1C0E8", "#A3C4F3", "#CFBAF0"))+
  guides(shape = guide_legend(override.aes = list(size = 0.5))) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
   theme(legend.key.size = unit(0.4, "cm"))
         


bpp <- bp + 
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif", tip.length = 0.01,
    bracket.nudge.y = -2
  )
  
  #If data is not normal:---------------------------------

stat.test5 <- mtt_table2 %>%
  group_by(Time) %>%
  wilcox_test(value ~ Groups, ref = "control") 

stat.test <- stat.test %>%
  add_xy_position(x = "Time", fun = "mean_se", dodge = 0.8)

bp + stat_compare_means(aes(group = Groups))



#save the file

ggsave("D:/university files/MSc/Thesis/SLC30A10-3/results/MTT.tiff", bpp, width = 15,height = 7, units = "cm",dpi = 300)
