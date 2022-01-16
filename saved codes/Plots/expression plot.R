library(ggpubr)
#SW48------------------------
ggbarplot(SW4, x = "Cell", y = "Count",
          fill = "Type",           # change fill color by mpg_level
          color = "Type",            # Set bar border colors to white
          palette = get_palette("jco", 5)[c(1,4)],            # jco journal color palett. see ?ggpar
                    # Sort the value in ascending order
               # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Relative expression (Log10)",
          xlab = FALSE,
          legend.title = "Groups",
          position = position_dodge(0.8),
          
          error.plot = "errorbar",
          add = "mean_sd"
          
) + labs_pubr() + 
  stat_pvalue_manual(
    p, x = "Cell", y.position = "y.position",
    label = "p = {p}",
    bracket.size = 0.3,
   tip.length = 0,
    bracket.nudge.y = -2
  ) + 
  geom_errorbar(aes(ymin = value, ymax = value + error, color = variable), data = SW3, position = position_dodge(0.8), width = 0.2)
  

rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }

  #+ scale_y_continuous(expand = expansion(mult = c(0, 0.1))) puts y and x at zero

SW3 <- data.frame(Cell = c(rep("SW48",6), rep("HT29", 6)),
                  Type = c(rep(c("Control", "Transfect" ),each = 3, 2)),
                  Count = c(round(rnorm2(3,325,25)),
                            round(rnorm2(3,980,45)),
                            round(rnorm2(3,256,20)),
                            round(rnorm2(3,893,32))))
t.test(SW3[1:3,3],SW3[4:6,3])                      
t.test(SW3[7:9,3],SW3[10:12,3])                      
SW4 <- SW3

SW4[,3] <- log10(SW3[,3]+1)

SW3 <- data.frame(Cell = c(rep("SW48",6), rep("HT29", 6)),
                  Type = c(rep(c("Control", "Transfect" ),each = 3, 2)),
                  Count = c(round(rnorm2(3,1,0)),
                            round(rnorm2(3,70.03479688,10)),
                            round(rnorm2(3,0,0)),
                            round(rnorm2(3,234.7530351,32))))


ggboxplot(With_patient_tumor_type, x = "type", y = "ddct-",
               fill = "type", palette = get_palette("jco", 5)[c(1,4)],
               add = "jitter", add.param = list(color = "Black"),shape = "type", width = 0.5)+ xlab("Samples") + ylab("Log2 Relative Expression \n of RNF43-SUPT4H1") + theme_pubr()+ labs_pubr()

ggplot(With_patient_tumor_type, aes(x = type, y = `ddct-`,color = factor(`type`),fill = factor(`type`))) + geom_dotplot(binaxis='y', stackdir='center',
                                                                                                  stackratio=1.5, dotsize=1.2)+labs_pubr() +stat_summary(fun.data=mean_se, mult=1,geom="pointrange", color="Black")+ylab("Log2 Relative Expression \n of RNF43-SUPT4H1") + theme_pubr()+ labs_pubr()+
  scale_color_manual(values = c("#0073C2FF", "#CD534CFF"))+scale_fill_manual(values = c("#0073C2FF", "#CD534CFF"))

library(ggpubr)

ggviolin(With_patient_tumor_type,
         x = "type",
         y = "ddct-",
         fill = "type",
         palette = get_palette("jco", 5)[c(1,4)],
         add = "boxplot",
         add.params = list(fill = "white"))+geom_jitter()+theme_pubr()+ labs_pubr()+ scale_y_continuous(breaks = seq(-6, 6, 2))+ylab("Log2 Relative Expression \n of ADAP1-NOC4L")



ggstripchart(With_patient_tumor_type, "type", "ddct-",  shape = "type",size = 3,
             color = "type", palette = c("#00AFBB", "#FC4E07"),
             add = "mean_se", add.params = list(color = "black"))+theme_pubr()+ labs_pubr()+ scale_y_continuous(breaks = seq(-6, 6, 2))+ylab("Log2 Relative Expression \n of RNF43-SUPT4H1")
             )
palette = c("#00AFBB", "#E7B800", "#FC4E07")