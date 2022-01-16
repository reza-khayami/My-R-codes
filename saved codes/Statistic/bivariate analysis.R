###Bivariate analysis
rm(list = ls())
require(car)
cor(Prestige$income, Prestige$education)
cor.test(Prestige$income, Prestige$education)
#Anova
# NULL : all population means are equal
#Alternate : at least one population mean is different
# F-ratio = mean between groups S of squares/ mean within ...
boxplot( income ~ type, data = Prestige)
anovamod <- aov(income ~ type, data = Prestige)
summary(anovamod)
# Mean Sq of type = variance in income between the groups of type variable
#Mean Sq of residuals = variances within the groups of type variables
#F statistics = ratio between/within 
TukeyHSD(anovamod, conf.level= 0.95 )
#Chisquare
Prestige$income_cat <- dplyr::ntile(Prestige$income, 4)
table(Prestige$income_cat, Prestige$type)
chisq.test(y = Prestige$income_cat, x = Prestige$type)
# challange
R <- aov(income ~ education + women + prestige + census, type , data = Prestige)
summary(R)
