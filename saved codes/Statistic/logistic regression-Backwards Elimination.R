setwd("H:/R/statsimperial/")
g <- read.csv("NVJ7Cw98Eem6Pgo4-YwqLg_35bb00c00f7c11e9903947c521ebe81a_final-diabetes-data-for-R-_csv_-_2_.csv")
g <- g[-1,]

##### Make the variables and run the models #####

dm <- as.factor(g[,"dm"]) 
insurance <- as.factor(g[,"insurance"])# let's say 0=none, 1=gov, 2=private 
fh <- as.factor(g[,"fh"]) # 1=FH, 0=no FH 
smoking <- as.factor(g[,"smoking"]) # 1,2,3 
chol <- g[,'chol'] 
hdl <- g[,'hdl'] 
ratio <- g[,'ratio'] 
location <- as.factor(g[,'location']) 
age <- g[,'age'] 
gender <- as.factor(g[,'gender']) 
frame <- as.factor(g[,'frame']) 
systolic <- g[,'bp.1s'] 
diastolic <- g[,'bp.1d']
height <- g[,"height"] 
weight <- g[,"weight"] 
height.si <- height*0.0254 
weight.si <- weight*0.453592 
bmi <- weight.si/height.si^2 

null_model <- glm(dm ~ 1, family = binomial(link = logit)) 
summary(null_model) 

model <- glm(dm ~ age + bmi + chol + hdl + systolic + diastolic, family = binomial(link = logit)) 

summary(model) 

anova(model, test = "Chisq") 

#It's clear that neither of the BP variables is significantly associated with the odds of being diagnosed with diabetes in this data set, but the other four variables were. If you drop the BP variables, you get this: 

model2 <- glm(formula = dm ~ age + bmi + chol + hdl, family = binomial(link = logit)) 
summary(model2)

#Have any of the coefficients for the four remaining variables changed? Not much, which is good. But why is blood pressure not significant here despite what the literature says? One way to find out is to see if it correlates with other variables. Here's the code to do that and the output.

cor.test(systolic, hdl)
cor.test(systolic, bmi)
cor.test(systolic, chol)
cor.test(systolic, age)
# So systolic BP correlates weakly (but statistically significantly) with cholesterol and moderately (and also statistically significantly) with age. Both of these results are entirely expected from what we know about physiology. 

model3 <- glm(formula = dm ~ bmi + chol + hdl, family = binomial(link = logit)) 
summary(model3)


model4 <- glm(dm ~ age + bmi + chol + hdl + systolic + diastolic + gender + location + frame + insurance + smoking, family = binomial(link = logit)) 

summary(model4) 

anova(model4, test = "Chisq") 

# The majority of what you included was not significantly associated with the outcome. There's potential for "clearing out the garbage" by removing the non-significant ones, but you would need to check that this doesn't alter the odds ratios for the variables that remain in the model. Also, with only 11 variables, the table of results isn't too big - and it's important to know which variables were not associated with the outcome in this data set. So-called "negative findings" are just as important to know about and report as "positive findings".
