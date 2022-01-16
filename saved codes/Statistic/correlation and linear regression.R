COPD <- read.csv("COPD.csv")
#corelation-------
hist(COPD$MWT1Best)
hist(COPD$MWT1Best, main="Histogram of MWT1Best", xlab="MWT1Best", breaks=12)

#there is quite a high value which is above 650.
#subset(dataframe, variable > x), so here: 
subset(COPD, MWT1Best > 650)

# If you now want to look at more than one value, let's say samples which have values of MWT1Best over 650 and under 150, you can use the following code:  
subset(COPD, MWT1Best > 600 | MWT1Best < 150) #where '|' stands for 'and'.

# To view descriptive statistics for MWT1Best, type the command:  
list("Summary" = summary(COPD$MWT1Best), "Mean" = mean(COPD$MWT1Best, na.rm=TRUE), "Standard Deviation" = sd(COPD$MWT1Best, na.rm=TRUE), "Range" = range(COPD$MWT1Best, na.rm=TRUE), "Inter-Quartile Range" = IQR(COPD$MWT1Best, na.rm=TRUE)) 

# And finally, let's calculate the correlation coefficient and have a look at the scatterplot of the two variables. The basic command for this is plot(x,y) and you can rename the axes in the same way you did in the histogram.
# So in this case: 
plot(COPD$FEV1, COPD$MWT1Best, xlab = "FEV1", ylab = "MWT1Best") 

# The basic command for a correlation test is cor.test(x,y) where you can specify which method you want to use using the command  method = "pearson" or "spearman". You need to remove missing values, otherwise you will have an error message. 
#To do this, use the command use#= "complete.obs" 

#linear regression--------
# The basic format of a linear regression is: 
#   
# Y = ?? + ??1*X + ?? 
# 
# Where: 
#   
# Y = outcome (i.e. dependent) variable 
# X = predictor (i.e. independent) variable 
# ?? and ?? are parameters of the regression, with ?? = intercept (average Y when X=0), and ?? = slope of the line (change in Y for a 1 unit increase in X). Note: ?? and ?? are unit specific, so you'll get different answers if you use distance in metres and in feet. 
# ?? is the random variation in Y, i.e. the residuals 
# To run a linear regression in R, the function you need to use is lm(). The basic format of this function is: 

lm(outcome ~ predictor, data =dataframe)
MWT1Best_FEV1 <- lm(MWT1Best~FEV1, data = COPD)
summary(MWT1Best_FEV1) 

# The (Intercept) indicates the regression constant ??.  
# FEV1 indicates the linear effect of lung function, i.e. ??. 

# To view 95% confidence intervals, use the command 
confint(MWT1Best_FEV1)

# To check model assumptions, you can graphically examine the linear regression model using the function plot(). This function allows to check for linearity, homoscedasticity, independence, and normality of your assumptions. Four plots are generated: 
#   
#   The first is a constant variance plot, which checks for the homogeneity of the variance and the linear relation. If you see no pattern in this graph, then your assumptions are met. 
# The second plot is a Q-Q plot, which checks that the residuals follow a normal distribution. The points should fall on a line if the normality assumption is met. 
# The third plot allows to detect heterogeneity of the variance. 
# The fourth plot allows for the detection of points that have a large impact on the regression coefficients. 
par(mfrow=c(2,2)) 
plot(MWT1Best_FEV1)




#Multiple Regression ------
# The basic format of a multiple linear regression is: 
#   
#   Y = ?? + ??1*X1 + ??2*X2 + ?? 
# 
# Where: 
#   
#   Y = outcome (i.e. dependent) variable. 
# 
# X1 = first predictor (i.e. independent) variable. 
# 
# X2 = second predictor (i.e. independent) variable. 
# 
# ?? = intercept (average Y when X1=X2=0). Note: ?? is unit specific. 
# 
# ??1 = slope of the line (change in Y for a 1 unit increase in X1 when X2 is held constant). Note: ??1 is unit specific. 
# 
# ??2 = slope of the line (change in Y for a 1 unit increase in X2 when X1 is held constant). Note: ??2 is unit specific. 
# 
# ?? is the random variation in Y, i.e. the residuals. 

# To run a multiple linear regression in R, the basic format of the function is very similar to that of a simple linear regression - the only difference is that you are adding two predictor variables instead of one: 
#   
# Model name <- lm(outcome ~ predictor1 + predictor2, data =dataframe) 
# 
# So, for example, if you want to run a multiple linear regression to check whether lung function and age are predictors of walking distance in COPD patients, the R command is: 
  
MWT1Best_FEV1_AGE <- lm(MWT1Best~FEV1+AGE, data = COPD) 
summary(MWT1Best_FEV1_AGE)
confint(MWT1Best_FEV1_AGE) 

###Simple Code------

dataframe <- COPD #type your dataframe
  
#single variable
reg <- formula(MWT1Best ~ AGE)


lmfit <- lm(reg, data =dataframe)
plot(reg, data =dataframe)
abline(lmfit)
summary(lmfit)
confint(lmfit)
par(mfrow=c(2,2)) 
plot(lmfit)
par(mfrow=c(1,1)) 

predict(lmfit)
residuals(lmfit)

#multiple regression
Mreg <- formula(MWT1Best~FEV1+AGE)

Mlmfit <- lm(Mreg, data = COPD) 
summary(Mlmfit)
confint(Mlmfit) 