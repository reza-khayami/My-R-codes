# Course Title: Mastering R Programming
# Author: Selva Prabhakaran
# URL: https://www.packtpub.com/big-data-and-business-intelligence/mastering-r-programming-video
# Contact: selva86@gmail.com
# Website: www.r-statistics.co


# Dataprep URL: https://goo.gl/ZTd9wf
# Alternate Dataprep URL: http://ow.ly/jXnR305qq2v

#OLS assumptions
# the regression model is linear in parameters
# the mean of residuals is zero, conditional on independant variables
# the predictors are linearly independant of each other
# the residuals are uncorrelated
# homoscdasticity (the variance of error is constant across the observations)
# normality of residuals


# Prep
data(Prestige)
set.seed(100)
train_rows <- sample(1:nrow(Prestige), size=0.7*nrow(Prestige))
training <- Prestige[train_rows, ]
test <- Prestige[-train_rows, ]

# Model
lmmod <- lm(prestige ~ income + education, data=training)
summary(lmmod)
par(mfrow=c(2,2))
plot(lmmod)

#residual vs fitted ->> if the error variance increses as the fitted values increase = hetroscdasticity
#this condition should be avoided to ensure the model does not bahave unpredictably  

#normal Q-Q plot is to make sure the residuals follow a normal distribution
#leverage is a meseare of how much each data point influences the regression

lmtest::bptest(lmmod)
#if p > 0.05 no hetroscdasticity becuase the error is constant

cooks.distance(lmmod) #shows how much the fiited values would change if an observation was removed
#points that are four times the mean are considered as influencial



car::influenceIndexPlot(lmmod, id.n = 5)

#which features are contributing to the curvature 
car::residualPlots(lmmod)

#there is a significant curvature in income so we use log(income) which is a commonly used transformation to dampen this effect
lmmod2 <- lm(prestige ~ log(income) + education, data=training)
car::residualPlots(lmmod2)
summary(lmmod2)
