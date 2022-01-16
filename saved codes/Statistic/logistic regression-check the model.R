############### 

# grouping cholesterol into categories 
chol_categorised <- ifelse(chol < 200, "healthy",  
                           ifelse(chol < 240, "borderline high", 
                                  ifelse(chol >= 240, "high", NA))) 


# plotting cholesterol categories vs log odds of diabetes 
############### 

# 1. make sure that it is treated as a factor/categorical variable and ordering the levels within the factor for the table 
chol_categorised <- factor(chol_categorised, levels = c("healthy", "borderline high", "high")) 

# 2. create a cross tabulation of cholesterol and diabetes status  
dm_by_chol_categorised <- table(chol_categorised, dm) # not including NA values because there aren't that many 

# 3. output the frequencies of diabetes status by cholesterol 
dm_by_chol_categorised_prop <- prop.table(dm_by_chol_categorised, margin = 1) 

# 4. calculate the odds of having diabetes 
odds_chol_categorised <- dm_by_chol_categorised_prop[, "yes"]/dm_by_chol_categorised_prop[, "no"] 

# 5. calculate the log odds 
logodds_chol_categorised <- log(odds_chol_categorised) 

# 6. plot the cholesterol found in the sample against the log odds of having diabetes 
dotchart(logodds_chol_categorised) 

############### 


# run a multiple logistic regression, generate odds of diabetes (95% CI) for each predictor variable using age, gender and BMI as predictor variables 
############### 

# 1. generate full multiple logistic regression 
full_model <- glm(dm ~ age + gender + bmi, family = binomial(link = logit)) 
summary(full_model) 
##  
## Call: 
## glm(formula = dm ~ age + gender + bmi, family = binomial(link = logit)) 
##  
## Deviance Residuals:  
## 	Min   	1Q   Median   	3Q  	Max   
## -1.6843  -0.5763  -0.3885  -0.2575   2.6991   
##  
## Coefficients: 
##           	Estimate Std. Error z value Pr(>|z|)     
## (Intercept)  -6.647817   0.961731  -6.912 4.77e-12 *** 
## age       	0.055454   0.009884   5.611 2.02e-08 *** 
## genderfemale -0.244852   0.322817  -0.758  0.44816     
## bmi       	0.073879   0.023310   3.169  0.00153 **  
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
## (Dispersion parameter for binomial family taken to be 1) 
##  
## 	Null deviance: 326.02  on 383  degrees of freedom 
## Residual deviance: 281.62  on 380  degrees of freedom 
##   (19 observations deleted due to missingness) 
## AIC: 289.62 
##  
## Number of Fisher Scoring iterations: 5 
# 2. exponentiate the confidence intervals around the log odds for each predictor variable to obtain the odds 
exp(confint(full_model)) 
## Waiting for profiling to be done... 
##                	2.5 %  	97.5 % 
## (Intercept)  0.000179159 0.007892707 
## age      	1.037353246 1.078493131 
## genderfemale 0.414992991 1.478718788 
## bmi      	1.028712126 1.127738696 
# generate the McFadden pseudo R-square for this multiple regression 
############### 

# 1. run a null model 
null_model <- glm(dm ~ 1, family = binomial(link = logit)) 

# 2. check 
summary(null_model) 
##  
## Call: 
## glm(formula = dm ~ 1, family = binomial(link = logit)) 
##  
## Deviance Residuals:  
##	Min  	1Q  Median  	3Q 	Max   
## -0.578  -0.578  -0.578  -0.578   1.935   
##  
## Coefficients: 
##         	Estimate Std. Error z value Pr(>|z|)     
## (Intercept)  -1.7047 	0.1403  -12.15   <2e-16 *** 
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
## (Dispersion parameter for binomial family taken to be 1) 
##  
## 	Null deviance: 334.87  on 389  degrees of freedom 
## Residual deviance: 334.87  on 389  degrees of freedom 
##   (13 observations deleted due to missingness) 
## AIC: 336.87 
##  
## Number of Fisher Scoring iterations: 3 
# 3. calculate McFadden's R-square 
R2 <- 1-logLik(full_model)/logLik(null_model) 

# 4. print it 
R2 
## 'log Lik.' 0.1590084 (df=4) 
############### 


# generate the c-statistic of the multiple regression 
############### 

# 1. install package 
install.packages("DescTools") 
## Installing package into 'C:/Users/eg1012/Documents/R/win-library/3.5' 
## (as 'lib' is unspecified) 
## package 'DescTools' successfully unpacked and MD5 sums checked 
##  
## The downloaded binary packages are in 
##  C:\Users\eg1012\AppData\Local\Temp\RtmpgXykjR\downloaded_packages 
# 2. load package 
require(DescTools) 
## Loading required package: DescTools 
## Warning: package 'DescTools' was built under R version 3.5.1
# 3. calculate C-statistic 
Cstat(full_model) 
## [1] 0.7787709 
############### 


# perform Hosmer-Lemeshow test 
############### 

# 1. install package "ResourceSelection" 
install.packages("ResourceSelection") 
## Installing package into 'C:/Users/eg1012/Documents/R/win-library/3.5' 
## (as 'lib' is unspecified) 
## package 'ResourceSelection' successfully unpacked and MD5 sums checked 
##  
## The downloaded binary packages are in 
##  C:\Users\eg1012\AppData\Local\Temp\RtmpgXykjR\downloaded_packages 
# 2. load package 
require(ResourceSelection) 
## Loading required package: ResourceSelection 
## Warning: package 'ResourceSelection' was built under R version 3.5.1 
## ResourceSelection 0.3-2   2017-02-28 

# 3. run Hosmer-Lemeshow test 
HL <- hoslem.test(x = full_model$y, y = fitted(full_model), g = 10) 
HL  
##  
##  Hosmer and Lemeshow goodness of fit (GOF) test 
##  
## data:  full_model$y, fitted(full_model) 
## X-squared = 15.826, df = 8, p-value = 0.04494 
# 4. plot the observed vs expected number of cases for each of the 10 groups 
plot(HL$observed[,"y1"], HL$expected[,"yhat1"]) 

############### 




# create multiple logistic regression with several predictor variables 
############### 

# 1. generate model 
model <- glm(dm ~ age + bmi + chol + hdl + systolic + diastolic + gender + location + frame + insurance + smoking, family = binomial(link = logit)) 

# 2. check model results 
summary(model) 
##  
## Call: 
## glm(formula = dm ~ age + bmi + chol + hdl + systolic + diastolic +  
## 	gender + location + frame + insurance + smoking, family = binomial(link = logit)) 
##  
## Deviance Residuals:  
## 	Min   	1Q   Median   	3Q  	Max   
## -1.4845  -0.5506  -0.3577  -0.1948   2.6399   
##  
## Coefficients: 
##             	Estimate Std. Error z value Pr(>|z|)     
## (Intercept)    -7.793101   1.988365  -3.919 8.88e-05 *** 
## age         	0.052432   0.012882   4.070 4.70e-05 *** 
## bmi         	0.054568   0.028911   1.887  0.05910 .   
## chol        	0.010785   0.003589   3.005  0.00265 **  
## hdl        	-0.028312   0.010874  -2.604  0.00922 **  
## systolic    	0.005573   0.009455   0.589  0.55560     
## diastolic   	0.002992   0.016633   0.180  0.85723     
## genderfemale   -0.159954   0.381967  -0.419  0.67539     
## locationLouisa -0.255176   0.330092  -0.773  0.43950     
## framelarge  	0.262753   0.969952   0.271  0.78647     
## framemedium 	0.275534   0.967371   0.285  0.77578     
## framesmall  	0.550051   1.027476   0.535  0.59241     
## insurance1 	-0.273335   0.391855  -0.698  0.48546     
## insurance2 	-0.529986   0.406492  -1.304  0.19230     
## smoking2   	-0.158872   0.369045  -0.430  0.66684     
## smoking3   	-0.179607   0.508335  -0.353  0.72385     
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
## (Dispersion parameter for binomial family taken to be 1) 
##  
## 	Null deviance: 324.38  on 378  degrees of freedom 
## Residual deviance: 260.93  on 363  degrees of freedom 
##   (24 observations deleted due to missingness) 
## AIC: 292.93 
##  
## Number of Fisher Scoring iterations: 5 
# 3. test significance of each variable 
anova(model, test = "Chisq") 
## Analysis of Deviance Table 
##  
## Model: binomial, link: logit 
##  
## Response: dm 
##  
## Terms added sequentially (first to last) 
##  
##  
##       	Df Deviance Resid. Df Resid. Dev  Pr(>Chi)     
## NULL                    	378 	324.38               
## age    	1   33.737   	377 	290.64 6.309e-09 *** 
## bmi    	1	9.295   	376     281.34  0.002298 **  
## chol   	1	7.949   	375     273.40  0.004812 **  
## hdl    	1	9.043   	374     264.35  0.002638 **  
## systolic   1	0.555   	373     263.80  0.456200     
## diastolic  1	0.000   	372 	263.80  0.985146     
## gender 	1	0.146   	371     263.65  0.702584     
## location   1	0.456   	370     263.20  0.499712     
## frame  	3	0.404   	367     262.79  0.939510     
## insurance  2	1.647   	365 	261.15  0.438976     
## smoking	2	0.213   	363     260.93  0.899082     
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
############### 


# test for correlation between variables age and systolic blood pressure 
cor.test(systolic, age) # extremely significant 
##  
##  Pearson's product-moment correlation 
##  
## data:  systolic and age 
## t = 9.8342, df = 396, p-value < 2.2e-16 
## alternative hypothesis: true correlation is not equal to 0 
## 95 percent confidence interval: 
##  0.3604404 0.5187477 
## sample estimates: 
##       cor  
## 0.4430412 
# see whether systolic blood pressure becomes significant with age removed from the model 
############### 

# 1. generate model without age 
model <- glm(dm ~ bmi + chol + hdl + systolic + diastolic + gender + location + frame + insurance + smoking, family = binomial(link = logit)) 

# 2. test significance of each variable 
anova(model, test = "Chisq") # systolic is now significant, so there is a lot of variation explained by age and systolic 
## Analysis of Deviance Table 
##  
## Model: binomial, link: logit 
##  
## Response: dm 
##  
## Terms added sequentially (first to last) 
##  
##  
##       	Df Deviance Resid. Df Resid. Dev  Pr(>Chi)     
## NULL                    	378 	324.38               
## bmi    	1   7.0294   	377     317.35  0.008018 **  
## chol       1  16.0244   	376 	301.32 6.253e-05 *** 
## hdl    	1   7.5857   	375     293.74  0.005883 **  
## systolic   1   8.5837   	374     285.15  0.003392 **  
## diastolic  1   2.4307   	373 	282.72  0.118982     
## gender 	1   0.3899   	372     282.33  0.532328     
## location   1   0.7944   	371     281.54  0.372763     
## frame  	3   0.6076   	368     280.93  0.894682     
## insurance  2   0.6933   	366 	280.24  0.707066     
## smoking	2   0.8738   	364     279.36  0.646024     
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
############### 