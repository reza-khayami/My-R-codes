#How to check your data in R

# The most common types of data are: 
  
# Numeric 	Decimal values 
# Integer 	Numeric value with no fraction 
# Character 	Text/String variables, enclosed by apostrophes 
# Logical 	Objects (TRUE or FALSE) that encode a logic 
# Factor 	Categorical/nominal variables, with levels of data 

COPD <- read.csv("COPD.csv")

class(COPD$MWT1Best)

class(COPD$copd)
# The output shows that both variables are saved as integers. However, only MWT1Best should be saved as an integer.  

# We therefore need to change the copd variable from an integer to a factor variable. We do this using the command factor()
COPD$copd <- factor(COPD$copd)
class(COPD$copd)
# To visualise the structure of the data in that variable, you can use the str() function: 

# Now, if you run your regression, you will get an output for three of the levels of the copd variable:  
  

lr1 <- lm(MWT1Best~ copd, data = COPD)
summary(lr1)

# there are only coefficients for three of the levels because the coefficients represent the comparison with the fourth level, which is often referred to as the reference group. 

##Changing the reference category of a categorical variable -----

#When running the above regression of walking distance on COPD severity, the linear regression model automatically uses COPD level 1 (i.e. 'mild') as the reference category:  

# If you wanted to use the 'severe' group as a reference category for COPD severity (i.e. level 3), you can use the relevel() function: 
  
  
COPD$copd <- relevel(COPD$copd, ref=3) 

# Create new variables from old ones ------
# You now want to create a new variable, which indicates the presence of at least one comorbidity or complete absence of comorbidities, based on the responses to the variables: Diabetes, muscular, hypertension, AtrialFib, and IHD. 
comorbid <- length(COPD$Diabetes) 

# comorbid will have a value of 1 if Diabetes = 1, OR muscular = 1, OR hypertension = 1, OR AtrialFib = 1, OR IHD = 1. 
# comorbid will have a value of 0 if ALL of the above variables = 0. 
comorbid[COPD$Diabetes == 1 | COPD$muscular == 1 | COPD$hypertension == 1 | COPD$IHD == 1 | COPD$AtrialFib == 1] <- 1
comorbid[is.na(comorbid)] <- 0
comorbid <- factor(comorbid)
COPD$comorbid <- comorbid

# Practice with R: Run a Good Practice Analysis----------------
#1. Inspect the dataset for missing values and outliers 
library('Hmisc')
describe(COPD)

#1a. For categorical variables, you can tabulate the data using the CrossTable() function from the 'gmodels' package, then use the sum(is.na()) functions to check for missing values. 
library('gmodels')

CrossTable(COPD$copd)
#OR
table(COPD$copd)
round(proportions(table(COPD$copd)), 3)

#1b. For continuous variables, such as MWT1Best, use the summary() command, which will allow you to look at the mean, median, minimum, maximum, 1st and 3rd quartiles, and the number of missing values (NAs):
summary(COPD)
hist(COPD$AGE) 

#If you do this for all the continuous variables in the dataset, you'll find at least one variable that has an extreme value. So you'll need to make a decision whether this value is an impossible value, maybe due to a coding error, in which case the value may need to be excluded, or whether you think this is it just an unusual value that should be left. In practice, if you saw a potential error, you would ask the researcher to go back to the original data and check, as part of the data cleaning prior to analysis.  

#2. Examine the relationship between your candidate predictor variables

#2a. At the start of this course, you learned how to calculate a correlation coefficient using the cor.test() command. To see the pairwise correlation coefficient only for continuous variables, you can use the cor() command. The output of this command will be Pearson's correlation coefficient by default. If you want Spearman's correlation coefficient, you can specify this by adding method = 'spearman' in the cor() command parentheses.  

mydata <- COPD[, c('AGE', 'PackHistory', 'FEV1', 'FEV1PRED', 'FVC', 'CAT', 'HAD', 'SGRQ')]
cor.matrix <- cor(mydata)
round(cor.matrix, 2)
#corplot
pairs(~AGE+ PackHistory+ FEV1+ FEV1PRED+ FVC+ CAT+ HAD+ SGRQ, data = mydata)

#2b. To examine associations between categorical variables, you can use cross tabulations. For this, you will again use the CrossTable() function from the 'gmodels' package, but adding your variables in the following format: CrossTable(mydata$myrowvar, mydata$mycolvar). 

  CrossTable(COPD$hypertension, COPD$IHD)

  #This cross tabulation shows you that 8 patients had IHD only, and 11 patients had hypertension only, one person had both hypertension and IHD and 81 had neither. So, you would not be concerned if you wanted to include both of these in the model.
  
#3. Fit a simple linear regression model 
  
  #And finally, it's useful to assess the relationship for each of variable in turn with the outcome. You can do this by fitting a regression model with just a single predictor in. Doing this allows an opportunity to spot anything unusual that may due to errors in either the data or coding of the variable, and also allows you to anticipate what you might expect to happen when you fit the multivariable model.    
  
  #outcome ~ each predictor
  
  model_name <- lm(outcome ~ predictor, data =dataframe) 
  
  summary(model_name) 
  
  confint(model_name)  
  
  Summary 
  
  mlr1 <- lm(MWT1Best~ FEV1 + AGE + factor(gender) + factor(COPDSEVERITY) + factor(comorbid), data = mydata)

  
  # 
  # The first thing is examining variable distributions using summary statistics, tabulations and graphs. 
  # 
  # The next stage is examining the relationship between candidate predictors using cross tabulations and correlations. 
  # 
  # And finally, get a feel for the relationships between each of the candidate predictors and the outcome by fitting a regression model for each variable in turn.  
  # 
  # This strategy will generate a fair amount of output for you to look through. My approach would be: 
  #   
  #   to set the statistical code up in advance 
  # run all the analyses
  # then grab a coffee - take some time to look over and absorb the results. 
  
  #check collinearity 
  library('mctest')
  imcdiag(model.matrix(mlr1)[,-1], mlr1$model[1], method = 'VIF')  