dimnames(Data)
insurance: 0=none, 1=government, 2=private

chol <- Data['chol'] # cholesterol is continuous, so it's easy
gender <- as.factor(Data[,"gender"]) # but gender isn't.
dm <- as.factor(Data[,"dm"]) # neither is dm

t <- table(gender) # store the tabulation for further manipulation
addmargins(t)
round(prop.table(t),digits=3)
round(100*prop.table(t),digits=1) # get %s rounded to 1dp

dm2 <- factor(dm, exclude=NULL) # make new factor from the old one

age <- ifelse(Data$age < 45, "U45",
              ifelse(Data$age >= 45 & Data$age < 65, "45-64",
                     ifelse(Data$age >= 65 & Data$age < 75, "65-74",
                            ifelse(Data$age >= 75, "75 over", NA))))
a <- table(age,gender)
round(prop.table(a) *100)


m <- glm(dm ~ 1, family=binomial (link=logit))
summary(m)

table(m$y)

m <- glm(dm ~ Data$age, family=binomial (link=logit))
summary(m)


# It's straightforward to include age as a single term in the model, but remember what I said in the video about assuming a linear relation with the outcome? More precisely, this assumes that the relation between age and the log odds of having diabetes is linear (more on this in detail in the next section). Is that reasonable? The easiest way is just to plot one against the other.
# create a cross tabulation of age and diabetes status  
dm_by_age <- table(Data$age, dm) 

# output the frequencies of diabetes status by age 
freq_table <- prop.table(dm_by_age, margin = 1) 

# calculate the odds of having diabetes 
odds <- freq_table[, "yes"]/freq_table[, "no"] 

# calculate the log odds 
logodds <- log(odds) 

# plot the ages found in the sample against the log odds of having diabetes 
plot(rownames(freq_table), logodds) 


#bmi 
height <- Data[,"height"] 
weight <- Data[,"weight"] 
height.si <- height*0.0254 
weight.si <- weight*0.453592 
bmi <- weight.si/height.si^2 


# categorising BMI 

bmi_categorised <- ifelse(bmi < 18.5, "underweight", 
                          ifelse(bmi >= 18.5 & bmi <= 25, "normal", 
                                 ifelse(bmi > 25 & bmi <= 30, "overweight", 
                                        ifelse(bmi > 30, "obese", NA))))

glmfit <- glm(dm~ age + gender + bmi_categorised, family = binomial(link = "logit"))
summary(glmfit)
plot(glmfit)

summary(glm(dm~gender,family = binomial(link = "logit")))
exp(glmfit$coefficients)
exp(confint(glmfit))   

