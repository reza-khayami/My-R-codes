#-------------
library(lattice)
library(mice)
library(SNPassoc)
library(mice) # Data imputation
library(dplyr) # Data manipulation
library(magrittr) # Flexible piping in R
library(purrr) # Flexible functional programming
set.seed(123)

#prepare the data rame---------------------------------
##Get an overview of the data by the summary() command:
temp3$Status <- droplevels(temp3$Status, exclude = "p-value")
caseT <- temp3[501:919,]
conT <- temp3[1:500,]
summary(caseT)
summary(conT)

a <- c(2,6,7,8,12)
caseT2 <- caseT[,-a]
b <- c(a,10)
conT2 <- conT[,-b]


#Inspect the missing data pattern---------------
md.pattern(caseT2)


imp <- mice(caseT2,
            m = 5,
            method = c("","pmm", "logreg", "logreg","logreg", "logreg", ""),
            maxit = 20,
            seed=123)
imp$pred
dis <- complete(imp)
imp <- mice(caseT2,
            m = 5,
            quickpred(caseT2, mincor = 0.3),
            maxit = 20,
            seed=123)
#Inspect the convergence of the algorithm
#should intermingle
plot(imp)
#Further diagnostic checking.
stripplot(imp)
stripplot(imp, Age~.imp, pch=20, cex=2)

fit <- with(imp, glm(Age ~ Addiction))
fit

summary(fit$analyses[[5]]) #for every iteration

pool.fit <- pool(fit)
summary(pool.fit)
pool.fit$pooled

R <- is.na(caseT2$Age) 
histogram(~ Age, data = caseT2)
histogram(~Addiction|R, data=caseT2)

#compare with incomplete
summary(caseT2)
summary(complete(imp))

#It is best to do a Fisher transformation before pooling the correlation estimates - and a backtransformation afterwards. Therefore we define the following two functions that allow us to transform and backtransform any value:
fisher.trans <- function(x) 1/2 * log((1 + x) / (1 - x))
fisher.backtrans <- function(x) (exp(2 * x) - 1) / (exp(2 * x) + 1)

cor <- imp %>%
  mice::complete("all") %>%
  map(select, -phb, -gen, -reg) %>%  
  map(stats::cor) %>%
  map(fisher.trans)
cor

# The object cor is a list over the m imputations where each listed index is a correlation matrix. To calculate the average over the correlation matrices, we can add the m listed indices and divide them by m:
  
  cor.rect <- Reduce("+", cor) / length(cor) # m is equal to the length of the list
cor.rect <- fisher.backtrans(cor.rect)



scope <- list(upper = ~  Age + Sex + Addiction + FHC + TL,
              lower = ~ 1)
expr <- expression(f1 <- lm(Genotype ~ 1),
                   f2 <- step(f1, 
                              scope = scope, 
                              direction = "forward",
                              trace = 0
                   ))
fit <- with(imp, expr)
