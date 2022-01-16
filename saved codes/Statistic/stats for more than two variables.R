rm(list = ls())
### multiple regression
data("USJudgeRatings")
reg1 <- lm(RTEN ~ CONT + INTG + DMNR + DILG + CFMG + DECI + PREP + FAMI + ORAL+ WRIT + PHYS, data = USJudgeRatings)
reg1
summary(reg1)
# More detail
anova(reg1)
coef(reg1)
confint(reg1)
resid(reg1)
hist(residuals(reg1))
#Possibility of stepwise variable selection
#(backwards and forwards); exercise caution!

#backwards stepwise regression
#repeating the first regression model, which contains
#all of the predictor variables and serves as the starting point
reg1 <- lm(RTEN ~ CONT + INTG + DMNR + DILG + CFMG + DECI + PREP + FAMI + ORAL+ WRIT + PHYS, data = USJudgeRatings)
regb <- step(reg1, 
             direction = "backward",
             trace = 0) # Don't print the steps
summary(regb)
#forward stepwise regression
#Start with model that has nothing but a constant
reg0 <- lm(RTEN ~ 1, data = USJudgeRatings) #Minimal model
reg0
regf <- step(reg0, #start with minimal model
             direction = "forward",
             scope = (~ CONT + INTG + DMNR + DILG + CFMG + DECI + PREP + FAMI + ORAL+ WRIT + PHYS), data = USJudgeRatings,
             trace = 0)
# more information "rms" package
### comapring means with a two factor anova
data("warpbreaks")
boxplot(breaks ~ wool*tension, data = warpbreaks)
#Model with interaction
aov1 <- aov(breaks ~ wool + tension + wool : tension,
           # or : wool*tension
           data = warpbreaks)
summary(aov1)
#additional info on model
model.tables(aov1)
model.tables(aov1, type = "means")
model.tables(aov1, type = "effects") #"effects is defualt"
#Post-hic test
TukeyHSD(aov1)
# Conducting a cluster analysis
data(mtcars)
mtcars1 <- mtcars[, c(1:4, 9:11)]
#Three major kinds of clustering:
# 1.split into set numbers of clusters (e.g., kmeans)
# 2. hierarchical : start separate and combine
# 3.dividing: start with a singke group and split
##hierarchical clustering
# Need distance matrix (dissmilarity matrix)
d <- dist(mtcars1)
# use distance matrix for clustering
c <- hclust(d)
c
#plot dendogram of clusters
plot(c)
# put observasions in groups
# need to specify either k = groups or h = height
g3 <- cutree(c, k = 3)
# Or do several levels of groups at once
# "gm" = "groups/multiple"
gm <- cutree(c, k = 2:5) # or k = c(2, 4)
# Draw boxes around clusters
rect.hclust(c, k = 2, border = "gray")
rect.hclust(c, k = 3, border = "blue")
rect.hclust(c, k = 4, border = "Red")
rect.hclust(c, k = 5, border = "green")
##k-means clustering
km <- kmeans(mtcars1,3)
km
# graph based on k-means
require(cluster)
clusplot(mtcars1, #data frame
         km$cluster, # cluster data
         color = TRUE, #color
         #shade = TRUE, # lines in clusters
         lines = 3, # lines connecting centroids
         labels = 2) # labels clusters and cases
### Conducting a principal components factor analysis
# components vs factor model : treatment of the variacnces for each item
# components are weighted composites of observed variables
#variables are weighed composites of the factor in the factor model
mtcars1 <- mtcars[, c(1:4, 6:7, 9:11)]
#Pricnipal components model using default method
# If using entire data frame: 
pc <- prcomp(mtcars1,
             center = TRUE, #center means to 0 (optional)
             scale = TRUE) #sets unit variance
# Or specify variables:
#pc <- prcomp(~ mpg + cyl + .... + carb, data = mtcars, scale TRUE)
?prcomp # generally preferred
?princomp # slightly different similar to S
summary(pc)
#screeplot
plot(pc)
#see how cases load on PCs
predict(pc)
#biplot
biplot(pc)

# Factor Analysis
# gives chi square test that number of factors
# is sufficient to match data (want p > .05)
# Also gives uniqueness values for valriables,
#variable loadings on factors, and variance statistics.
factanal(mtcars1, 1)
factanal(mtcars1, 2)
factanal(mtcars1, 3)
factanal(mtcars1, 4) # first p>0.05 with four factors the observations do not significantly differ from the models
factanal(mtcars1, 5)
###excersice
# Load data
scd <- read.csv("~/Desktop/R/StateClusterData.csv", header = TRUE)
rownames(scd) <- scd[,1]  # Use state names for row names
scd[,1] <- NULL  # Remove state names as variable
scd[1:5, ]

# We'll use hierarchical clustering
d <- dist(scd)  # Distance matrix
c <- hclust(d)  # Get clusters
plot(c)  # Dendrogram

# Draw boxes around clusters
rect.hclust(c, k = 2, border = "gray")
rect.hclust(c, k = 3, border = "blue")
rect.hclust(c, k = 4, border = "green4")
rect.hclust(c, k = 5, border = "darkred")

rm(list = ls())  # Clean up