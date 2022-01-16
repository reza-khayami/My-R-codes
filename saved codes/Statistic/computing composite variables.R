rm(list = ls())
#create variable rn1 with 1 million random normal values 
#will vary from one time to another
rn1 <- rnorm(1000000)
hist(rn1)
summary(rn1)
#create variable rn2 with 1 million random normal values 
rn2 <- rnorm(1000000)
hist(rn2)
summary(rn2)
#average scores across two variables
rn.mean <- (rn1+ rn2)/2
hist(rn.mean)
#multiply scores across two variables
rn.prod <- rn1*rn2
hist(rn.prod)
summary(rn.prod)
#kurtosis comparisans
kurtosi(rn1)
