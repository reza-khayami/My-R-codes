h <- c(1,2,NA,4,5)
which(is.na(h))
mean(h, na.rm=T)
h2 <- h
h2[is.na(h)] <- 0
h
#second way
h3 <- ifelse(is.na(h), 0, h)
#packages
mi
mice
imputation