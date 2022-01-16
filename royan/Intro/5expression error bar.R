###### expression
#SD DE HLXb9
sd(x[4:6, 1])
se <- function(x){
  sd(x)/sqrt(length(x))
}
# because the function is just one line you can write it as follows:
# se <- function(x) sd(x)/sqrt(length(x))
#Another way:
#se <- function(x){
#e = (sd(x)/sqrt(length(x)))
#return(e)
#}
a <- se(x[4:6, 1])
length(x[4:6, 1])
options(digits = 4 )
#mean
head(x)
x$Sample <- rownames(x)
x$Sample <- substr(x$Sample, 1, nchar(x$Sample)-2)
#nchar returns lenghth of the charachter
#nchar("DE") = 2
#substr("Hello", 4, 5)
x.m <- aggregate(. ~ Sample, data = x, FUN = mean)
# . means everything
head(x.m)
ggplot(x.m, aes(Sample, HLXb9, fill = Sample))+ geom_bar(stat = "identity")
#### khodam error bar keshidam!
x.d <- aggregate(. ~ Sample, data = x, FUN = sd)
x.d1 <- x.d[,-1]
x.d2 <- x.d1/sqrt(length(x.d1))
ggplot(x.m, aes(Sample, HLXb9, fill = Sample))+ geom_bar(stat = "identity") + geom_errorbar(aes( ymin = x.m[,2]-1.96*x.d2[, 1],                                    ymax = x.m[,2]+1.96*x.d2[, 1]),
                                                                                            width=.2, position=position_dodge(.9))