#dataframe

products <- c("p1","p2","p3")
products
unitprice <- c(20,15,40)
unitprice
monthlydemands <- c(1500,2000,850)
monthlydemands
df<-data.frame("products","unitprice","monthlydemands")
df

x <- mtcars
x
class(x)
dim(x)
head(x)
tail(x)
x$mpg
x[x$mpg<20,]
