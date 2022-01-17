#basic
getwd()
setwd()
#vector

b <- c(1,2,3,4,5)
d<- 1:10
d
d <-  seq(1,10,by=2)
d
d<- seq(1,10,length=4)
d
d1<- rep(5,5)
d1
length(d1)
d1[4]
d1[c(3,4)]

sort(d,decreasing = TRUE)

q <- 2*1:5
q
w <- 2*1:5-1
w

z<- 2.143
z
t<- -1.345
t
abs(t)
sign(z)
sqrt()
round(z)
log()
sin()
factorial()
r<-1:5
r
cumsum(r)
cumprod(r)


v1<- c(1,4,5,7,9,13)
v1
v2<- c(10,37,53,68,81,118)
v2
min()
mean()
max()


var()
sd(v2)
cor(v1,v2)
median(v1)
sum(v1)

v1==5
which(v1==5)
any(v1>7)
all(v1<=7)
!v1==5
& 
|
help(c)
x <- 1:120
x
x[x%%2==0]

