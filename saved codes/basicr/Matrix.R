#matrix

A<- matrix(c(1:6), nrow = 2, ncol=3,byrow = TRUE)
A
dim(A)
A[1,3]
A[1,]
A[,2]
A[,c(2,3)]    A<-[,-1]
A[A%%2==0]
rownames(A)<-c("r1","r2")
colnames(A)<- c("c1","c2","c3")
A
t(A)
B <- matrix(c(7,8,9,19),nrow = 2,ncol=2,byrow = TRUE)
B
cbind(A,B)
rbind(A,B)
rbind(t(A),B)
A[1,3] <- 10
A
A[A<4] <- 0
A
A+5
c(A)

