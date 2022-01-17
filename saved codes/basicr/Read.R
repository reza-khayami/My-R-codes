#read

scan()
readline()
scan(what=" ")
#exampel
my.name <- readline(prompt = "enter name:")
my.age <- readline(prompt = "enter age:")
my.age <- as.integer (my.age)
class(my.age)
print(paste("hi",my.name,"next year you will be",my.age +1,"years old"))


#read
getwd()
setwd("Desktop/")
read.csv()
data <- read.csv("Data.csv",header =T)
head(data)
read.table()
data <- read.table("Data.csv",sep= ",",header=T)
data
data1<- read.table("Desktop/Data.csv",sep= ",",header=T,nrow=50)
data2 <- read.table("Data.txt",sep=",",header=T)
data2
d <- read.delim("Desktop/Data.csv",sep=",")
d
