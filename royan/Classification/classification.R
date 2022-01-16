library(reshape2)
library(ggplot2)
library(MASS)
library(e1071)
install.packages("e1071")
irism <- melt(iris)
ggplot(irism, aes(variable, value, fill = Species)) +
  geom_boxplot()

N <- nrow(iris)
#Randomly distribute half of data in two data sets one for train
trainset <- sample(seq(N), N/2)

# if we want to get 25 for each group:
trainset <- sapply(unique(iris$Species), function(s) sample(
  which(iris$Species == s), 25))
trainset <- as.numeric(trainset)
testset <- setdiff(seq(N), trainset)
length(testset)
any(testset %in% trainset)
intersect(trainset, testset)

model.lda <- lda(Species ~ ., iris[trainset, ])
model.lda

predict.lda <- predict(model.lda, iris[testset, 1:4]) 
# we should not give Species the model should guess it!
#although the model won't cheat anyway!
predict.lda
names(predict.lda)
predict.lda$class
#test to see how many Species the model has  succesfully differentiated
sum(predict.lda$class == iris[testset, "Species"])
table(predict.lda$class, iris[testset, "Species"])
which(predict.lda$class != iris[testset, "Species"])

#measuring accuracy
cat("Accuracy is: ", sum(predict.lda$class == iris[testset, "Species"])/length(testset))
    

#SVM

model.svm <- svm(Species ~ ., iris[trainset, ])
predict.svm <- predict(model.svm, iris[testset, ]) 
cat("Accuracy is: ", sum(predict.svm == iris[testset, "Species"])/length(testset))

#better thing to do is to run the model lets say 1000 times and then get accuracy mean

find.accuracy <- function(train, test, algorithm) {
  model <- algorithm(Species ~ ., train)
  pr <- predict(model, test)
  sum(pr == test$Species)/ nrow(test)
}

find.accuracy(iris[trainset,], iris[testset,], svm) #remember lda has class and this function returns 0
debug(find.accuracy)
undebug(find.accuracy)

data.partition <- function(data, algorithm) {
  N = nrow(data)
  trainset <- sample(seq(N), N/2)
  testset <- setdiff(seq(N), trainset)
  find.accuracy(data[trainset,], data[testset,], algorithm)
}
result <- sapply(1:1000, function(i) data.partition(iris, svm))
mean(result)
hist(result)                 

#we can run both algorithms :
result <- lapply(c(svm, randomForest), function(alg) sapply(1:1000, function(i) data.partition(iris, alg)))
result <- do.call(cbind, result) #do.call runs a function on a list 