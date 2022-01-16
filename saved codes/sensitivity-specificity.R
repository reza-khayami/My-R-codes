pred.summary <- function(data, Predict, Golden, Positive){
  true.positive <- sum(data[,Predict] == Positive & m[,Golden] == Positive)
  true.negative <- sum(data[,Predict] != Positive & m[,Golden] != Positive)
  false.positive <- sum(data[,Predict] == Positive & m[,Golden] != Positive)
  false.negative <- sum(data[,Predict] != Positive & m[,Golden] == Positive)
  accuracy <- round(true.positive + true.negative)/nrow(data), 4)
ppv <- round(true.positive/sum(data[,Predict] != Positive),4)
npv <- round(true.negative/sum(data[,Predict] == Positive),4)
sensitivity <- round(true.positive/sum(data[,Golden] == Positive),4)
specifity <- round(true.negative/sum(data[,Golden] != Positive),4)
list(true.positive = true.positive, true.negative = true.negative,
     false.positive = false.positive, false.negative = false.negative,
     accuracy = accuracy, ppv = ppv, npv = npv, sensitivity = sensitivity,
     specifity = specifity)
}
combn[1:5,2] #returns all dual combinations
?seq
