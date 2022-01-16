install.packages("hflights")
data("hflights" , package = "hflights")
head(hflights)
hflights$delayed <- ifelse(hflights$DepDelay > 0, 1, 0) 
  carrier_delay_count <- table(hflights$UniqueCarrier, hflights$delayed)
  carrier_delay_count <- as.data.frame.matrix(carrier_delay_count)
  carrier_delay_count$perc_delayed <- carrier_delay_count$`1`/as.numeric(carrier_delay_count$`1` + carrier_delay_count$`0`)
  max(carrier_delay_count)
  rownames(carrier_delay_count)[which.max(carrier_delay_count$`1`)]
  carrier_delay_count[carrier_delay_count == max(carrier_delay_count)]
  b <- (c(hflights$UniqueCarrier ~ hflights$DepTime))
  head(b)
h <- as.data.frame(cbind(hflights$UniqueCarrier, hflights$DepTime))
g <- hflights[,c(5,7)]
?aggregate()
e <- rep(g, g$UniqueCarrier )
b <- aggregate(hflights$DayofMonth ~ hflights$Month   , FUN = sum)
View(a)      
c <- aggregate(hflights$DayofMonth, by=list(hflights$Month, hflights$DayOfWeek), FUN=sum)
install.packages("lubridate")
require(lubridate)
hflights$Date <- paste(hflights$Year, hflights$Month, 
                       hflights$DayofMonth, sep = "-")
?paste
b <- paste(hflights$Year, hflights$Month)
head(b)
hflights$Date <- ymd(hflights$Date, tz= "UTC")
departures <- table(hflights$Date)
departures[departures == max(departures)]
departures[which.max(departures)]
a <- hflights[,c(7,9)]
carrirer_tailnum <- (as.data.frame.matrix(table(hflights$UniqueCarrier, hflights$TailNum)))
fleetcount <- as.data.frame(ifelse(carrirer_tailnum > 0, 1, 0) )
rowSums(fleetcount)
a <- as.data.frame.matrix(table(hflights$UniqueCarrier, hflights$TailNum))
b <- as.data.frame(ifelse(a >0 , 1, 0))
rowSums(b)
b[b == max(b)]
a1<- aggregate(cbind(ActualElapsedTime, AirTime, ArrDelay,DepDelay, Distance, TaxiIn, TaxiOut) ~ FlightNum, data = hflights, FUN = function(x)round(mean(x, na.rm = T), 2 ))
a2 <- aggregate( Year ~ FlightNum,data = hflights, FUN = length)
a3 <- aggregate(Cancelled ~ FlightNum,data=hflights[hflights$Cancelled==1,], FUN = length)
require("reshape2")
        gb <- dcast(UniqueCarrier ~ Month, data = hflights, value.var="DepDelay", fun.aggregate = function(x)round(mean(x, na.rm=T),2))
   ccgb <- gb[,-1]
