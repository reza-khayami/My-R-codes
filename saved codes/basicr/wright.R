#wright
wright.csv
wright.table
writeLines()
#example
products <- c("p1","p2","p3")
price <- c(20,15,40)
demand <- c(1500,2000,850)
df <- data.frame("products","price","demand")
write.table(df,"",sep = " ",row.names = FALSE,col.names=TRUE)