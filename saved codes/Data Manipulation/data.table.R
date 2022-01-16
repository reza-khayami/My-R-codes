rm(list = ls())
install.packages("data.table")
require(data.table)

url <- "https://raw.githubusercontent.com/selva86/datasets/master/mtcars.csv"
#setDT() for big files data.table for small
#setDT() is not a replacement for data.table(). It's a more efficient replacement for as.data.table() which can be used with certain types of objects.

#mydata <- as.data.table(mydata) will copy the object behind mydata, convert the copy to a data.table, then change the mydata symbol to point to the copy.
#setDT(mydata) will change the object behind mydata to a data.table. No copying is done.
#So what's a realistic situation to use setDT()? When you can't control the class of the original data. For instance, most packages for working with databases give data.frame output
mt <- fread(url)
mt
mtcars$carname <- rownames(mtcars)
df <- mtcars
mtcars_dt <- data.table(df)
setDT(df) #datatable to dataframe
setDF(df)
#data.frame ==> mtcars[mtcars$cyl == 6 & mtcars$gear == 4]
#data.table ==> mtcars_dt[cyl==6 & gear==4,]
#dataframe mtcars[, 1]  === data table mtcars_dt[, 1, with = F] 
#or mtcars_dt[, mpg] 
#multple columns
mtcars_dt[, list(mpg, cyl, gear)] 
# or mtcars_dt[, .(mpg, cyl, gear)]
a <- airquality
setDT(a)
a2 <- a[!is.na(Ozone), ]
a2[, list(Ozone, Solar.R, Wind, Temp)]

###syntax
#i = where, j = select,  by = group by

mtcars$carname <- rownames(mtcars)
mtcars_dt <- data.table(mtcars)
mtcars_dt$cyl_gear <- mtcars_dt$cyl + mtcars_dt$gear
mtcars_dt
#OR
mtcars_dt[, cyl_gear2 := cyl + gear]
# multiple new columns :
mtcars_dt[, ':='(cyl_gear3 = cyl * gear, cyl_gear4 = cyl - gear)]
# if you want to add and select only new ones :
mtcars_dt[, .(cyl_gear3 = cyl * gear, cyl_gear4 = cyl - gear)]
#name of column is in another charachter vector
myvar <- c('var1')
mtcars_dt[, myvar:=1] #didnt want this
mtcars_dt[, c(myvar):= 1]
mtcars_dt[, c("cyl_gear2", "cyl_gear3", "cyl_gear4", "cyl_gear", "myvar", "var1") := NULL] #delete
mtcars_d <- mtcars_dt[, scaled_mpg := round(mpg/mean(mpg), 2)]
###chaining, functions and .SD
mtcars$carname <- rownames(mtcars)
mtcars_dt <- data.table(mtcars)
dt1 <- mtcars_dt[, .(mean_mpg = mean(mpg), mean_disp = mean(disp),
mean_wt = mean(wt), mean_qsec = mean(qsec)), by = cyl]
output <- dt1[order(cyl), ]
#chain
output <- mtcars_dt[, .(mean_mpg = mean(mpg),
                                 mean_disp = mean(disp),
                                 mean_wt = mean(wt),
                                 mean_qsec = mean(qsec)),
                             by = cyl][order(cyl), ]
 mtcars_dt[, .SD, by = cyl]
 output <- mtcars_dt[, lapply(.SD[, 1:10, with = F], mean), by = cyl] #OR
output <- mtcars_dt[, lapply(.SD, mean), by = cyl,
                    .SDcol = c("mpg", "disp", "hp", "drat",
                               "wt", "qsec")]
#challange
set.seed(100)
rand <- sample(1:10, 100000, replace = T)
mymat <- matrix(rand, ncol = 5)
mymat <- data.table(mymat)
ss <- function(x){(x['V1']^2)+(x['V2']^2)+(x['V3']^2)+(x['V4']^2)+(x['V5']^2)}
mymat[, ss := apply( .SD, 1, ss)]
mtcars$carname <- rownames(mtcars)
mtcars_dt <- data.table(mtcars, key = "carname" ) # sorts by key
key(mtcars_dt) #or
mtcars_dt <- data.table(mtcars)
setkey(mtcars_dt, carname)
#filter the table
mtcars_dt["Merc 230"]
#joining
dt1 <- mtcars_dt[1:10, .(carname, mpg, cyl)]
dt2 <- mtcars_dt[6:15, .(carname, gear)]
dt1[dt2] #caues both have the same key
dt2[dt1]
dt2[dt1, nomatch = 0] #retain only the common rows inner join
k <- union(dt2$carname, dt1$carname) # outer join 
dt1[dt2[.(k)]]
# you can just use merge
merge(dt1, dt2, all = T) #FULL
merge(dt1, dt2, all = F) #inner
merge(dt1, dt2, all.x = T) #left
merge(dt1, dt2, all.y = T) #right
#multiple keys 
#safer to enclose them in a . or j function when they are numeric
setkey(mtcars_dt, cyl, gear)
key(mtcars_dt) 
mtcars_dt[c(8, 3)] #wrong data
mtcars_dt[.(8, 3)]
mtcars_dt[J(8, 3)]
#remove keys
setkey(mtcars_dt, NULL)
mtcars_dt
#aggregating
output <- mtcars_dt[, .(mean_mpg = mean(mpg),
                        mean_disp = mean(disp),
                        mean_wt = mean(wt),
                        mean_qsec = mean(qsec)), by = cyl][order(cyl),]
#or
output <- mtcars_dt[, .(mean_mpg = mean(mpg),
                        mean_disp = mean(disp),
                        mean_wt = mean(wt),
                        mean_qsec = mean(qsec)), keyby = cyl]
output

#set function [works in data frames as well]
M <- matrix(1, nrow = 10000, ncol = 2)
dt <- as.data.table(M)
#time fir doing a for loop
system.time({
  for(i in 1:nrow(dt)){dt[i, V1 := 0]
  }
  })
dt <- as.data.table(M)
system.time({
  for(i in 1:nrow(dt)){set(dt, i, 2L, 0)
  }
})
#challange
M <- matrix(1, nrow = 10000, ncol = 10000)
df <- as.data.frame(M)
args(set)
for(i in 1:nrow(df)){
  set(df, i, i, i)
}

df[1:5, 1:5]
