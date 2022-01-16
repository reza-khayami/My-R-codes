rm(list=ls())
df1 <- mtcars[mtcars$hp > 150,]
df2 <- aggregate( . ~ cyl, data = df1, FUN = mean)
install.packages("magrittr")
?magrittr
require("magrittr")
rnorm(10) %T>% plot(main = "rnorm") %>% sum
  
data(mtcars)
mtcars %>% subset (cyl == 6) %$% cor(mpg, wt) 
a <- mtcars[mtcars$cyl == 6,]
cor(a$mpg, a$wt)
#dplyr
require("hflights")
install.packages("dplyr")
require(dplyr)
hf_tbl <- as.tbl(hflights)
hf_tbl
?tbl
dplyr::glimpse(hf_tbl)
filter(hf_tbl, Month == 1, Year == 2011)
hf_tbl %>% filter(Month == 1) %>% slice(1:5)
select(hf_tbl, Year, DayOfWeek)
select(hf_tbl, Year:DayOfWeek)
select(hf_tbl, 1,2)
select(hf_tbl, c(1,2))
select(hf_tbl, ends_with("Time"))
select(hf_tbl, contains("Time"))
select(hf_tbl, starts_with("Day"))
hf_tbl_1 <- rename(hf_tbl, mth = Month)
hf_tbl_2 <- mutate(hf_tbl, yearmonth = paste0(Year, Month) %>% as.numeric)
hf_tbl_3 <- transmute(hf_tbl, yearmonth = paste0(Year, Month) %>% as.numeric)
arrange(hf_tbl, Year, Month) #sort
arrange(hf_tbl, Year, desc(Month)) #sort by decending order

#challange
data(mtcars)
mtcars_tbl <- as.tbl(mtcars)
mtcars_tbl <- mutate(mtcars_tbl, make = rownames(mtcars))
mtcars_tbl <- mutate(mtcars_tbl, kpl = mpg*0.425)      
mtcars_tbl <- mutate(mtcars_tbl,kpl_int = round(kpl))
mtcars_tbl <- arrange(mtcars_tbl, kpl)
mtcars_tbl <- filter(mtcars_tbl, vs== 1)
#Answer
data(mtcars)
mtcars_tbl <- mtcars %>% mutate( make = rownames(mtcars)) %>%
   mutate(kpl = 0.425 * mpg, kpl_int = round(kpl)) %>% arrange (kpl) %>% filter(vs == 1)
###Aggregaton and Special Functions
  require("hflights")
require("dplyr")
hf_tbl <- as.tbl(hflights)
hf_tbl_gp <- group_by(hf_tbl, UniqueCarrier)
summarise(hf_tbl_gp, mean_dist = mean (Distance), mea_delay = mean(DepDelay, na.rm = T))
unique_carrier_agg <- hflights %>% group_by(UniqueCarrier) %>% summarise(delay = mean(DepDelay, na.rm = T), num_obs =n(#gets the n of obs by the group variable
  ), num_distinct_obs = n_distinct(DepDelay), first_obs = first(DepDelay), second_obs = nth(DepDelay , 2), last_obs = last (DepDelay) )
hf_tbl_gp <- ungroup(hf_tbl_gp)
#challange
hf_tbl <- as.tbl(hflights)
Origin_agg <- hflights %>% group_by(Origin) %>% summarise(delay = mean(DepDelay, na.rm = T), num_obs =n(#gets the n of obs by the group variable
), num_distinct_obs = n_distinct(DepDelay), first_obs = first(DepDelay), second_obs = nth(DepDelay , 2), last_obs = last (DepDelay) )
###Two table verbs
mtcars$carname <- mtcars %>% rownames
mtcars_tbl <- mtcars %>% as.tbl
mtcars_tbl_1 <- mtcars_tbl %>% select(carname, mpg:drat)
mtcars_tbl_2 <- mtcars_tbl %>% select(carname, wt:am)%>% sample_frac(0.4)
mtcars_tbl_3 <- rename(mtcars_tbl_2, car = carname)
mtcars_tbl_1 %>% left_join(mtcars_tbl_2)
mtcars_tbl_1 %>% left_join(mtcars_tbl_3, c("carname" = "car"))
right_join()   inner_join()   full_join()
###Working with databases
