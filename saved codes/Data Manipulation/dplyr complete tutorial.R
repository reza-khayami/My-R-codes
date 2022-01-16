# https://rpubs.com/justmarkham/dplyr-tutorial


# load packages----------------------------------------------------------

suppressMessages(library(dplyr))
library(hflights)


# explore data-----------------------------------------------------------
data(hflights)
head(hflights)


# convert to local data frame--------------------------------------------

flights <- tbl_df(hflights)     #tbl_df creates a "local data frame"
                                # Local data frame is simply a wrapper                                  for a data frame that prints nicely

# printing only shows 10 rows and as many columns as can fit on your screen
flights

# you can specify that you want to see more rows
print(flights, n = 20)

# convert to a normal data frame to see all of the columns
data.frame(head(flights))


#### filter: Keep rows matching criteria-------------------------------------

# Command structure (for all dplyr verbs):
#   first argument is a data frame
# return value is a data frame
# nothing is modified in place

# base R approach to view all flights on January 1
flights[flights$Month==1 & flights$DayofMonth==1, ]

# dplyr approach
# note: you can use comma or ampersand to represent AND condition
filter(flights, Month==1, DayofMonth==1)

# use pipe for OR condition
filter(flights, UniqueCarrier=="AA" | UniqueCarrier=="UA")

# you can also use %in% operator
filter(flights, UniqueCarrier %in% c("AA", "UA"))


####select: Pick columns by name ----------------------------------------

# base R approach to select DepTime, ArrTime, and FlightNum columns
flights[, c("DepTime", "ArrTime", "FlightNum")]

# dplyr approach
select(flights, DepTime, ArrTime, FlightNum)

# use colon to select multiple contiguous columns, and use `contains` to match columns by name
# note: `starts_with`, `ends_with`, and `matches` (for regular expressions) can also be used to match columns by name
select(flights, Year:DayofMonth, contains("Taxi"), contains("Delay"))


####"Chaining" or "Pipelining" -------------------------------------------

# nesting method to select UniqueCarrier and DepDelay columns and filter for delays over 60 minutes
filter(select(flights, UniqueCarrier, DepDelay), DepDelay > 60)

# chaining method
flights %>%
  select(UniqueCarrier, DepDelay) %>%
  filter(DepDelay > 60)

#Chaining increases readability significantly when there are many commands

# create two vectors and calculate Euclidian distance between them
x1 <- 1:5; x2 <- 2:6
sqrt(sum((x1-x2)^2))

# chaining method
(x1-x2)^2 %>% sum() %>% sqrt()