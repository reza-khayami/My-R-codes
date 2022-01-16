install.packages(c("reshape2", "ggplot2"))
a <- dcast(french_fries, treatment ~ subject,
           value.var = "potato", fun.aggregate = function(x){mean(x, na.rm = T)})
