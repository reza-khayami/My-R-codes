rm(list = ls())
hist(islands, breaks = 16)
islands.z <- scale(islands)
?scale
hist(islands.z, breaks = 16)
mean(islands.z)
round(mean(islands.z), 2)
attr(islands.z, "scaled:center")
attr(islands.z, "scaled:scale")
island.z <- as.numeric(islands.z)
islands.ln <- log(islands)
hist(islands.ln)
i <- rank(islands.z)
i
?ifelse
continent <- ifelse(islands > 1000, 1, 0)
islands
