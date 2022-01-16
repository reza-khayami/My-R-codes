data("USJudgeRatings")
boxplot(USJudgeRatings, 
        horizontal = TRUE,
        xlab = "Lawyers' Ratings",
        las = 1,
        notch = FALSE,
        ylim = c(0,10),
        col = "#9FB6CD",
       boxwex = 0.5,
       whisklty = 1,
       staplelty = 0,
       outpch = 16,
       outcol = "slategray3",
       main = "Box plot\nReza"
       )
?boxplot
?width
par(mar=c(5, 4, 4, 5))
?par
