require("datasets")
?chickwts
?datasets
library(help = "datasets")
View(chickwts)
?data
plot(chickwts$feed, col = palette(heat (6)) )
rm(list = ls())
feeds <- table(chickwts$feed)
feeds
barplot(feeds, col = "#AEAA55")
?barplot
?order
barplot(feeds[order(feeds, decreasing = TRUE)])
barplot(feeds[order(feeds, decreasing = TRUE)], col = c("#EEFF55", "#11FFEE"))
?par
par(oma = c(1,3,5,1))
par(mar = c(4,5,2,1))
barplot(feeds[order(feeds)],
        horiz  = TRUE,
        las    = 1,
        col    = c("beige", "blanchedalmond", "bisque1", "bisque2","bisque3", "bisque4"),
        border = NA,
        main   = "Frequencies of Different Feeds\nin chickwts Dataset",
        xlab   ="Number of Chicks")
?par