library(StatProg)
library(ggplot2)


galaxies <- as.data.frame(galaxies)
names(galaxies) <- "km"

ggplot(galaxies, aes(x = km)) +
    geom_density()

