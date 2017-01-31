library(tidyverse)
library(cowplot)

counts <- read_csv("test2.csv", col_names = F)
colnames(counts) <- c("iteration", "time", "species", "count")

ggplot(counts, aes(x = time, y = count, group = species, color = species)) + geom_line(size = 1)




