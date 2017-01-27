library(tidyverse)
library(cowplot)

times <- read_csv("times.csv", col_names = F)
colnames(times) <- c("name", "time")

ggplot(times, aes(x = time)) + geom_histogram()




