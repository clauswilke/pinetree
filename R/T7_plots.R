library(tidyverse)
library(cowplot)
library(stringr)
library(ggrepel)

col_types = c(
  X1 = col_integer(),
  X2 = col_double(),
  X3 = col_character(),
  X4 = col_integer()
)

data <- read_delim("runs/full_T7_test_high_ribo.csv", ",", col_names = FALSE, col_types = col_types)
colnames(data) <- c("iteration", "time", "species", "count")

data$time <- as.numeric(data$time)
data$count <- as.numeric(data$count)

data <- filter(data, !str_detect(species, "_total$"), species != " rbs", !str_detect(species, "phi")) %>% filter(species != " rbs", !str_detect(species, "promoter"), !str_detect(species, "ecolipol"), !str_detect(species, "kinase"))

ggplot(data, aes(x = time, y = count, group = species, color = species)) + 
  geom_line() +
  geom_text_repel(
    data = filter(data, time == max(time)),
    aes(label = species),
    size = 4,
    nudge_x = 45,
    point.padding = unit(0.5, "lines")
  ) +
  coord_cartesian(xlim = c(0, max(data$time) + 150)) +
  theme(legend.position = "none")
