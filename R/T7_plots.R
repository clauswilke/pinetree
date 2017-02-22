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

plot <- ggplot(data, aes(x = time, y = count, group = species, color = species)) + 
  geom_line() +
  geom_label_repel(
    data = filter(data, time == max(time)),
    aes(label = species, fill = species),
    size = 4,
    color = 'white',
    nudge_x = 45,
    point.padding = unit(0.5, "lines"),
    segment.color = 'grey50',
    force = 5
  ) +
  coord_cartesian(xlim = c(0, max(data$time) + 300)) +
  theme(legend.position = "none")

save_plot("runs/T7.pdf", plot, base_width = 12, base_height = 9)
