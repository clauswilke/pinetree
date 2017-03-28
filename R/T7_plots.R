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

data <- read_delim("runs/run_032517.out.txt", ",", col_names = FALSE, col_types = col_types)
colnames(data) <- c("iteration", "time", "species", "count")

data$time <- as.numeric(data$time)
data$count <- as.numeric(data$count)
data$species <- trimws(data$species)

classes <- read_csv("runs/classes.csv")

data <- filter(data, !str_detect(species, "_total$"), species != " rbs", !str_detect(species, "phi")) %>% filter(species != "rbs", !str_detect(species, "promoter"))

data <- left_join(data, classes)

plot <- ggplot(data, aes(x = time, y = count, group = species, color = factor(class, levels = c("ribo", 1, 2, 3)))) + 
  geom_line(size = 0.8) +
  geom_label_repel(
    data = filter(data, time == max(time)),
    aes(label = species, fill = factor(class)),
    size = 2,
    color = 'white',
    nudge_x = 45,
    point.padding = unit(0.5, "lines"),
    segment.color = 'grey50',
    force = 5
  ) +
  coord_cartesian(xlim = c(0, max(data$time) + 300)) +
  theme(legend.position = "none")

save_plot("runs/T7_032517.pdf", plot, base_height = 12, base_width = 14)
