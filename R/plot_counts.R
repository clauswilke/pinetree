library(tidyverse)
library(cowplot)
library(ggrepel)

counts <- read_csv("runs/three_genes_out.csv", col_names = F)
colnames(counts) <- c("iteration", "time", "species", "count")
counts <- filter(counts, species %in% c("ribosome", "proteinX", "proteinY", "rnapol"))

plot1 <- ggplot(counts, aes(x = time, y = count, group = species, color = species)) + 
  geom_line(size = 0.8) +
  geom_label_repel(
    data = filter(counts, time == max(time)),
    aes(label = species, fill = species),
    size = 4,
    color = 'white',
    nudge_x = 5,
    point.padding = unit(0.2, "lines"),
    segment.color = 'grey80',
    force = 5
  ) +
  coord_cartesian(xlim = c(0, max(counts$time) + 10)) +
  theme(legend.position = "none")

save_plot("./runs/three_genes.pdf", plot1, base_aspect_ratio = 1.4)

counts <- read_csv("runs/dual_promoter_out.csv", col_names = F)
colnames(counts) <- c("iteration", "time", "species", "count")
counts <- filter(counts, species %in% c("ribosome", "proteinX", "proteinY", "rnapol"))

plot2 <- ggplot(counts, aes(x = time, y = count, group = species, color = species)) + 
  geom_line(size = 0.8) +
  geom_label_repel(
    data = filter(counts, time == max(time)),
    aes(label = species, fill = species),
    size = 4,
    color = 'white',
    nudge_x = 5,
    point.padding = unit(0.2, "lines"),
    segment.color = 'grey80',
    force = 5
  ) +
  coord_cartesian(xlim = c(0, max(counts$time) + 10)) +
  theme(legend.position = "none")

save_plot("./runs/dual_promoters.pdf", plot2, base_aspect_ratio = 1.4)


counts <- read_csv("runs/readthrough_out.csv", col_names = F)
colnames(counts) <- c("iteration", "time", "species", "count")
counts <- filter(counts, species %in% c("ribosome", "proteinX", "proteinY", "rnapol"))

plot3 <- ggplot(counts, aes(x = time, y = count, group = species, color = species)) + 
  geom_line(size = 0.8) +
  geom_label_repel(
    data = filter(counts, time == max(time)),
    aes(label = species, fill = species),
    size = 4,
    color = 'white',
    nudge_x = 5,
    point.padding = unit(0.2, "lines"),
    segment.color = 'grey80',
    force = 5
  ) +
  coord_cartesian(xlim = c(0, max(counts$time) + 10)) +
  theme(legend.position = "none")

save_plot("./runs/readthrough.pdf", plot3, base_aspect_ratio = 1.4)

counts <- read_csv("runs/dual_polymerases_out.csv", col_names = F)
colnames(counts) <- c("iteration", "time", "species", "count")
counts <- filter(counts, species %in% c("ecoli", "rnapol", "ribosome", "proteinX", "proteinY"))

plot4 <- ggplot(counts, aes(x = time, y = count, group = species, color = species)) + 
  geom_line(size = 0.8) +
  geom_label_repel(
    data = filter(counts, time == max(time)),
    aes(label = species, fill = species),
    size = 4,
    color = 'white',
    nudge_x = 5,
    point.padding = unit(0.2, "lines"),
    segment.color = 'grey80',
    force = 5
  ) +
  coord_cartesian(xlim = c(0, max(counts$time) + 10)) +
  theme(legend.position = "none")

save_plot("./runs/dual_polymerases.pdf", plot4, base_aspect_ratio = 1.4)


counts <- read_tsv("/Users/ben/Box Sync/scratch/pysinthe_tests/recoded_test_counts.tsv", col_names = F)
colnames(counts) <- c("time", "species", "count")
counts <- filter(counts, species %in% c("ribosome", "proteinX", "proteinY", "rnapol"))

plot5 <- ggplot(counts, aes(x = time, y = count, group = species, color = species)) + 
  geom_line(size = 0.8) +
  geom_label_repel(
    data = filter(counts, time == max(time)),
    aes(label = species, fill = species),
    size = 4,
    color = 'white',
    nudge_x = 5,
    point.padding = unit(0.2, "lines"),
    segment.color = 'grey80',
    force = 5
  ) +
  coord_cartesian(xlim = c(0, max(counts$time) + 10), ylim=c(0, 310)) +
  theme(legend.position = "none")

save_plot("~/Desktop/recoded.pdf", plot5, base_aspect_ratio = 1.4)
