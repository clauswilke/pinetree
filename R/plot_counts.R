library(tidyverse)
library(cowplot)

counts <- read_csv("runs/three_genes_out.csv", col_names = F)
colnames(counts) <- c("iteration", "time", "species", "count")
counts <- filter(counts, species %in% c("ribosome", "proteinX", "capsid"))

ggplot(counts, aes(x = time, y = count, group = species, color = species)) + geom_line()

counts <- read_csv("runs/dual_promoter_out.csv", col_names = F)
colnames(counts) <- c("iteration", "time", "species", "count")
counts <- filter(counts, species %in% c("ribosome", "proteinX", "capsid"))

ggplot(counts, aes(x = time, y = count, group = species, color = species)) + geom_line()

counts <- read_csv("runs/dual_polymerases_out.csv", col_names = F)
colnames(counts) <- c("iteration", "time", "species", "count")
counts <- filter(counts, species %in% c("ecoli", "ribosome", "proteinX", "capsid"))

ggplot(counts, aes(x = time, y = count, group = species, color = species)) + geom_line()




