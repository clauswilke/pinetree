library(cowplot)
library(tidyverse)
library(stringr)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

three_genes_data <- read_tsv("three_genes_counts.tsv") %>% 
  filter(!str_detect(species, "^__"), species != "p1") %>%
  gather(type, value, -time, -species)

recoded_data <- read_tsv("three_genes_recoded_counts.tsv") %>% 
  filter(!str_detect(species, "^__"), species != "p1") %>%
  gather(type, value, -time, -species)

recoded <- ggplot(recoded_data %>% filter(type != "ribo_density"), 
                   aes(x = time, y = value, group = interaction(species, type), color = species)) + 
  geom_line(aes(linetype = type)) +
  ylab("count") + 
  scale_color_manual("species\n(unbound)", values = cbbPalette[c(2,3,4,8)])

protein_recoded <- recoded_data %>% 
  filter(species == "proteinY", type == "ribo_density") %>% 
  mutate(condition = "recoded")
protein_wt <- three_genes_data %>% 
  filter(species == "proteinY", type == "ribo_density") %>% 
  mutate(condition = "wildtype")

protein_y <- bind_rows(protein_recoded, protein_wt)

protein_y_plot <- ggplot(protein_y, aes(x = time, y = value, color = condition)) + 
  geom_line() +
  ylab("ribosomes per transcript") +
  scale_color_manual(values = cbbPalette[c(6,7)])

recoded_plot_all <- plot_grid(recoded, protein_y_plot)
save_plot("plots.pdf", recoded_plot_all, base_width = 9, base_height = 3.3)

degrade_data <- read_tsv("three_genes_rnase.tsv") %>% 
  filter(!str_detect(species, "^__"), species != "p1", species != "p2") %>%
  gather(type, value, -time, -species) %>% 
  filter(type == "transcript", species != "rnapol")

degrade_plot <- ggplot(degrade_data, aes(x=time, y=value, group=species, color=species)) +
  geom_line() + xlab("time (s)") + ylab("transcript count")

save_plot("degrade_plot.png", degrade_plot, base_aspect_ratio = 1.3)


