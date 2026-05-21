setwd("/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/01_data")


library(tidyverse)
library(readxl)
library(ggplot2)

meta <- read_excel("whwo_Meta.xlsx")
meta <- meta %>%
  filter(!is.na(SequenceID), !is.na(Ecoregion), !is.na(Subspecies))
meta$Ecoregion <- as.factor(meta$Ecoregion)
meta$Subspecies <- as.factor(meta$Subspecies)

custom_cols <- c(
  "#00A087",
  "#902F21",
  "#3C5488",
  "#4DBBD5",
  "#9DAAC4",
  "#E64B35",
  "#0F1522",
  "#F2A59A"
)

ecoregion_cols <- setNames(custom_cols[1:length(levels(meta$Ecoregion))],levels(meta$Ecoregion))
subspecies_shapes <- setNames(c(16, 17), levels(meta$Subspecies))


pca_dir <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/02_analysis/02_PCA'
eigenvec <- read_table(file.path(pca_dir, "whwo_pruned.eigenvec"), col_types = cols())
eigenval <- read_table(file.path(pca_dir, "whwo_pruned.eigenval"), col_names = "eigenval",col_types = cols())
percent_var <- eigenval$eigenval / sum(eigenval$eigenval) * 100


pca_data <- eigenvec %>%
  left_join(meta, by = c("IID" = "SequenceID"))


pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Ecoregion, shape = Subspecies)) +
  geom_point(size = 3, alpha = 0.95) +
  scale_color_manual(values = ecoregion_cols) +
  scale_shape_manual(values = subspecies_shapes) +
  labs(x = paste0("PC1 (", round(percent_var[1], 2), "%)"), y = paste0("PC2 (", round(percent_var[2], 2), "%)")) +
  theme_bw() +
  theme(legend.position = "right", legend.title = element_blank(), axis.title = element_text(size = 13), axis.text = element_text(size = 11), panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA))


pca_plot

ggsave(file.path(pca_dir, "leuconotopicus_PCA_PC1_PC2.png"), pca_plot, height = 6, width = 7, dpi = 600)
ggsave(file.path(pca_dir, "leuconotopicus_PCA_PC1_PC2.pdf"), pca_plot, height = 6, width = 7, useDingbats = FALSE)
