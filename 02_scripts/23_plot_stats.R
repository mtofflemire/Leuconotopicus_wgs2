############################################################
## Libraries
############################################################

library(tidyverse)

############################################################
## Paths
############################################################

BASE <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs"

PI_DIR <- paste0(BASE, "/03_results/05_genome_wide_variation_subspecies/pi")
TAJ_DIR <- paste0(BASE, "/03_results/05_genome_wide_variation_subspecies/tajimasD")
FST_DIR <- paste0(BASE, "/03_results/05_genome_wide_variation_subspecies/pairwise_fst")

############################################################
## Labels + order
############################################################

sub_labels <- c(
  "albolarvatus_CA_samps" = "ALB (CA)",
  "albolarvatus_WA_OR_samps" = "ALB (WA+OR)",
  "gravirostris_samps" = "GRA"
)

order_levels <- c(
  "albolarvatus_CA_samps",
  "albolarvatus_WA_OR_samps",
  "gravirostris_samps"
)

############################################################
## PI
############################################################

files <- list.files(PI_DIR, pattern = "windowed\\.pi$", full.names = TRUE)

pi_all <- map_df(files, function(f) {
  x <- read.table(f, header = TRUE)
  x$Population <- basename(f) |> sub("\\.windowed\\.pi$", "", x = _)
  x
})

pi_all <- pi_all %>% filter(!is.na(PI))
pi_all$Population <- factor(pi_all$Population, levels = order_levels)


p_pi <- ggplot(pi_all, aes(Population, PI)) +
  geom_boxplot(fill = "white", linewidth = 1, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, quantile(pi_all$PI, 0.99, na.rm = TRUE))) +
  scale_x_discrete(labels = sub_labels) +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "gray40", linewidth = 0.5),
    panel.grid.minor = element_line(color = "gray40", linewidth = 0.35)
  ) +
  labs(title = "Nucleotide Diversity", x = "Subspecies", y = expression(pi))

print(p_pi)

ggsave(
  filename = paste0(PI_DIR, "/nucleotide_diversity_pi_boxplot.pdf"),
  plot = p_pi,
  width = 3.6,
  height = 6
)

ggsave(
  filename = paste0(PI_DIR, "/nucleotide_diversity_pi_boxplot.png"),
  plot = p_pi,
  width = 6,
  height = 8,
  dpi = 300
)

############################################################
## Tajima D
############################################################

files <- list.files(TAJ_DIR, pattern = "Tajima\\.D$", full.names = TRUE)

taj_all <- map_df(files, function(f) {
  x <- read.table(f, header = TRUE)
  x$Population <- basename(f) |> sub("\\.Tajima\\.D$", "", x = _)
  x
})

taj_all <- taj_all %>% filter(!is.na(TajimaD))
taj_all$Population <- factor(taj_all$Population, levels = order_levels)

p_taj <- ggplot(taj_all, aes(Population, TajimaD)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40", linewidth = 0.7) +
  geom_boxplot(fill = "white", linewidth = 1, outlier.shape = NA) +
  coord_cartesian(ylim = quantile(taj_all$TajimaD, c(0.01, 0.99), na.rm = TRUE)) +
  scale_x_discrete(labels = sub_labels) +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "gray40", linewidth = 0.5),
    panel.grid.minor = element_line(color = "gray40", linewidth = 0.35)
  ) +
  labs(title = "Tajima's D", x = "Subspecies", y = "Tajima's D")

print(p_taj)

ggsave(
  filename = paste0(TAJ_DIR, "/tajimasD_boxplot.pdf"),
  plot = p_taj,
  width = 3.6,
  height = 6
)

ggsave(
  filename = paste0(TAJ_DIR, "/tajimasD_boxplot.png"),
  plot = p_taj,
  width = 6,
  height = 5,
  dpi = 300
)

############################################################
## FST
############################################################

files <- list.files(FST_DIR, pattern = "windowed\\.weir\\.fst$", full.names = TRUE)

fst_all <- map_df(files, function(f) {
  x <- read.table(f, header = TRUE)
  x$Comparison <- basename(f) |> sub("\\.windowed\\.weir\\.fst$", "", x = _)
  x
})

fst_all <- fst_all %>% filter(!is.na(MEAN_FST), MEAN_FST >= 0)

fst_all$Comparison <- fst_all$Comparison |>
  str_replace_all("albolarvatus_CA_samps", "ALB (CA)") |>
  str_replace_all("albolarvatus_WA_OR_samps", "ALB (WA+OR)") |>
  str_replace_all("gravirostris_samps", "GRA")

order_levels <- c(
  "ALB (CA)_vs_ALB (WA+OR)",
  "ALB (CA)_vs_GRA",
  "ALB (WA+OR)_vs_GRA"
)





fst_all$Facet_group <- fst_all$Comparison |>
  str_replace_all("_vs_", "-") |>
  str_replace_all("ALB \\(WA\\+OR\\)", "ALB(WA/OR)") |>
  str_replace_all("ALB \\(CA\\)", "ALB(CA)")

fst_summary <- fst_all %>%
  group_by(Facet_group, Comparison) %>%
  summarise(mean_fst = mean(MEAN_FST, na.rm = TRUE), .groups = "drop")

p_fst <- ggplot(fst_all, aes(Comparison, MEAN_FST)) +
  geom_violin(
    fill = "gray80",
    color = "gray20",
    linewidth = 0.8,
    trim = FALSE,
    scale = "width",
    width = 0.85
  ) +
  geom_boxplot(
    width = 0.16,
    fill = "white",
    color = "gray20",
    linewidth = 0.8,
    outlier.shape = NA
  ) +
  geom_point(
    data = fst_summary,
    aes(Comparison, mean_fst),
    inherit.aes = FALSE,
    size = 2.6,
    color = "gray20"
  ) +
  coord_flip(
    ylim = c(0, quantile(fst_all$MEAN_FST, 0.99, na.rm = TRUE))
  ) +
  facet_grid(Facet_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 18),
    plot.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.1),
    panel.grid.major.x = element_line(color = "gray45", linewidth = 0.7),
    panel.grid.minor.x = element_line(color = "gray70", linewidth = 0.45),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background.y = element_rect(fill = "black", color = "black", linewidth = 1.1),
    strip.text.y.left = element_text(color = "white", face = "bold", size = 15, angle = 90)
  ) +
  labs(x = NULL, y = expression(F[ST]))

print(p_fst)

