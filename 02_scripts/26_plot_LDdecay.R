library(tidyverse)


LD_DIR <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/03_results/06_ld_decay"

files <- c(
  "ALB (CA)" = file.path(LD_DIR, "leuconotopicus_subspecies_lddecay.ALB_CA"),
  "ALB (WA+OR)" = file.path(LD_DIR, "leuconotopicus_subspecies_lddecay.ALB_WA_OR"),
  "GRA" = file.path(LD_DIR, "leuconotopicus_subspecies_lddecay.GRA")
)

ld_all <- imap_dfr(files, function(f, pop) {
  read.table(f, header = TRUE) %>%
    mutate(Population = pop)
})

names(ld_all)




library(tidyverse)

LD_DIR <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/03_results/06_ld_decay"

files <- c(
  "ALB (CA)" = file.path(LD_DIR, "leuconotopicus_subspecies_lddecay.ALB_CA"),
  "ALB (WA+OR)" = file.path(LD_DIR, "leuconotopicus_subspecies_lddecay.ALB_WA_OR"),
  "GRA" = file.path(LD_DIR, "leuconotopicus_subspecies_lddecay.GRA")
)

ld_all <- imap_dfr(files, function(f, pop) {
  read.table(f, header = FALSE, fill = TRUE) %>%
    rename(
      Distance_bp = V1,
      Mean_r2 = V2
    ) %>%
    mutate(
      Distance_kb = Distance_bp / 1000,
      Population = pop
    )
})

p_ld <- ggplot(ld_all, aes(x = Distance_kb, y = Mean_r2, color = Population)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c(
    "ALB (CA)" = "#0072B2",
    "ALB (WA+OR)" = "#D55E00",
    "GRA" = "#009E73"
  )) +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 18, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    panel.grid.major = element_line(color = "gray85", linewidth = 0.5),
    panel.grid.minor = element_line(color = "gray92", linewidth = 0.35),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "LD Decay",
    x = "Distance (kb)",
    y = expression(r^2)
  )

print(p_ld)

ggsave(file.path(LD_DIR, "leuconotopicus_subspecies_lddecay_ggplot.pdf"), p_ld, width = 7, height = 5)
ggsave(file.path(LD_DIR, "leuconotopicus_subspecies_lddecay_ggplot.png"), p_ld, width = 7, height = 5, dpi = 300)