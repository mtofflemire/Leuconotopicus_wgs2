# load libraries
library(tidyverse)
library(gtools)
library(patchwork)

# file paths
base_dir <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/02_analysis/03_ADMIXTURE/seed43'
prefix <- "whwp_SNPdata1"

famfile <- file.path(base_dir, paste0(prefix, ".fam"))
metadata_file <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/01_data/whwo_Meta.xlsx'

k_values <- 2:10
stacked_k_values <- 3:10

ecoregion_order <- c(
  "Eastern Cascades forests", "Blue Mountains forests", 
  "Central-Southern Cascades Forests", 
  "Klamath-Siskiyou forests", "Sierra Nevada forests",
  "California montane chaparral and woodlands"
)

ecoregion_abbrev <- c(
  "ECF", "BM", "CSCF", "KSF", "SNF", "CMCW"
)

ecoregion_labels <- paste0(ecoregion_abbrev, " [", seq_along(ecoregion_abbrev), "]")

cluster_palette <- c(
  "ancestral1" = "salmon",
  "ancestral2" = "cyan3",
  "ancestral3" = "green",
  "ancestral4" = "purple",
  "ancestral5" = "orange",
  "ancestral6" = "blue",
  "ancestral7" = "magenta",
  "ancestral8" = "gold",
  "ancestral9" = "gray30",
  "ancestral10" = "brown"
)

# load shared data
fam <- read.table(famfile, header = FALSE)
metadata <- read_excel(metadata_file)
metadata
sample_metadata <-
  tibble(SequenceID = fam$V2) %>%
  left_join(metadata %>% select(SequenceID, Ecoregion), by = "SequenceID") %>%
  filter(!is.na(Ecoregion)) %>%
  mutate(
    Ecoregion = factor(Ecoregion, levels = ecoregion_order),
    EcoAbbrev = factor(
      ecoregion_labels[match(Ecoregion, ecoregion_order)],
      levels = ecoregion_labels
    )
  )

# establish K = 2 sample order and cluster identity
q2_raw <-
  read.table(file.path(base_dir, paste0(prefix, ".2.Q")), header = FALSE) %>%
  as_tibble() %>%
  setNames(c("ancestral1", "ancestral2")) %>%
  mutate(SequenceID = fam$V2) %>%
  left_join(sample_metadata, by = "SequenceID") %>%
  filter(!is.na(EcoAbbrev)) %>%
  arrange(EcoAbbrev, desc(ancestral1))

sample_order <- q2_raw$SequenceID

q2_matrix <-
  q2_raw %>%
  select(SequenceID, ancestral1, ancestral2) %>%
  arrange(match(SequenceID, sample_order))

# align K clusters to K = 2 clusters
align_to_k2 <- function(K) {
  
  qfile <- file.path(base_dir, paste0(prefix, ".", K, ".Q"))
  original_names <- paste0("raw", 1:K)
  
  qk_raw <-
    read.table(qfile, header = FALSE) %>%
    as_tibble() %>%
    setNames(original_names) %>%
    mutate(SequenceID = fam$V2) %>%
    filter(SequenceID %in% sample_order) %>%
    arrange(match(SequenceID, sample_order))
  
  if (K == 2) {
    return(
      qk_raw %>%
        transmute(SequenceID, ancestral1 = raw1, ancestral2 = raw2)
    )
  }
  
  q2_values <- as.matrix(q2_matrix[, c("ancestral1", "ancestral2")])
  qk_values <- as.matrix(qk_raw[, original_names])
  
  perms <- permutations(n = K, r = 2, v = seq_len(K))
  
  scores <- apply(perms, 1, function(p) {
    cor(q2_values[, 1], qk_values[, p[1]]) +
      cor(q2_values[, 2], qk_values[, p[2]])
  })
  
  best <- perms[which.max(scores), ]
  extra <- setdiff(seq_len(K), best)
  
  ordered_cols <- c(best, extra)
  aligned_names <- paste0("ancestral", seq_len(K))
  
  qk_raw %>%
    select(SequenceID, all_of(original_names[ordered_cols])) %>%
    setNames(c("SequenceID", aligned_names))
}

# prepare plotting data
make_plot_data <- function(K) {
  
  cluster_names <- paste0("ancestral", 1:K)
  
  align_to_k2(K) %>%
    left_join(sample_metadata, by = "SequenceID") %>%
    mutate(
      SequenceID = factor(SequenceID, levels = sample_order),
      K_label = factor(paste0("K = ", K), levels = paste0("K = ", k_values))
    ) %>%
    arrange(SequenceID) %>%
    pivot_longer(
      cols = all_of(cluster_names),
      names_to = "cluster",
      values_to = "ancestry"
    ) %>%
    mutate(cluster = factor(cluster, levels = rev(cluster_names)))
}

# individual plot function
plot_admixture_individual <- function(K) {
  
  cluster_names <- paste0("ancestral", 1:K)
  plot_data <- make_plot_data(K)
  
  ggplot(plot_data, aes(SequenceID, ancestry, fill = cluster)) +
    geom_col(color = "gray30", linewidth = 0.08, width = 1) +
    facet_grid(~ EcoAbbrev, switch = "x", scales = "free_x", space = "free_x") +
    labs(title = paste0("K = ", K), y = "Ancestry", x = NULL) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 0.7)) +
    scale_fill_manual(values = cluster_palette[cluster_names], guide = "none") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, size = 10, face = "bold"),
      plot.title.position = "plot",
      panel.spacing.x = unit(0.08, "lines"),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.placement = "outside",
      strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 7),
      axis.title.y = element_text(size = 8),
      axis.text.y = element_text(size = 7)
    )
}

# stacked plot data
stacked_plot_data <-
  map_dfr(stacked_k_values, make_plot_data) %>%
  mutate(K_label = factor(K_label, levels = paste0("K = ", stacked_k_values)))

# stacked plot
stacked_plot <-
  ggplot(stacked_plot_data, aes(SequenceID, ancestry, fill = cluster)) +
  geom_col(color = "gray30", linewidth = 0.05, width = 1) +
  facet_grid(K_label ~ EcoAbbrev, switch = "x", scales = "free_x", space = "free_x") +
  labs(y = "Ancestry", x = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 0.7)) +
  scale_fill_manual(values = cluster_palette, guide = "none") +
  theme_minimal() +
  theme(
    panel.spacing.x = unit(0.08, "lines"),
    panel.spacing.y = unit(1, "lines"),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    strip.text.y = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 8)
  )

# save individual plots
for (K in k_values) {
  
  p <- plot_admixture_individual(K)
  
  print(p)
  
  ggsave(
    file.path(base_dir, paste0(prefix, "_ADMIXTURE_K", K, "_by_ecoregion_K2_aligned.pdf")),
    p,
    width = 10,
    height = 3,
    useDingbats = FALSE
  )
  
  ggsave(
    file.path(base_dir, paste0(prefix, "_ADMIXTURE_K", K, "_by_ecoregion_K2_aligned.png")),
    p,
    width = 10,
    height = 3,
    dpi = 600
  )
}

# save stacked K4-K10 plot
stacked_plot

ggsave(
  file.path(base_dir, paste0(prefix, "_ADMIXTURE_K4-K10_stacked_K2_aligned.pdf")),
  stacked_plot,
  width = 8,
  height = 10,
  useDingbats = FALSE
)

ggsave(
  file.path(base_dir, paste0(prefix, "_ADMIXTURE_K4-K10_stacked_K2_aligned.png")),
  stacked_plot,
  width = 8,
  height = 10,
  dpi = 600
)




