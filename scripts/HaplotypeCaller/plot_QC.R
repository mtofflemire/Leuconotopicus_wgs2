library(readr)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)
library(gridExtra)
library(data.table)

qc_dir <- "/usr/scratch2/userdata2/mtofflemire/projects/leuconotopicus/reseq_WD/04_vcf"
out_dir <- file.path(qc_dir, "site_qc_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("[", Sys.time(), "] Starting plot_QC.R")
message("[", Sys.time(), "] QC directory: ", qc_dir)
message("[", Sys.time(), "] Output directory: ", out_dir)

message("[", Sys.time(), "] Reading invariant QC files")
inv_imiss <- read_delim(file.path(qc_dir, "leuconotopicus.raw.invariant.imiss"), delim = "\t", show_col_types = FALSE)
inv_idepth <- read_delim(file.path(qc_dir, "leuconotopicus.raw.invariant.idepth"), delim = "\t", show_col_types = FALSE)
inv_lmiss <- read_delim(file.path(qc_dir, "leuconotopicus.raw.invariant.lmiss"), delim = "\t", show_col_types = FALSE)
inv_ldepth <- read_delim(file.path(qc_dir, "leuconotopicus.raw.invariant.ldepth.mean"), delim = "\t", show_col_types = FALSE)

message("[", Sys.time(), "] Reading variant QC files")
var_imiss <- read_delim(file.path(qc_dir, "leuconotopicus.raw.variants.imiss"), delim = "\t", show_col_types = FALSE)
var_idepth <- read_delim(file.path(qc_dir, "leuconotopicus.raw.variants.idepth"), delim = "\t", show_col_types = FALSE)
var_lmiss <- read_delim(file.path(qc_dir, "leuconotopicus.raw.variants.lmiss"), delim = "\t", show_col_types = FALSE)
var_ldepth <- read_delim(file.path(qc_dir, "leuconotopicus.raw.variants.ldepth.mean"), delim = "\t", show_col_types = FALSE)
var_lqual <- read_delim(file.path(qc_dir, "leuconotopicus.raw.variants.lqual"), delim = "\t", show_col_types = FALSE)
var_het <- read_delim(file.path(qc_dir, "leuconotopicus.raw.variants.het"), delim = "\t", show_col_types = FALSE)

message("[", Sys.time(), "] Sampling 100000 annotation rows before loading into R")
ann_file <- file.path(qc_dir, "leuconotopicus.raw.variants.annotations.table")
sample_cmd <- sprintf(
  "{ head -n 1 %s; tail -n +2 %s | shuf -n 100000; }",
  shQuote(ann_file),
  shQuote(ann_file)
)

ann <- fread(
  cmd = sample_cmd,
  sep = "\t",
  na.strings = c("NA", "")
) %>%
  as_tibble() %>%
  select(QUAL, QD, DP, MQ, MQRankSum, FS, ReadPosRankSum, SOR)

message("[", Sys.time(), "] Annotation rows kept: ", nrow(ann))

message("[", Sys.time(), "] Cleaning numeric columns")
inv_imiss <- inv_imiss %>%
  mutate(F_MISS = as.numeric(F_MISS)) %>%
  filter(is.finite(F_MISS))

inv_idepth <- inv_idepth %>%
  mutate(MEAN_DEPTH = as.numeric(MEAN_DEPTH)) %>%
  filter(is.finite(MEAN_DEPTH))

inv_lmiss <- inv_lmiss %>%
  mutate(F_MISS = as.numeric(F_MISS)) %>%
  filter(is.finite(F_MISS))

inv_ldepth <- inv_ldepth %>%
  mutate(MEAN_DEPTH = as.numeric(MEAN_DEPTH)) %>%
  filter(is.finite(MEAN_DEPTH))

var_imiss <- var_imiss %>%
  mutate(F_MISS = as.numeric(F_MISS)) %>%
  filter(is.finite(F_MISS))

var_idepth <- var_idepth %>%
  mutate(MEAN_DEPTH = as.numeric(MEAN_DEPTH)) %>%
  filter(is.finite(MEAN_DEPTH))

var_lmiss <- var_lmiss %>%
  mutate(F_MISS = as.numeric(F_MISS)) %>%
  filter(is.finite(F_MISS))

var_ldepth <- var_ldepth %>%
  mutate(MEAN_DEPTH = as.numeric(MEAN_DEPTH)) %>%
  filter(is.finite(MEAN_DEPTH))

var_lqual <- var_lqual %>%
  mutate(QUAL = as.numeric(QUAL)) %>%
  filter(is.finite(QUAL))

var_het <- var_het %>%
  mutate(F = as.numeric(F)) %>%
  filter(is.finite(F))

ann <- ann %>%
  mutate(
    QUAL = as.numeric(QUAL),
    QD = as.numeric(QD),
    DP = as.numeric(DP),
    MQ = as.numeric(MQ),
    MQRankSum = as.numeric(MQRankSum),
    FS = as.numeric(FS),
    ReadPosRankSum = as.numeric(ReadPosRankSum),
    SOR = as.numeric(SOR)
  )

message("[", Sys.time(), "] Writing missingness/depth summary tables")
site_summary <- tibble(
  metric = c("invariant_site_missingness", "invariant_site_depth", "variant_site_missingness", "variant_site_depth"),
  n = c(nrow(inv_lmiss), nrow(inv_ldepth), nrow(var_lmiss), nrow(var_ldepth)),
  mean = c(
    mean(inv_lmiss$F_MISS, na.rm = TRUE),
    mean(inv_ldepth$MEAN_DEPTH, na.rm = TRUE),
    mean(var_lmiss$F_MISS, na.rm = TRUE),
    mean(var_ldepth$MEAN_DEPTH, na.rm = TRUE)
  ),
  median = c(
    median(inv_lmiss$F_MISS, na.rm = TRUE),
    median(inv_ldepth$MEAN_DEPTH, na.rm = TRUE),
    median(var_lmiss$F_MISS, na.rm = TRUE),
    median(var_ldepth$MEAN_DEPTH, na.rm = TRUE)
  ),
  q01 = c(
    quantile(inv_lmiss$F_MISS, 0.01, na.rm = TRUE),
    quantile(inv_ldepth$MEAN_DEPTH, 0.01, na.rm = TRUE),
    quantile(var_lmiss$F_MISS, 0.01, na.rm = TRUE),
    quantile(var_ldepth$MEAN_DEPTH, 0.01, na.rm = TRUE)
  ),
  q05 = c(
    quantile(inv_lmiss$F_MISS, 0.05, na.rm = TRUE),
    quantile(inv_ldepth$MEAN_DEPTH, 0.05, na.rm = TRUE),
    quantile(var_lmiss$F_MISS, 0.05, na.rm = TRUE),
    quantile(var_ldepth$MEAN_DEPTH, 0.05, na.rm = TRUE)
  ),
  q25 = c(
    quantile(inv_lmiss$F_MISS, 0.25, na.rm = TRUE),
    quantile(inv_ldepth$MEAN_DEPTH, 0.25, na.rm = TRUE),
    quantile(var_lmiss$F_MISS, 0.25, na.rm = TRUE),
    quantile(var_ldepth$MEAN_DEPTH, 0.25, na.rm = TRUE)
  ),
  q75 = c(
    quantile(inv_lmiss$F_MISS, 0.75, na.rm = TRUE),
    quantile(inv_ldepth$MEAN_DEPTH, 0.75, na.rm = TRUE),
    quantile(var_lmiss$F_MISS, 0.75, na.rm = TRUE),
    quantile(var_ldepth$MEAN_DEPTH, 0.75, na.rm = TRUE)
  ),
  q95 = c(
    quantile(inv_lmiss$F_MISS, 0.95, na.rm = TRUE),
    quantile(inv_ldepth$MEAN_DEPTH, 0.95, na.rm = TRUE),
    quantile(var_lmiss$F_MISS, 0.95, na.rm = TRUE),
    quantile(var_ldepth$MEAN_DEPTH, 0.95, na.rm = TRUE)
  ),
  q99 = c(
    quantile(inv_lmiss$F_MISS, 0.99, na.rm = TRUE),
    quantile(inv_ldepth$MEAN_DEPTH, 0.99, na.rm = TRUE),
    quantile(var_lmiss$F_MISS, 0.99, na.rm = TRUE),
    quantile(var_ldepth$MEAN_DEPTH, 0.99, na.rm = TRUE)
  )
)

sample_summary <- tibble(
  metric = c("invariant_sample_missingness", "invariant_sample_depth", "variant_sample_missingness", "variant_sample_depth"),
  n = c(nrow(inv_imiss), nrow(inv_idepth), nrow(var_imiss), nrow(var_idepth)),
  mean = c(
    mean(inv_imiss$F_MISS, na.rm = TRUE),
    mean(inv_idepth$MEAN_DEPTH, na.rm = TRUE),
    mean(var_imiss$F_MISS, na.rm = TRUE),
    mean(var_idepth$MEAN_DEPTH, na.rm = TRUE)
  ),
  median = c(
    median(inv_imiss$F_MISS, na.rm = TRUE),
    median(inv_idepth$MEAN_DEPTH, na.rm = TRUE),
    median(var_imiss$F_MISS, na.rm = TRUE),
    median(var_idepth$MEAN_DEPTH, na.rm = TRUE)
  ),
  q01 = c(
    quantile(inv_imiss$F_MISS, 0.01, na.rm = TRUE),
    quantile(inv_idepth$MEAN_DEPTH, 0.01, na.rm = TRUE),
    quantile(var_imiss$F_MISS, 0.01, na.rm = TRUE),
    quantile(var_idepth$MEAN_DEPTH, 0.01, na.rm = TRUE)
  ),
  q05 = c(
    quantile(inv_imiss$F_MISS, 0.05, na.rm = TRUE),
    quantile(inv_idepth$MEAN_DEPTH, 0.05, na.rm = TRUE),
    quantile(var_imiss$F_MISS, 0.05, na.rm = TRUE),
    quantile(var_idepth$MEAN_DEPTH, 0.05, na.rm = TRUE)
  ),
  q25 = c(
    quantile(inv_imiss$F_MISS, 0.25, na.rm = TRUE),
    quantile(inv_idepth$MEAN_DEPTH, 0.25, na.rm = TRUE),
    quantile(var_imiss$F_MISS, 0.25, na.rm = TRUE),
    quantile(var_idepth$MEAN_DEPTH, 0.25, na.rm = TRUE)
  ),
  q75 = c(
    quantile(inv_imiss$F_MISS, 0.75, na.rm = TRUE),
    quantile(inv_idepth$MEAN_DEPTH, 0.75, na.rm = TRUE),
    quantile(var_imiss$F_MISS, 0.75, na.rm = TRUE),
    quantile(var_idepth$MEAN_DEPTH, 0.75, na.rm = TRUE)
  ),
  q95 = c(
    quantile(inv_imiss$F_MISS, 0.95, na.rm = TRUE),
    quantile(inv_idepth$MEAN_DEPTH, 0.95, na.rm = TRUE),
    quantile(var_imiss$F_MISS, 0.95, na.rm = TRUE),
    quantile(var_idepth$MEAN_DEPTH, 0.95, na.rm = TRUE)
  ),
  q99 = c(
    quantile(inv_imiss$F_MISS, 0.99, na.rm = TRUE),
    quantile(inv_idepth$MEAN_DEPTH, 0.99, na.rm = TRUE),
    quantile(var_imiss$F_MISS, 0.99, na.rm = TRUE),
    quantile(var_idepth$MEAN_DEPTH, 0.99, na.rm = TRUE)
  )
)

write.table(
  site_summary,
  file = file.path(out_dir, "site_missingness_depth_summary.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sample_summary,
  file = file.path(out_dir, "sample_missingness_depth_summary.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("[", Sys.time(), "] Calculating annotation quantiles")
quantile_table <- tibble(
  statistic = c("QUAL", "DP", "QD", "FS", "MQ", "SOR", "MQRankSum", "ReadPosRankSum"),
  q01 = c(
    quantile(ann$QUAL, 0.01, na.rm = TRUE),
    quantile(ann$DP, 0.01, na.rm = TRUE),
    quantile(ann$QD, 0.01, na.rm = TRUE),
    quantile(ann$FS, 0.01, na.rm = TRUE),
    quantile(ann$MQ, 0.01, na.rm = TRUE),
    quantile(ann$SOR, 0.01, na.rm = TRUE),
    quantile(ann$MQRankSum, 0.01, na.rm = TRUE),
    quantile(ann$ReadPosRankSum, 0.01, na.rm = TRUE)
  ),
  q05 = c(
    quantile(ann$QUAL, 0.05, na.rm = TRUE),
    quantile(ann$DP, 0.05, na.rm = TRUE),
    quantile(ann$QD, 0.05, na.rm = TRUE),
    quantile(ann$FS, 0.05, na.rm = TRUE),
    quantile(ann$MQ, 0.05, na.rm = TRUE),
    quantile(ann$SOR, 0.05, na.rm = TRUE),
    quantile(ann$MQRankSum, 0.05, na.rm = TRUE),
    quantile(ann$ReadPosRankSum, 0.05, na.rm = TRUE)
  ),
  q95 = c(
    quantile(ann$QUAL, 0.95, na.rm = TRUE),
    quantile(ann$DP, 0.95, na.rm = TRUE),
    quantile(ann$QD, 0.95, na.rm = TRUE),
    quantile(ann$FS, 0.95, na.rm = TRUE),
    quantile(ann$MQ, 0.95, na.rm = TRUE),
    quantile(ann$SOR, 0.95, na.rm = TRUE),
    quantile(ann$MQRankSum, 0.95, na.rm = TRUE),
    quantile(ann$ReadPosRankSum, 0.95, na.rm = TRUE)
  ),
  q99 = c(
    quantile(ann$QUAL, 0.99, na.rm = TRUE),
    quantile(ann$DP, 0.99, na.rm = TRUE),
    quantile(ann$QD, 0.99, na.rm = TRUE),
    quantile(ann$FS, 0.99, na.rm = TRUE),
    quantile(ann$MQ, 0.99, na.rm = TRUE),
    quantile(ann$SOR, 0.99, na.rm = TRUE),
    quantile(ann$MQRankSum, 0.99, na.rm = TRUE),
    quantile(ann$ReadPosRankSum, 0.99, na.rm = TRUE)
  )
)

write.table(
  quantile_table,
  file = file.path(out_dir, "variant_hardfilter_quantiles_table.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

q_qual <- quantile(ann$QUAL, c(0.01, 0.99), na.rm = TRUE)
q_dp <- quantile(ann$DP, c(0.01, 0.99), na.rm = TRUE)
q_qd <- quantile(ann$QD, c(0.01, 0.99), na.rm = TRUE)
q_fs <- quantile(ann$FS, c(0.01, 0.99), na.rm = TRUE)
q_mq <- quantile(ann$MQ, c(0.01, 0.99), na.rm = TRUE)
q_sor <- quantile(ann$SOR, c(0.01, 0.99), na.rm = TRUE)
q_mqrank <- quantile(ann$MQRankSum, c(0.01, 0.99), na.rm = TRUE)
q_rprs <- quantile(ann$ReadPosRankSum, c(0.01, 0.99), na.rm = TRUE)

dens_layer <- function(...) {
  geom_density(
    fill = "dodgerblue1",
    colour = "black",
    alpha = 0.3,
    adjust = 2.5,
    n = 4096,
    linewidth = 0.5,
    ...
  )
}

message("[", Sys.time(), "] Building plots")
p1 <- ggplot(inv_imiss, aes(x = F_MISS, y = reorder(INDV, F_MISS))) +
  geom_point(color = "firebrick2", size = 1.8, alpha = 0.85) +
  theme_light() +
  labs(x = "Fraction missing", y = "Sample", title = "Invariant sample missingness")

p2 <- ggplot(inv_idepth, aes(x = MEAN_DEPTH, y = reorder(INDV, MEAN_DEPTH))) +
  geom_point(color = "steelblue", size = 1.8, alpha = 0.85) +
  theme_light() +
  labs(x = "Mean depth", y = "Sample", title = "Invariant sample mean depth")

p3 <- ggplot(inv_lmiss, aes(x = F_MISS)) +
  dens_layer() +
  theme_light() +
  labs(x = "Fraction missing", y = "Density", title = "Invariant site missingness")

p4 <- ggplot(inv_ldepth, aes(x = MEAN_DEPTH)) +
  dens_layer() +
  coord_cartesian(xlim = c(0, 30)) +
  theme_light() +
  labs(x = "Mean site depth", y = "Density", title = "Invariant site mean depth")

p5 <- ggplot(var_imiss, aes(x = F_MISS, y = reorder(INDV, F_MISS))) +
  geom_point(color = "firebrick2", size = 1.8, alpha = 0.85) +
  theme_light() +
  labs(x = "Fraction missing", y = "Sample", title = "Variant sample missingness")

p6 <- ggplot(var_idepth, aes(x = MEAN_DEPTH, y = reorder(INDV, MEAN_DEPTH))) +
  geom_point(color = "steelblue", size = 1.8, alpha = 0.85) +
  theme_light() +
  labs(x = "Mean depth", y = "Sample", title = "Variant sample mean depth")

p7 <- ggplot(var_lmiss, aes(x = F_MISS)) +
  dens_layer() +
  theme_light() +
  labs(x = "Fraction missing", y = "Density", title = "Variant site missingness")

p8 <- ggplot(var_ldepth, aes(x = MEAN_DEPTH)) +
  dens_layer() +
  coord_cartesian(xlim = c(0, 30)) +
  theme_light() +
  labs(x = "Mean site depth", y = "Density", title = "Variant site mean depth")

p9 <- ggplot(var_lqual, aes(x = QUAL)) +
  dens_layer() +
  coord_cartesian(xlim = c(0, 1000)) +
  theme_light() +
  labs(x = "Site QUAL", y = "Density", title = "Variant site quality")

p10 <- ggplot(var_het, aes(x = F)) +
  dens_layer() +
  theme_light() +
  labs(x = "Inbreeding coefficient (F)", y = "Density", title = "Variant sample heterozygosity")

p11 <- ggplot(filter(ann, is.finite(QUAL), QUAL > 0), aes(x = QUAL)) +
  dens_layer() +
  geom_vline(xintercept = q_qual, colour = "blue", linetype = "dashed") +
  coord_cartesian(xlim = c(0, q_qual[2])) +
  theme_light() +
  labs(x = "QUAL", y = "Density", title = "QUAL")

p12 <- ggplot(filter(ann, is.finite(DP), DP > 0), aes(x = DP)) +
  dens_layer() +
  geom_vline(xintercept = q_dp, colour = "blue", linetype = "dashed") +
  coord_cartesian(xlim = c(0, q_dp[2])) +
  theme_light() +
  labs(x = "DP", y = "Density", title = "DP")

p13 <- ggplot(filter(ann, is.finite(QD)), aes(x = QD)) +
  dens_layer() +
  geom_vline(xintercept = 2, colour = "red") +
  geom_vline(xintercept = q_qd, colour = "blue", linetype = "dashed") +
  theme_light() +
  labs(x = "QD", y = "Density", title = "QD")

p14 <- ggplot(filter(ann, is.finite(FS), FS >= 0), aes(x = FS)) +
  dens_layer() +
  geom_vline(xintercept = 60, colour = "red") +
  geom_vline(xintercept = q_fs, colour = "blue", linetype = "dashed") +
  coord_cartesian(xlim = c(0, min(q_fs[2], 100))) +
  theme_light() +
  labs(x = "FS", y = "Density", title = "FS")

p15 <- ggplot(filter(ann, is.finite(MQ)), aes(x = MQ)) +
  dens_layer() +
  geom_vline(xintercept = 40, colour = "red") +
  geom_vline(xintercept = q_mq, colour = "blue", linetype = "dashed") +
  theme_light() +
  labs(x = "MQ", y = "Density", title = "MQ")

p16 <- ggplot(filter(ann, is.finite(SOR), SOR >= 0), aes(x = SOR)) +
  dens_layer() +
  geom_vline(xintercept = 3, colour = "red") +
  geom_vline(xintercept = q_sor, colour = "blue", linetype = "dashed") +
  coord_cartesian(xlim = c(0, min(q_sor[2], 10))) +
  theme_light() +
  labs(x = "SOR", y = "Density", title = "SOR")

p17 <- ggplot(filter(ann, is.finite(MQRankSum)), aes(x = MQRankSum)) +
  dens_layer() +
  geom_vline(xintercept = -12.5, colour = "red") +
  geom_vline(xintercept = q_mqrank, colour = "blue", linetype = "dashed") +
  theme_light() +
  labs(x = "MQRankSum", y = "Density", title = "MQRankSum")

p18 <- ggplot(filter(ann, is.finite(ReadPosRankSum)), aes(x = ReadPosRankSum)) +
  dens_layer() +
  geom_vline(xintercept = -8, colour = "red") +
  geom_vline(xintercept = q_rprs, colour = "blue", linetype = "dashed") +
  theme_light() +
  labs(x = "ReadPosRankSum", y = "Density", title = "ReadPosRankSum")

variant_panel <- (p11 + p12) / (p13 + p14) / (p15 + p16) / (p17 + p18)

message("[", Sys.time(), "] Saving plots")
ggsave(file.path(out_dir, "invariant_sample_missingness.png"), p1, width = 8, height = 10, dpi = 300)
ggsave(file.path(out_dir, "invariant_sample_depth.png"), p2, width = 8, height = 10, dpi = 300)
ggsave(file.path(out_dir, "invariant_site_missingness.png"), p3, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "invariant_site_depth.png"), p4, width = 8, height = 6, dpi = 300)

ggsave(file.path(out_dir, "variant_sample_missingness.png"), p5, width = 8, height = 10, dpi = 300)
ggsave(file.path(out_dir, "variant_sample_depth.png"), p6, width = 8, height = 10, dpi = 300)
ggsave(file.path(out_dir, "variant_site_missingness.png"), p7, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "variant_site_depth.png"), p8, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "variant_site_quality.png"), p9, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "variant_sample_heterozygosity.png"), p10, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "variant_annotation_summaries.png"), variant_panel, width = 14, height = 18, dpi = 300)

png(file.path(out_dir, "variant_hardfilter_quantiles_table.png"), width = 1200, height = 500, bg = "white")
grid.table(quantile_table)
dev.off()

message("[", Sys.time(), "] plot_QC.R finished successfully")

