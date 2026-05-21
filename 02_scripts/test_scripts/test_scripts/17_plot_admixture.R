# ----- Paths -----
base   <- "/usr/scratch2/userdata2/mtofflemire/projects/picoides_genotyping/out2"
famfn  <- "/usr/scratch2/userdata2/mtofflemire/projects/picoides_genotyping/out2/picoides.fam"
Qfn    <- "/home/mtofflemire/picoides.2.Q"   # change 2 → whichever K
outpng <- "/usr/scratch2/userdata2/mtofflemire/projects/picoides_genotyping/out2/picoides.K2.by_county.png"
outcsv <- "/usr/scratch2/userdata2/mtofflemire/projects/picoides_genotyping/out2/picoides.K2.scores.by_county.csv"
missfn <- file.path(base, "tokens_missing_county.txt")

# ----- Read ADMIXTURE results -----
Q <- read.table(Qfn, header = FALSE)  # no header in .Q
fam <- read.table(famfn, header = FALSE)
ids <- fam$V2                       # sample IDs

# attach sample IDs to Q matrix
Q$SampleID <- ids

# ----- Extract museum token from IID (same as PCA script) -----
Q$token <- sub(".*-([A-Za-z]+\\d+).*", "\\1", ids, perl = TRUE)

# ----- Token -> County map (reuse from PCA) -----
token_to_county <- c(
  "DMNS47995"  = "Crook",
  "LACM109434" = "Tulare",
  "LACM113646" = "Los Angeles",
  "SDSU2586"   = "Riverside",
  "SDSU2588"   = "Mendocino",
  "SDSU2589"   = "Mendocino",
  "SDSU2590"   = "Mendocino",
  "SDSU2591"   = "Trinity",
  "SDSU2592"   = "Trinity",
  "SDSU2593"   = "Shasta",
  "SDSU2594"   = "Shasta",
  "SDSU2595"   = "Shasta",
  "SDSU2596"   = "Shasta",
  "SDSU2601"   = "El Dorado",
  "SDSU2602"   = "El Dorado",
  "SDSU2603"   = "El Dorado",
  "SDSU2604"   = "El Dorado",
  "SDSU2606"   = "Kern",
  "SDSU2607"   = "Kern",
  "SDSU2608"   = "Kern",
  "SDSU2609"   = "Kern",
  "SDSU2610"   = "Kern",
  "SDSU2612"   = "Fresno",
  "SDSU2643"   = "San Bernardino",
  "SDSU2646"   = "San Bernardino",
  "SDSU2647"   = "San Bernardino",
  "SDSU2669"   = "Fresno",
  "SDSU2671"   = "Siskiyou",
  "SDSU2673"   = "Siskiyou",
  "SDSU2680"   = "Modoc",
  "UWBM51058"  = "Yakima",
  "UWBM57122"  = "Yakima"
)

# Assign counties; warn if any missing
Q$county <- unname(token_to_county[Q$token])
missing <- is.na(Q$county)
if (any(missing)) {
  writeLines(sort(unique(Q$token[missing])), missfn)
  message("WARNING: some tokens not in map. Wrote list to: ", missfn)
  Q$county[missing] <- "Unknown"
}

# ----- Reshape for plotting -----
K <- ncol(Q) - 3  # number of clusters (exclude SampleID, token, county)
colnames(Q)[1:K] <- paste0("Cluster", 1:K)

library(tidyr)
library(dplyr)
library(ggplot2)

dat <- Q %>%
  select(SampleID, county, starts_with("Cluster")) %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Cluster", values_to = "Ancestry")

# ----- Plot -----
p <- ggplot(dat, aes(x = SampleID, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(~ county, scales = "free_x", space = "free_x") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    strip.text = element_text(size = 9)
  ) +
  ylab("Ancestry proportion") +
  xlab("Individuals")

ggsave(outpng, p, width = 12, height = 4, dpi = 300)

# ----- Save CSV -----
write.csv(dat, outcsv, row.names = FALSE)
message("Wrote plot: ", outpng)
message("Wrote CSV:  ", outcsv)

