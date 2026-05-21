library(ggplot2)
library(scales)
library(tidyverse)

# paths
evec <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Picoides_wgs/out/1_pca/picoides_pca.eigenvec'
eval <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Picoides_wgs/out/1_pca/picoides_pca.eigenval'
meta <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Picoides_wgs/data/picoides_Meta.csv'


# --- read PCA + variance ---
eig  <- read.table(evec, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
vals <- scan(eval, quiet = TRUE)
pcvar <- 100 * vals / sum(vals)

# --- build join key ---
eig$ID_base <- sub("(_S\\d+)?(_L\\d+)?$", "", eig$IID)

# --- metadata ---
m <- read.csv(meta, check.names = FALSE, stringsAsFactors = FALSE)
m2 <- m[, c("sampleID","ecoregion","subspecies")]
names(m2) <- c("ID_base","ECO","SUBS")

# --- merge ---
df <- merge(eig, m2, by = "ID_base", all.x = TRUE)
df <- df[!is.na(df$ECO) & nzchar(trimws(df$ECO)) & df$ECO != "Unknown", , drop = FALSE]

df$ECO  <- droplevels(factor(df$ECO))
df$SUBS <- droplevels(factor(df$SUBS))

eco_levels  <- levels(df$ECO)
subs_levels <- levels(df$SUBS)

# labels
df$ID_label <- sub(".*-", "", df$ID_base)

# ggplot default hues
cols <- scales::hue_pal()(length(eco_levels))
names(cols) <- eco_levels

# axis labels
xlab <- sprintf("PC1 (%.2f%%)", pcvar[1])
ylab <- sprintf("PC2 (%.2f%%)", pcvar[2])

# ======================
# PCA Plot
# ======================
gg2d <- ggplot(df, aes(x = PC1, y = PC2)) +
  
  geom_point(
    aes(fill = ECO, shape = SUBS),
    size = 6,
    color = "black",
    stroke = 0.5
  ) +
  
  
  
  scale_fill_manual(values = cols, drop = FALSE) +
  
  scale_shape_manual(
    values = c(21, 24, 25)[seq_along(subs_levels)],
    drop = FALSE
  ) +
  
  labs(x = xlab, y = ylab, fill = "Ecoregion", shape = "Subspecies") +
  
  theme_light(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # 🔥 tight legend spacing
    legend.key.height = unit(0.2, "cm"),
    legend.spacing.y  = unit(0.08, "cm"),
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    plot.margin = margin(t = 5, r = 30, b = 15, l = 5)
  ) +
  
  guides(
    fill = guide_legend(
      title = "Ecoregion",
      ncol = 1,
      byrow = TRUE,
      override.aes = list(
        shape = 21,
        size = 2.5,
        stroke = 0.5,
        color = "black"
      )
    ),
    shape = guide_legend(
      title = "Subspecies",
      ncol = 1,
      byrow = TRUE,
      override.aes = list(
        fill = "grey80",
        size = 2.5,
        stroke = 0.5,
        color = "black"
      )
    )
  ) +
  
  coord_cartesian(clip = "off")

# show
gg2d


