# ----- Paths -----
base   <- "/usr/scratch2/userdata2/mtofflemire/projects/picoides_genotyping/out2"
evec   <- file.path(base, "picoides.pca.eigenvec")   # tab-delimited with header
evals  <- file.path(base, "picoides.pca.eigenval")
outpng <- file.path(base, "picoides.pca.PC1_PC2.by_county.png")
outcsv <- file.path(base, "picoides.pca.scores.by_county.csv")
missfn <- file.path(base, "tokens_missing_county.txt")

# ----- Read PCA -----
ev  <- read.table(evec, header = TRUE, sep = "\t", check.names = FALSE)
lam <- scan(evals, quiet = TRUE)
pvar <- lam / sum(lam) * 100
pc1_lab <- sprintf("PC1 (%.2f%%)", pvar[1])
pc2_lab <- sprintf("PC2 (%.2f%%)", pvar[2])

# ----- Extract museum token from IID (e.g., ...-SDSU2643_... -> SDSU2643) -----
ev$token <- sub(".*-([A-Za-z]+\\d+).*", "\\1", ev$IID, perl = TRUE)

# ----- Token -> County map (from your list) -----
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

# Assign counties; report any tokens not found
ev$county <- unname(token_to_county[ev$token])
missing <- is.na(ev$county)
if (any(missing)) {
  writeLines(sort(unique(ev$token[missing])), missfn)
  message("WARNING: some tokens were not in the map. Wrote list to: ", missfn)
  ev$county[missing] <- "Unknown"
}

# ----- Colors -----
uniq <- sort(unique(ev$county))
pal  <- grDevices::hcl.colors(length(uniq), "Dark 3")
names(pal) <- uniq
cols <- pal[ev$county]

# ----- Plot (legend bottom-right) -----
png(outpng, width = 1800, height = 1400, res = 220)
par(mar = c(4.6, 4.6, 1, 1))
plot(ev$PC1, ev$PC2, pch = 19, cex = 1.25, col = cols,
     xlab = pc1_lab, ylab = pc2_lab)
grid()
legend("bottomright",
       legend = names(pal),
       col    = unname(pal),
       pch    = 19,
       cex    = 0.95,
       pt.cex = 1.2,
       inset  = 0.02,
       bty    = "o",
       bg     = "white",
       box.lwd= 0.6)
dev.off()

# ----- Save CSV with county -----
write.csv(ev, outcsv, row.names = FALSE)
message("Wrote plot: ", outpng)
message("Wrote CSV:  ", outcsv)

