plot_dir <- "/Users/michaeltofflemire/softs/eems_clean/plotting/rEEMSplots/R"
files <- list.files(plot_dir, pattern = "\\.R$", full.names = TRUE)
invisible(lapply(files, source))

library(rEEMSplots)
getLoadedDLLs()
library(rEEMSplots)
is.loaded("_rEEMSplots_tiles2contours_standardize")


eems_results <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/03_analysis/04_eems'

plot_prefix <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/03_analysis/04_eems/figures'



outline_file <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/02_scripts/EEMS/eems.Leuconotopicus.outer'

library(raster)
install.packages(c(
  "sf",
  "rEEMSplots",
  "ggplot2",
  "fields",
  "sp"
))

eems.plots(
    mcmcpath    = eems_results,
    plotpath    = plot_prefix,
    longlat     = TRUE,
    add.outline = TRUE,
    out.png     = FALSE
)


install.packages(c("rworldmap", "rworldxtra", "sp", "raster"))



library(sp)
library(raster)
library(rworldmap)
library(rworldxtra)
install.packages("sf")   # only once
library(sf)



projection_longlat <- "+proj=longlat +datum=WGS84 +no_defs"
projection_mercator <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

eems.plots(
    mcmcpath      = eems_results,
    plotpath      = plot_prefix,
    longlat       = TRUE,
    
    add.map       = TRUE,
    col.map       = "black",
    lwd.map       = 1.5,
    add.demes = TRUE,
    
    projection.in  = projection_longlat,
    projection.out = projection_mercator,
    
    out.png       = FALSE
)
