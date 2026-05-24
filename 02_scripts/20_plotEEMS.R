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















library(sp)
library(raster)
library(rworldmap)
library(rworldxtra)
library(sf)
library(rEEMSplots)

projection_longlat <- "+proj=longlat +datum=WGS84 +no_defs"
projection_mercator <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

eems.plots(
  mcmcpath      = eems_results,
  plotpath      = plot_prefix,
  
  longlat       = TRUE,
  
  add.map       = TRUE,
  col.map       = "white",
  lwd.map       = 2,
  
  add.demes     = TRUE,
  
  projection.in  = projection_longlat,
  projection.out = projection_mercator,
  
  out.png       = FALSE,
  
  add.grid      = TRUE,
  col.grid      = "gray40",
  lwd.grid      = 0.5,
  max.cex.demes = 2
)








library(sp)
library(raster)
library(rworldmap)
library(rworldxtra)
library(sf)
library(rEEMSplots)

projection_longlat <- "+proj=longlat +datum=WGS84 +no_defs"
projection_mercator <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

eems.plots(
  mcmcpath      = eems_results,
  plotpath      = plot_prefix,
  
  longlat       = TRUE,
  
  add.map       = TRUE,
  col.map       = "white",
  lwd.map       = 2,
  
  add.demes     = TRUE,
  
  projection.in  = projection_longlat,
  projection.out = projection_mercator,
  
  out.png       = FALSE,
  
  add.grid      = TRUE,
  col.grid      = "gray40",
  lwd.grid      = 0.5,
  max.cex.demes = 2
)








eems.plots(
  mcmcpath      = eems_results,
  plotpath      = plot_prefix,
  longlat       = TRUE,
  
  add.map       = TRUE,
  col.map       = "black",
  lwd.map       = 3,
  
  add.demes     = TRUE,
  
  projection.in = projection_longlat,
  
  out.png       = FALSE,
  
  add.grid      = TRUE,
  col.grid      = "gray40",
  lwd.grid      = 0.5,
  max.cex.demes = 2
)














library(rEEMSplots)
library(maps)

projection_longlat <- "+proj=longlat +datum=WGS84 +no_defs"
projection_mercator <- "+proj=merc +datum=WGS84"

eems.plots(
  mcmcpath = eems_results,
  plotpath = plot_prefix,
  
  longlat = TRUE,
  
  projection.in = projection_longlat,
  projection.out = projection_mercator,
  
  add.demes = TRUE,
  
  add.grid = TRUE,
  col.grid = "gray40",
  lwd.grid = 0.5,
  max.cex.demes = 2,
  
  m.plot.xy = {
    map(
      "state",
      col = "black",
      lwd = 1,
      add = TRUE
    )
  },
  
  q.plot.xy = {
    map(
      "state",
      col = "black",
      lwd = 1,
      add = TRUE
    )
  },
  
  out.png = FALSE
)











library(rEEMSplots)
library(maps)

projection_longlat <- "+proj=longlat +datum=WGS84 +no_defs"
projection_mercator <- "+proj=merc +datum=WGS84"

eems.plots(
  mcmcpath = eems_results,
  plotpath = plot_prefix,
  
  longlat = TRUE,
  
  # KEEP COUNTRY / COAST BORDERS
  add.map = TRUE,
  col.map = "black",
  lwd.map = 1.2,
  
  # ADD STATE BORDERS
  m.plot.xy = {
    map(
      "state",
      col = "black",
      lwd = 0.6,
      add = TRUE
    )
  },
  
  q.plot.xy = {
    map(
      "state",
      col = "black",
      lwd = 0.6,
      add = TRUE
    )
  },
  
  add.demes = TRUE,
  
  projection.in = projection_longlat,
  projection.out = projection_mercator,
  
  out.png = FALSE,
  
  add.grid = TRUE,
  col.grid = "gray40",
  lwd.grid = 0.5,
  max.cex.demes = 2
)






library(rEEMSplots)
library(maps)

projection_longlat <- "+proj=longlat +datum=WGS84 +no_defs"

eems.plots(
  mcmcpath = eems_results,
  plotpath = plot_prefix,
  
  longlat = TRUE,
  
  add.map = TRUE,
  col.map = "black",
  lwd.map = 1.7,
  
  add.demes = TRUE,
  
  projection.in = projection_longlat,
  
  out.png = FALSE,
  
  add.grid = FALSE,
  col.grid = "gray40",
  max.cex.demes = 3,
  
  m.plot.xy = {
    map(
      "state",
      col = "black",
      lwd = 1,
      add = TRUE
    )
  },
  
  q.plot.xy = {
    map(
      "state",
      col = "black",
      lwd = 1,
      add = TRUE
    )
  }
)


















projection_longlat <- "+proj=longlat +datum=WGS84 +no_defs"

projection_longlat <- "+proj=longlat +datum=WGS84 +no_defs"

eems.plots(
  mcmcpath = eems_results,
  plotpath = plot_prefix,
  
  longlat = TRUE,
  
  # FIGURE SIZE
  plot.width = 7,
  plot.height = 9,
  add.title=FALSE,
  
  add.map = TRUE,
  col.map = "black",
  lwd.map = 2.2,
  
  add.demes = TRUE,
  
  projection.in = projection_longlat,
  
  out.png = FALSE,
  
  add.grid = FALSE,
  col.grid = "gray40",
  
  max.cex.demes = 3,
  
  m.plot.xy = {
    map(
      "state",
      col = "black",
      lwd = 1.5,
      add = TRUE
    )
  },
  
  q.plot.xy = {
    map(
      "state",
      col = "black",
      lwd = 1.5,
      add = TRUE
    )
  }
)







library(rEEMSplots)
library(maps)

projection_longlat <- "+proj=longlat +datum=WGS84 +no_defs"

eems.plots(
  mcmcpath = eems_results,
  plotpath = plot_prefix,
  
  longlat = TRUE,
  
  add.title = FALSE,
  
  plot.width = 6,
  plot.height = 9,
  
  add.map = TRUE,
  col.map = "black",
  lwd.map = 3.2,
  
  add.demes = TRUE,
  
  projection.in = projection_longlat,
  
  out.png = FALSE,
  
  # Gray grid
  add.grid = TRUE,
  col.grid = "gray70",
  lwd.grid = 1.4,
  
  max.cex.demes = 3
)
