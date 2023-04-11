#########################################################################
## This R scrip contains all nessassery code for 
## plotting results of LLE and Cr in each SLAs 


## created 2023-3-31 by Yuxin Huang 
##########################################################################

## loading pekages
library(haven)
library(dplyr)
library(ggplot2)
library(spData)
library(spdep)
library(rgdal)
library(sf)
library(terra)
library(rgeos)
library(maptools)
library(sp)
library(scales)
library(gridExtra)
## loading data
lle <- read_dta("C:/Users/n11117761/Work/Data/FinalResult/lle_est_gamma05.dta")

cif <- read_dta("C:/Users/n11117761/Work/Data/FinalResult/cif_est_gamma05.dta")

## prepare data
# calculate the average estimation in each SLA
lle_avg <- lle %>%
  group_by(grid) %>%
  summarise_at(c("exp","lle_2","lle_50","lle_97","obs_2","obs_50","obs_97"), mean) %>%
  rename(id=grid)


cif_avg <- cif %>%
  group_by(grid) %>%
  summarise_at(c("cr_2","netm_2","cr_50","netm_50","cr_97","netm_97"), mean) %>%
  rename(id=grid)

## loading Shapfile: QLD SLA 2006
map_qld_06  <- readOGR("C:/Users/n11117761/Work/Data/MapInfo file/SLA_QLD_06/SLA_QLD_06.shp"
)

## prepare map data  for plotting
map_df <- fortify(map_qld_06)            # as data frame object(for merging with estimation)
map_df$id <- as.numeric(map_df$id) + 1

map.border <- unionSpatialPolygons(map_qld_06, IDs = rep(1,478)) #map of border (for ggplot) 
map.border <- fortify(map.border)

## creat sub-map if Brisbane

map.proj <- proj4string(map_qld_06)

# User-defined function for creating inset windows
get.inset <- function(x1, x2, y1, y2, map_qld_06, map.proj){
  Inset.window <- as(raster::extent(x1, x2, y1, y2), "SpatialPolygons")
  proj4string(Inset.window) <- map.proj
  map.inset <- gIntersection(map_qld_06, Inset.window, byid = TRUE, drop_lower_td = TRUE, 
                             id = sapply(map_qld_06@polygons, function(x) x@ID))
  return(map.inset)
}

map_brisbane <- get.inset(152.6, 153.6, -28, -27, map_qld_06, map.proj)

map.border.brisbane <- unionSpatialPolygons(map_brisbane, IDs = rep(1, length(map_brisbane)))
map_brisbane_df <- fortify(map_brisbane)
map_brisbane_df$id <- as.numeric(map_brisbane_df$id) + 1

################################ map of median LLE ######################################

## define color and value scales
Fill.colours.lle <- c("#d4e6f5", "#abcfeb", "#96c3e6", "#81b7e1", "#276da1", "#1d5178","#0e2b3a")

# Values corresponding to the fill colours

Fill.values.lle.p <- c(0.3,10,15,20,30,40)
Fill.values.lle.p.r <- rescale(Fill.values.lle.p, from = range(as.numeric(lle_avg$lle_50)))

#rescale values
Fill.values.lle.exp.r <- rescale(Fill.values.lle.exp, from = range(as.numeric(lle_avg$exp)))


## merge values with shapfile dataframe
lle_grid <- inner_join(map_df,lle_avg,by = "id") 
lle_brisbane <- inner_join(map_brisbane_df, lle_avg, by = "id")

## base map
gg.base <- ggplot(data = NULL, aes(x = long, y = lat, group = group)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() + 
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(barwidth = 15)) +
  coord_map() 

## layer map
gg.lle50 <- gg.base + 
  geom_polygon(data = map.border, fill = "grey30", color = "black", size = 0.8) +
  geom_polygon(data = lle_grid, aes(fill = lle_50), color = NA) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.lle, 
                       values = Fill.values.lle.p.r,
                       limits = range(lle_avg$lle_50)
  ) +
  ggtitle("Loss of Life Expectancy") +
  annotate("rect",
           xmin = 152.6 - 0.01, xmax = 153.6 + 0.01,
           ymin = -28 - 0.01, ymax = -27 + 0.01,
           colour = "red", fill = NA,size = 1)

## sub-map of Brisbane
gg.base.inset <- gg.base + guides(fill = FALSE, alpha = FALSE) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.lle, 
                       values = Fill.values.lle.p.r,
                       limits = range(lle_avg$lle_50))
gg.inset.brisbane50 <- gg.base.inset + 
  geom_polygon(data = map.border.brisbane, fill = "grey90", color = "black", size = 0.8) +
  geom_polygon(data = lle_brisbane, aes(fill = lle_50), color = NA) +
  ggtitle("Inner city, 50%")

## arrange
png(filename = "C:/Users/n11117761/Work/Data/FinalResult/lle50_R.png")
grid.arrange(
  gg.lle50,
  gg.inset.brisbane50,
  nrow = 1,
  ncol = 2,
  widths = c(10,5)
)
dev.off()

############################### map of Cr #######################################
Fill.colours.cif <- c("#fff7c4", "#ffffc4","#ffffd4","#fee391", "#fec44f","#ff9e33","#d95f0e", "#be3d00", "#993404")

Fill.values.cif.p <- c(0.1,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.85)
Fill.values.cif.p50.r <- rescale(Fill.values.cif.p, from = range(as.numeric(cif_avg$cr_50)))

cif_grid <- inner_join(map_df,cif_avg,by = "id") 

gg.cif50 <- gg.base + 
  geom_polygon(data = map.border, fill = "grey90", color = "black", size = 0.8) +
  geom_polygon(data = cif_grid, aes(fill = cr_50), color = NA) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.cif, 
                       values = Fill.values.cif.p50.r,
                       limits = range(cif_avg$cr_50)
  ) +
  ggtitle("Crude Mortality, 50%") +
  annotate("rect",
           xmin = 152.6 - 0.01, xmax = 153.6 + 0.01,
           ymin = -28 - 0.01, ymax = -27 + 0.01,
           colour = "black", fill = NA,size = 1)

gg.base.inset.cif50 <- gg.base + guides(fill = FALSE, alpha = FALSE) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.cif, 
                       values = Fill.values.cif.p50.r,
                       limits = range(cif_avg$cr_50))

gg.inset.brisbane50 <- gg.base.inset.cif50 + 
  geom_polygon(data = map.border.brisbane, fill = "grey90", color = "black", size = 0.8) +
  geom_polygon(data = cif_brisbane, aes(fill = cr_50), color = NA) +
  ggtitle("Inner city")

png(filename = "C:/Users/n11117761/Work/Data/FinalResult/cif50_R.png")
grid.arrange(
  gg.cif50,
  gg.inset.brisbane50,
  nrow = 1,
  ncol = 2,
  widths = c(10,5)
)
dev.off()

