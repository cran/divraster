## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE, 
  message = FALSE
)

## -----------------------------------------------------------------------------
library(divraster)
library(terra)
library(dplyr)
library(sf)

## -----------------------------------------------------------------------------
# Loading data
# Presence-absence SpatRaster
bin1 <- terra::rast(system.file("extdata", 
                                "ref_frugivor.tif", 
                                package = "divraster"))
bin2 <- terra::rast(system.file("extdata",
                                "fut_frugivor.tif",
                                package = "divraster"))

# Change extension to process faster
terra::ext(bin1)
e <- c(-40, -39, -14, -13)
bin1 <- terra::crop(bin1, e)
bin2 <- terra::crop(bin2, e)

# Species traits
traits <- read.csv(system.file("extdata", 
                               "traits_frugivor.csv", 
                               package = "divraster"), 
                   sep = ";", 
                   row.names = 1)

# Phylogenetic tree
tree <- ape::read.tree(system.file("extdata", 
                                   "tree_frugivor.tre", 
                                   package = "divraster"))

# Alpha TD calculation for scenario 1
divraster::spat.alpha(bin1)

# Alpha TD calculation for scenario 2
divraster::spat.alpha(bin2)

## ----eval = FALSE-------------------------------------------------------------
# divraster::spat.alpha(bin1, traits)

## -----------------------------------------------------------------------------
# Alpha PD calculation
divraster::spat.alpha(bin1, tree)

## ----eval = FALSE-------------------------------------------------------------
# # SES FD calculation
# divraster::spat.rand(x = bin1,
#                      tree = traits,
#                      aleats = 3,
#                      random = "site")

## ----eval = FALSE-------------------------------------------------------------
# # SES PD calculation
# divraster::spat.rand(x = bin1,
#                      tree = tree,
#                      aleats = 3,
#                      random = "site")

## ----eval = FALSE-------------------------------------------------------------
# # Beta spatial TD calculation
# divraster::spat.beta(bin1)

## ----eval = FALSE-------------------------------------------------------------
# # Beta spatial FD calculation
# divraster::spat.beta(bin1, traits)

## ----eval = FALSE-------------------------------------------------------------
# # Beta spatial PD calculation
# divraster::spat.beta(bin1, tree)

## -----------------------------------------------------------------------------
# Beta temporal TD calculation
divraster::temp.beta(bin1, bin2)

## -----------------------------------------------------------------------------
# Beta temporal FD calculation
divraster::temp.beta(bin1, bin2, traits)

## -----------------------------------------------------------------------------
# Beta temporal PD calculation
divraster::temp.beta(bin1, bin2, tree)

## -----------------------------------------------------------------------------
# Average traits calculation
# Scenario 1
divraster::spat.trait(bin1, traits)

# Scenario 2
divraster::spat.trait(bin2, traits)

## -----------------------------------------------------------------------------
# Suitability change between climate scenarios
change <- divraster::suit.change(bin1[[1:4]], bin2[[1:4]])

# Visualization
# Initialize an empty list to store color mappings for each layer
change.col <- list()

# Loop through each layer in the 'change' raster
for(i in 1:terra::nlyr(change)){
  # Get unique values, omitting NA values, and sort them
  change.col[[i]] <- sort(na.omit(unique(terra::values(change[[i]]))))
  
  # Map numeric values to specific colors
  change.col[[i]][change.col[[i]] == 1] <- "blue"   # Gain
  change.col[[i]][change.col[[i]] == 2] <- "red"    # Loss
  change.col[[i]][change.col[[i]] == 3] <- "grey"   # No change
  change.col[[i]][change.col[[i]] == 4] <- "white"  # Unsuitable
  
  # Plot the change in suitability
  terra::plot(change[[i]], col = change.col[[i]], 
              main = names(change)[i],
              cex.main = 1,          
              font.main = 4)  
}

## -----------------------------------------------------------------------------
# Climate suitable area scenario 1
divraster::area.calc(bin1[[1:4]])

# Climate suitable area scenario 2
divraster::area.calc(bin2[[1:4]])

## -----------------------------------------------------------------------------
# Difference in species richness between climate scenarios
divraster::differ.rast(divraster::spat.alpha2(bin1),
                       divraster::spat.alpha2(bin2), perc = FALSE)

## -----------------------------------------------------------------------------
set.seed(1)

occurrences <- data.frame(
species = c(
rep("Species_A", 3),
rep("Species_B", 4),
"Species_C"
),
lon = c(-40, -41, -42, -43, -44, -45, -46, -47),
lat = c(-20, -21, -22, -23, -24, -25, -26, -27)
)

avg_dist <- occ.avg.dist(occurrences)
avg_dist

## -----------------------------------------------------------------------------
# Continuous raster (e.g. environmental suitability)
r_env <- rast(
ncol = 50, nrow = 50,
xmin = 0, xmax = 10,
ymin = 0, ymax = 10)
values(r_env) <- runif(ncell(r_env), 0, 1)
names(r_env) <- "suitability"

# Two polygons with IDs
p1 <- vect("POLYGON ((0 0, 5 0, 5 10, 0 10, 0 0))", crs = "EPSG:4326")
p2 <- vect("POLYGON ((5 0, 10 0, 10 10, 5 10, 5 0))", crs = "EPSG:4326")
polys <- rbind(p1, p2)
polys$poly_id <- c("P1", "P2")

# Mean suitability per polygon
poly_mean <- rast.by.polys(
x = r_env,
polygons = polys,
id_col = "poly_id"
)
poly_mean

## -----------------------------------------------------------------------------
poly_stats <- rast.by.polys(
x = r_env,
polygons = polys,
id_col = "poly_id",
fun = function(v, ...) c(
mean = mean(v, ...),
min = min(v, ...),
max = max(v, ...)
),
na.rm = TRUE
)
poly_stats

## -----------------------------------------------------------------------------
#Continuous raster
r_cont <- rast(
ncol = 50, nrow = 50,
xmin = 0, xmax = 10,
ymin = 0, ymax = 10)
values(r_cont) <- runif(ncell(r_cont), 0, 1)
names(r_cont) <- "suitability"

#Binary footprint (circle around centre)
r_bin <- rast(r_cont)
xy <- terra::xyFromCell(r_bin, 1:ncell(r_bin))
dist_center <- sqrt((xy[, 1] - 5)^2 + (xy[, 2] - 5)^2)
values(r_bin) <- ifelse(dist_center <= 3, 1, 0)
names(r_bin) <- "footprint"

#Crop continuous raster to footprint
r_cropped <- bin2crop(r_bin = r_bin, r_cont = r_cont)

r_cropped

par(mfrow = c(1, 3))
plot(r_cont, main = "Original raster")
plot(r_bin, main = "Binary footprint")
plot(r_cropped, main = "Cropped raster")
par(mfrow = c(1, 1))

