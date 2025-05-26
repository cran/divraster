## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE, 
  message = FALSE
)

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

