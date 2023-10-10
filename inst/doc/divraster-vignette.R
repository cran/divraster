## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE, 
  message = FALSE
)

## ----fig.height = 4, fig.width = 4, fig.align = 'center'----------------------
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
e <- c(-41, -39, -15, -13)
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
alpha.td <- divraster::spat.alpha(bin1)
alpha.td
terra::plot(alpha.td, main = paste0(names(alpha.td), "_sce1"))

# Alpha TD calculation for scenario 2
alpha.td2 <- divraster::spat.alpha(bin2)
alpha.td2
terra::plot(alpha.td2, main = paste0(names(alpha.td2), "_sce2"))

# Difference in Alpha TD between scenarios
alpha.td2-alpha.td
terra::plot(alpha.td2-alpha.td, main = "Delta Alpha TD")

## ----fig.height = 4, fig.width = 4, fig.align = 'center'----------------------
alpha.fd <- divraster::spat.alpha(bin1, traits)
alpha.fd
terra::plot(alpha.fd, main = names(alpha.fd))

## ----fig.height = 4, fig.width = 4, fig.align = 'center'----------------------
# Alpha PD calculation
alpha.pd <- divraster::spat.alpha(bin1, tree)
alpha.pd
terra::plot(alpha.pd, main = names(alpha.pd))

## ----fig.height = 5, fig.width = 6, fig.align = 'center'----------------------
# SES FD calculation
ses.fd <- divraster::spat.rand(x = bin1, 
                               tree = traits, 
                               aleats = 3, 
                               random = "site")
ses.fd
terra::plot(ses.fd, main = names(ses.fd))

## ----fig.height = 5, fig.width = 6, fig.align = 'center'----------------------
# SES PD calculation
ses.pd <- divraster::spat.rand(x = bin1, 
                               tree = tree, 
                               aleats = 3, 
                               random = "site")
ses.pd
terra::plot(ses.pd, main = names(ses.pd))

## ----fig.height = 5, fig.width = 6, fig.align = 'center'----------------------
# Beta spatial TD calculation
beta.td <- divraster::spat.beta(bin1)
beta.td
terra::plot(beta.td, main = names(beta.td))

## ----fig.height = 5, fig.width = 6, fig.align = 'center'----------------------
# Beta spatial FD calculation
beta.fd <- divraster::spat.beta(bin1, traits)
beta.fd
terra::plot(beta.fd, main = names(beta.fd))

## ----fig.height = 5, fig.width = 6, fig.align = 'center'----------------------
# Beta spatial PD calculation
beta.pd <- divraster::spat.beta(bin1, tree)
beta.pd
terra::plot(beta.pd, main = names(beta.pd))

## ----fig.height = 5, fig.width = 6, fig.align = 'center'----------------------
# Beta temporal TD calculation
betatemp.td <- divraster::temp.beta(bin1, bin2)
betatemp.td
terra::plot(betatemp.td, main = names(betatemp.td))

## ----fig.height = 5, fig.width = 6, fig.align = 'center'----------------------
# Beta temporal FD calculation
betatemp.fd <- divraster::temp.beta(bin1, bin2, traits)
betatemp.fd
terra::plot(betatemp.fd, main = names(betatemp.fd))

## ----fig.height = 5, fig.width = 6, fig.align = 'center'----------------------
# Beta temporal PD calculation
betatemp.pd <- divraster::temp.beta(bin1, bin2, tree)
betatemp.pd
terra::plot(betatemp.pd, main = names(betatemp.pd))

## ----fig.height = 4, fig.width = 4, fig.align = 'center'----------------------
# Average traits calculation
# Scenario 1
avg.traits1 <- divraster::spat.trait(bin1, traits)
avg.traits1[[4]]
terra::plot(avg.traits1[[4]], main = paste0(names(avg.traits1[[4]]), "_sce1"))

# Scenario 2
avg.traits2 <- divraster::spat.trait(bin2, traits)
avg.traits2[[4]]
terra::plot(avg.traits2[[4]], main = paste0(names(avg.traits2[[4]]), "_sce2"))

# Percentage of change
change.traits <- (avg.traits2 - avg.traits1) / avg.traits1 * 100
change.traits[[4]]
terra::plot(change.traits[[4]], main = paste0(names(change.traits[[4]]), "_%"))

