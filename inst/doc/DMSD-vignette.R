## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  require(devtools)
#  devtools::install_github("flaviomoc/divraster", build_vignettes = TRUE)

## ---- fig.height = 5, fig.width = 8, fig.align = 'center'---------------------
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

## ---- fig.height = 4, fig.width = 5, fig.align = 'center'---------------------
# Taxonomic
alpha.td <- divraster::spat.alpha(bin1)
alpha.td
terra::plot(alpha.td, main = names(alpha.td))

# Functional
alpha.fd <- divraster::spat.alpha(bin1, traits)
alpha.fd
terra::plot(alpha.fd, main = names(alpha.fd))

# Phylogenetic
alpha.pd <- divraster::spat.alpha(bin1, tree)
alpha.pd
terra::plot(alpha.pd, main = names(alpha.pd))

## ---- fig.height = 4, fig.width = 5, fig.align = 'center'---------------------
avg.traits1 <- divraster::spat.trait(bin1, traits)
avg.traits1
terra::plot(avg.traits1, main = names(avg.traits1))

## ---- fig.height = 4, fig.width = 5, fig.align = 'center'---------------------
# Functional
ses.fd <- divraster::spat.rand(bin1, traits, 3, "site")
ses.fd
terra::plot(ses.fd, main = names(ses.fd))

# Phylogenetic
ses.pd <- divraster::spat.rand(bin1, tree, 3, "site")
ses.pd
terra::plot(ses.pd, main = names(ses.pd))

## ---- fig.height = 4, fig.width = 5, fig.align = 'center'---------------------
# Taxonomic
beta.td <- divraster::spat.beta(bin1)
beta.td
terra::plot(beta.td, main = names(beta.td))

# Functional
beta.fd <- divraster::spat.beta(bin1, traits)
beta.fd
terra::plot(beta.fd, main = names(beta.fd))

# Phylogenetic
beta.pd <- divraster::spat.beta(bin1, tree)
beta.pd
terra::plot(beta.pd, main = names(beta.pd))

## ---- fig.height = 4, fig.width = 5, fig.align = 'center'---------------------
# Taxonomic
betatemp.td <- divraster::temp.beta(bin1, bin2)
betatemp.td
terra::plot(betatemp.td, main = names(betatemp.td))

# Functional
betatemp.fd <- divraster::temp.beta(bin1, bin2, traits)
betatemp.fd
terra::plot(betatemp.fd, main = names(betatemp.fd))

# Phylogenetic
betatemp.pd <- divraster::temp.beta(bin1, bin2, tree)
betatemp.pd
terra::plot(betatemp.pd, main = names(betatemp.pd))

