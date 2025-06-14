#' Check if objects are valid
#'
#' @param bin1 A SpatRaster with presence-absence data (0 or 1)
#' for a set of species.
#' @param bin2 A SpatRaster with presence-absence data (0 or 1)
#' for a set of species. Species names in 'bin2' and 'bin1' must
#' match!
#' @param tree It can be a 'data.frame' with species traits or a
#' 'phylo' with a rooted phylogenetic tree. Species names in 'tree',
#' 'bin1', and 'bin2' must match!
#'
#' @return Either a success message or an error.
inputs_chk <- function(bin1, bin2, tree) {
  success_message <- "Awesome! All objects are compatible with
  the package."  # Initialize success message

  # Check if 'bin1' is a valid SpatRaster object with geographic
  # coordinates and at least 2 layers
  if (is.null(bin1) || !inherits(bin1, "SpatRaster") ||
      !terra::is.lonlat(bin1) || terra::nlyr(bin1) < 2) {
    stop("'bin1' must be a valid SpatRaster with geographic
         coordinates and at least 2 layers.")
  }

  # Check if 'tree' is missing
  if (missing(tree)) {
    # Check if 'bin2' is provided
    if (!missing(bin2)) {
      # Check if 'bin2' is a valid SpatRaster object with
      # geographic coordinates and at least 2 layers
      if (is.null(bin2) || !inherits(bin2, "SpatRaster") ||
          !terra::is.lonlat(bin2) || terra::nlyr(bin2) < 2) {
        stop("'bin2' must be a valid SpatRaster with geographic
             coordinates and at least 2 layers.")
      }

      # Check if species names in 'bin1' and 'bin2' match
      if (!identical(names(bin1), names(bin2))) {
        stop("Species names in 'bin1' and 'bin2' must match.")
      }
    }
  } else {
    # Check if 'tree' is a data.frame or a phylo object
    if (!inherits(tree, c("data.frame", "phylo"))) {
      stop("'tree' must be a data.frame or a phylo object.")
    }

    # Check if 'bin2' is missing
    if (missing(bin2)) {
      # Check if 'tree' is a data.frame
      if (inherits(tree, "data.frame")) {
        # Check if 'bin1' and 'tree' have matching species names
        if (!identical(names(bin1), rownames(tree))) {
          stop("Species names in 'bin1' and 'tree' must match.")
        }
      } else if (inherits(tree, "phylo")) {
        # Check if species names in 'bin1' and 'tree' match
        if (!identical(sort(names(bin1)), sort(tree$tip.label))) {
          stop("Species names in 'bin1' and 'tree' must match.")
        }
      }
    } else {
      # Check if 'bin2' is a valid SpatRaster object with
      # geographic coordinates and at least 2 layers
      if (is.null(bin2) || !inherits(bin2, "SpatRaster") ||
          !terra::is.lonlat(bin2) || terra::nlyr(bin2) < 2) {
        stop("'bin2' must be a valid SpatRaster with geographic
             coordinates and at least 2 layers.")
      }

      # Check if species names in 'bin1', 'bin2', and 'tree' match
      if (inherits(tree, "phylo")) {
        if (!identical(sort(names(bin1)), sort(names(bin2))) ||
            !identical(sort(tree$tip.label), sort(names(bin1)))) {
          stop("Species names in 'bin1', 'bin2', and 'tree' must
               match.")
        }
      } else if (inherits(tree, "data.frame")) {
        if (!identical(sort(names(bin1)), sort(names(bin2))) ||
            !identical(sort(rownames(tree)), sort(names(bin1)))) {
          stop("Species names in 'bin1', 'bin2', and 'tree'
               must match.")
        }
      }
    }
  }

  # If no errors occurred, return the success message
  return(success_message)
}
