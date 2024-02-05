##' Generate polygon based on alpha hulls
##' 
##' Generates an apha hull polygon, where the alpha parameter is determined by
##' the spatial distribution of the coordinates.
##' 
##' 
##' From a set of coordinates, this function will create an alpha hull with
##' \code{alpha = initialAlpha}, and will then increase \code{alpha} by
##' \code{alphaIncrement} until both the \code{fraction} and \code{partCount}
##' conditions are met.
##' 
##' If the conditions cannot be satisfied, then a minimum convex hull is
##' returned.
##' 
##' If \code{clipToCoast} is set to "terrestrial" or "aquatic", the resulting
##' polygon is clipped to the coastline, using a basemap from naturalearth. 
##' The first time this function is run, this basemap will be downloaded. 
##' Subsequent calls will use the downloaded map.
##' 
##' @param x dataframe of coordinates in decimal degrees, with a minimum of 3
##' rows. 
##' @param fraction the minimum fraction of occurrences that must be included
##' in polygon.
##' @param partCount the maximum number of disjunct polygons that are allowed.
##' @param buff buffering distance in meters
##' @param initialAlpha the starting value for alpha
##' @param coordHeaders the column names for the longitude and latitude
##' columns, respectively.  If x has two columns, these are assumed to be
##' longitude and latitude, and \code{coordHeaders} is ignored.
##' @param clipToCoast Either "no" (no clipping), "terrestrial" (only
##' terrestrial part of range is kept) or "aquatic" (only non-terrestrial part
##' is clipped). See Details.
##' @param alphaIncrement the amount to increase alpha with each iteration
##' @param verbose prints the alpha value to the console, intended for
##' debugging.
##' @param alphaCap Max alpha value before function aborts and returns a
##' minimum convex hull.
##' @return a list with 2 elements: \item{hull}{ a sf polygon object }
##' \item{alpha}{ the alpha value that was found to satisfy the criteria.  If a
##' convex hull was returned, this will list MCH.  }
##' @author Pascal Title
##' @seealso Alpha hulls are created with \code{\link{ahull}}.
##' @examples
##' 
##' data(crotalus)
##' 
##' # create a polygon range for Crotalus atrox
##' x <- crotalus[which(crotalus$genSp == 'Crotalus_atrox'),]
##' x <- x[sample(1:nrow(x), 50),]
##' 
##' range <- getDynamicAlphaHull(x, coordHeaders=c('decimallongitude','decimallatitude'), 
##' 	clipToCoast = 'no')
##' 
##' plot(range[[1]], col=transparentColor('dark green', 0.5), border = NA)
##' points(x[,c('decimallongitude','decimallatitude')], cex = 0.5, pch = 3)
##' 
##' # to add a basic coastline, you can use the internal map
##' # world <- rangeBuilder:::loadWorldMap()
##' # plot(world, add = TRUE, lwd = 0.5)
##' 
##' @export
getDynamicAlphaHull <- function (x, fraction = 0.95, partCount = 3, buff = 10000, initialAlpha = 3, 
                                 coordHeaders = c("Longitude", "Latitude"), clipToCoast = "terrestrial", 
                                 alphaIncrement = 1, verbose = FALSE, alphaCap = 400) {
  if (clipToCoast == FALSE) 
    clipToCoast <- "no"
  clipToCoast <- match.arg(clipToCoast, c("no", "terrestrial","aquatic"))
  if (ncol(x) == 2) {
    coordHeaders <- c(1, 2)
  }
  x <- x[!duplicated(x[, coordHeaders]), coordHeaders]
  x <- x[stats::complete.cases(x), ]
  if (nrow(x) < 3) {
    stop("This function requires a minimum of 3 unique coordinates (after removal of duplicates).")
  }
  while ((x[1, coordHeaders[1]] == x[2, coordHeaders[1]] & 
          x[2, coordHeaders[1]] == x[3, coordHeaders[1]]) | 
         (x[1,2] == x[2, coordHeaders[2]] & x[2, coordHeaders[2]] == x[3, coordHeaders[2]])) {
    x <- x[sample(1:nrow(x), size = nrow(x)), ]
  }
  
  x <- sf::st_as_sf(as.data.frame(x), coords = 1:2, crs = 4326)
  if (nrow(x) < 3) {
    stop("This function requires a minimum of 3 unique coordinates.")
  }
  
  alpha <- initialAlpha
  problem <- FALSE
  if (verbose) 
    cat("\talpha:", alpha, "\n")
  hull <- try(alphahull::ahull(sf::st_coordinates(x), alpha = alpha), 
              silent = TRUE)
  dropPt <- c()
  while (inherits(hull, "try-error") & any(grepl("duplicate points", hull))) {
    ptDist <- sf::st_distance(x, x)
    diag(ptDist) <- NA
    units(ptDist) <- NULL
    closest <- which(ptDist == min(ptDist, na.rm = TRUE), 
                     arr.ind = TRUE)
    hull <- try(alphahull::ahull(sf::st_coordinates(x)[-closest[1,1]], alpha = alpha),
                silent = TRUE)
    if (inherits(hull, "ahull")) {
      dropPt <- closest[1, 1]
    }
  }
  if (length(dropPt) > 0) {
    x <- x[-dropPt, ]
  }
  hull <- try(alphahull::ahull(sf::st_coordinates(x), alpha = alpha), 
              silent = TRUE)
  while (inherits(hull, "try-error")) {
    if (verbose) {
      cat("\talpha:", alpha, "\n")
    }
    # add a jitter to the coordinates to avoid the duplicated points error message
    x<- sf::st_jitter(x)
    alpha <- alpha + alphaIncrement
    hull <- try(alphahull::ahull(sf::st_coordinates(x), alpha = alpha), 
                silent = TRUE)
    if (alpha > alphaCap) {
      problem <- TRUE
      break
    }
  }
  if (!problem) {
    hull <- try(AlphaHullRangeModeller::ah2sf(hull),silent=TRUE)
    validityCheck <- function(hull) {
      if (!is.null(hull) & !inherits(hull, "try-error")) {
        if (!all(sf::st_is_valid(hull))) {
          TRUE
        }
        else {
          FALSE
        }
      }
      else {
        FALSE
      }
    }
    while (is.null(hull) | inherits(hull, "try-error") | 
           validityCheck(hull)) {
      alpha <- alpha + alphaIncrement
      if (verbose) {
        cat("\talpha:", alpha, "\n")
      }
      hull <- try(AlphaHullRangeModeller::ah2sf(alphahull::ahull(sf::st_coordinates(x), 
                                         alpha = alpha)), silent = TRUE)
      if (alpha > alphaCap) {
        problem <- TRUE
        break
      }
    }
    pointWithin <- sf::st_intersects(x, hull)
    alphaVal <- alpha
    buffered <- FALSE
    buff <- units::set_units(buff, "m")
    while (any(length(hull) > partCount, (sum(lengths(pointWithin))/nrow(x)) < 
               fraction, !all(sf::st_is_valid(hull)))) {
      alpha <- alpha + alphaIncrement
      if (verbose) {
        cat("\talpha:", alpha, "\n")
      }
      hull <- try(alphahull::ahull(sf::st_coordinates(x), 
                                   alpha = alpha), silent = TRUE)
      while (inherits(hull, "try-error") & alpha <= alphaCap) {
        alpha <- alpha + alphaIncrement
        hull <- try(alphahull::ahull(sf::st_coordinates(x), 
                                     alpha = alpha), silent = TRUE)
      }
      if (!inherits(hull, "try-error")) {
        hull <- AlphaHullRangeModeller::ah2sf(hull)
        hull <- sf::st_transform(hull, crs = "+proj=eqearth")
        if (all(sf::st_is_valid(hull))) {
          hull <- sf::st_buffer(hull, dist = buff)
          hull <- sf::st_transform(hull, crs = 4326)
          buffered <- TRUE
          pointWithin <- sf::st_intersects(x, hull)
        }
      }
      alphaVal <- alpha
      if (alpha > alphaCap) {
        hull <- sf::st_convex_hull(sf::st_union(x))
        hull <- sf::st_transform(hull, crs = "+proj=eqearth")
        hull <- sf::st_buffer(hull, dist = buff)
        hull <- sf::st_transform(hull, crs = 4326)
        buffered <- TRUE
        alphaVal = "MCH"
        break
      }
    }
  }
  else {
    hull <- sf::st_convex_hull(sf::st_union(x))
    hull <- sf::st_transform(hull, crs = "+proj=eqearth")
    hull <- sf::st_buffer(hull, dist = buff)
    hull <- sf::st_transform(hull, crs = 4326)
    buffered <- TRUE
    alphaVal = "MCH"
  }
  if (!buffered) {
    hull <- sf::st_transform(hull, crs = "+proj=eqearth")
    hull <- sf::st_buffer(hull, dist = buff)
    hull <- sf::st_transform(hull, crs = 4326)
  }
  if (clipToCoast != "no") {
    world <- rangeBuilder:::loadWorldMap()
    if (clipToCoast == "terrestrial") {
      hull <- sf::st_intersection(hull, world)
    }
    else {
      hull <- sf::st_difference(hull, world)
    }
  }
  return(list(hull, alpha = paste("alpha", alphaVal, sep = "")))
}


