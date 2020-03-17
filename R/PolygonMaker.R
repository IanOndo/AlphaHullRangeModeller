#' Generate an alpha-hull determine by the spatial distribution of points.
#'
#' Make alpha-hulls from a set of coordinates. PolygonMaker will try maxIter times to create an alpha hull with alpha = initialAlpha, and will then increase alpha by alphaIncrement
#  until both the fraction and partCount conditions are met.If the conditions cannot be satisfied, then a minimum convex hull is returned. If this process fails, it will create circular buffers around each occurrence points.
#'
#' @param x A two-column data.frame of occurrence records coordinates of the species.
#' @param fraction A numeric between 0 and 1 specifying the minimum fraction of occurrences that must be included in the polygon.
#' @param partCount A numeric integer specifying the maximum number of disjunct polygons that are allowed.
#' @param buffer A numeric specifying the buffering distance in meters around each point.
#' @param initialAlpha A numeric specifying the starting value for the parameter alpha.
#' @param coordHeaders A character string vector of length 2 specifying the names of the columns indicating the longitude and latitude respectively.
#' @param clipToCoast A logical. Should the terrestrial \code{clipToCoast='terrestrial'}, the aquatic part \code{clipToCoast='aquatic'} or none part \code{clipToCoast='no'} part of the range be kept ?
#' @param proj A character string specifying the projection of the coordinates
#' @param alphaIncrement A numeric specifying the amount to increase alpha at each iteration
#' @param alphaDecrement A numeric specifying the amount to decrease alpha if the function fails after the first iterations
#' @param maxIter A numeric integer specifying the maximum number of iterations before trying to build circular buffers.
#' @param other_buffers A numeric value or a vector of values of length \code{maxIter} specifying alternative buffer sizes (in meters) if the function fails after the first iterations.
#' @export
PolygonMaker <- function(x, coordHeaders = NULL, fraction=0.95, partCount=10, buffer, initialAlpha=2, alphaIncrement=1, alphaDecrement=1, maxIter=2, other_buffers=c(2000), clipToCoast="terrestrial", coastline=NULL, proj='+proj=longlat +datum=WGS84'){

  require(rangeBuilder)
  if(!inherits(x,"data.frame"))
    stop("Argument x must be a data.frame")
  if(buffer<=0)
    stop("Argument 'default_buffer' must be > 0")
  do.clipping = !is.null(coastline) && inherits(coastline, "sf")
  if(do.clipping){
    # The default clipping is provided by the Global Self-Consistent Hierarchical High-resolution Geography Database (GSHHG) via the package 'rangeBuilder')
    # But, high land basemap resolution may help to avoid crashes of the script
    clipToCoast="no"
  }
  #---------------------------------------------------------------------------
  #= 1. get coords (with a regular expression to look for longitude, latitude)
  #---------------------------------------------------------------------------
  if(is.null(coordHeaders)|length(coordHeaders)<2){
    x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(x),value=TRUE)[1]
    y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(x),value=TRUE)[1]
    coordHeaders <- c(x_lon,y_lat)
  }
  x  <- x[, coordHeaders]
  ll <- sf::st_as_sf(x, coords=c(1,2), crs=sf::st_crs(proj))
  # Do we have at least 3 points?
  if (lengths(ll) >= 3) {
    # Try making polygons
    made.polygons	<- tryCatch({
      made.polygons <- rangeBuilder::getDynamicAlphaHull(
        x,
        coordHeaders = coordHeaders,
        fraction = fraction,
        partCount = partCount,
        buff = buffer,
        initialAlpha = initialAlpha,
        alphaIncrement = alphaIncrement,
        proj = proj,
        clipToCoast = clipToCoast
      )
      sf::st_as_sf(made.polygons[[1]])
    }, error = function(err) return(err))
    iter=0
    lowerAlpha = max(initialAlpha-alphaDecrement,0)
    while(inherits(made.polygons,"try-error") && lowerAlpha!=0 && iter < maxIter){
      # Try making polygons with lower initial alpha value and other buffer radius
      made.polygons <- tryCatch({
        made.polygons <-rangeBuilder::getDynamicAlphaHull(
          x,
          coordHeaders = coordHeaders,
          fraction = fraction,
          partCount = partCount,
          buff = tryCatch(other_buffers[iter+1], error=function(err) return(other_buffers[length(other_buffers)])),
          initialAlpha = lowerAlpha,
          alphaIncrement = 0.5,
          proj = proj,
          clipToCoast = clipToCoast
        )
        sf::st_as_sf(made.polygons[[1]])
      }, error = function(err) { # If this fails return NULL
        return(NULL)
      })
      lowerAlpha = max(lowerAlpha-alphaDecrement,0)
      iter = iter + 1
    }
    if (is.null(made.polygons)) { # If we have returned NULL try making buffered circles
      made.polygons <- tryCatch({
        circ		<- dismo::circles(p = ll, d = other_buffers[length(other_buffers)])
        polygons	<- circ@polygons
        sf::st_as_sf(polygons)
      },error = function(err){ # If this fails then return NA
        return(list(NA))
      })
    }
  }else{ # If we have fewer points try making buffered circles
    made.polygons <- tryCatch({
      circ		<- dismo::circles(p=ll, d=other_buffers[length(other_buffers)])
      polygons	<- circ@polygons
      sf::st_as_sf(polygons)
    },
    error=function(err){ # If this fails return NA
      return(list(NA))
    })
  }

  suppressWarnings(
    if(inherits(made.polygons,"sf") && do.clipping && clipToCoast=="no"){ # If we don't have NA then clip polygon to coastline (suppress warnings)
      # Reproject to same projection as coastline
      made.polygons = sf::st_transform(made.polygons, sf::st_crs(proj))
      # Intersect
      made.polygons = sf::st_intersection(made.polygons, sf::st_geometry(coastline))

    }
  )
  return(made.polygons)
}

