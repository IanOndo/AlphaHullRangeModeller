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
#' @param maxIter A numeric integer specifying the maximum number of iterations before trying to build circular buffers. Default is 2.
#' @param other_buffers A numeric value or a vector of values of length \code{maxIter} specifying alternative buffer sizes (in meters) if the function fails after the first iterations.
#' @export
PolygonMaker <- function(x, coordHeaders = NULL, fraction=0.95, partCount=10, buffer, initialAlpha=2, alphaIncrement=1, alphaDecrement=1, maxIter=2, other_buffers=rep(200000, maxIter), clipToCoast="terrestrial", coastline=NULL, proj='+proj=longlat +datum=WGS84'){

  if(!inherits(x,"data.frame"))
    stop("Argument x must be a data.frame")
  if(buffer<=0)
    stop("Argument 'default_buffer' must be > 0")
  do.clipping = !is.null(coastline) && inherits(coastline, c("sf","sfc")) && clipToCoast!="no"
  #---------------------------------------------------------------------------
  #= 1. get coords (with a regular expression to look for longitude, latitude)
  #---------------------------------------------------------------------------
  if(is.null(coordHeaders)|length(coordHeaders)<2){
    x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|^[Xx]$",x = names(x),value=TRUE)[1]
    y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|^[Yy]$",x = names(x),value=TRUE)[1]
    coordHeaders <- c(x_lon,y_lat)
  }
  x  <- x[, coordHeaders]
  ll <- sf::st_as_sf(x, coords=c(1,2), crs=sf::st_crs(proj))
  
  # transform to lon/lat for getDynamicAlphaHull to run
  is_lonlat <- sf::st_is_longlat(ll)
  if(!is_lonlat){
    x <- sf::st_transform(ll, crs=4326) %>%
      sf::st_coordinates() %>%
      as.data.frame()
  }

  # Do we have at least 3 points?
  if (length(sf::st_geometry(ll)) >= 3){
    # Try making polygons
    made.polygons <- tryCatch(AlphaHullRangeModeller::getDynamicAlphaHull(
          x,
          coordHeaders = coordHeaders,
          fraction = fraction,
          partCount = partCount,
          buff = buffer,
          initialAlpha = initialAlpha,
          alphaIncrement = alphaIncrement,
          clipToCoast = "no"
        ),
        error=function(err){
          # if we have 3 collinear points
          if(length(sf::st_geometry(ll)) == 3){
            iter = 0
            while(iter<maxIter){
              # sample non-collinear new coordinates within the buffer #dismo::circles
              # circ		<- lapply(c(nrow(x)-1,nrow(x)), function(row){
              #   sf::st_buffer(
              #     sf::st_transform(sf::st_as_sf(x[row,],
              #                                   coords=c(1,2),
              #                                   crs=4326),
              #                      crs="+proj=eqearth"),
              #     dist = buffer)}
              #   )
              # list.coords <- lapply(circ, function(buf){
              #   as.vector(
              #     sf::st_coordinates(sf::st_sample(buf, size=1, type="random"))
              #     )
              #   })
              # x.new   <- rbind(x[-c(nrow(x)-1,nrow(x)),], setNames(data.frame(do.call(rbind, list.coords)), colnames(x)))
              x.new <- sf::jitter(x) %>%
                sf::st_coordinates() %>%
                as.data.frame()
              
              made.polygons <- try(AlphaHullRangeModeller::getDynamicAlphaHull(
                x.new,
                coordHeaders = c("x","y"),#coordHeaders,
                fraction = fraction,
                partCount = partCount,
                buff = buffer,
                initialAlpha = initialAlpha,
                alphaIncrement = alphaIncrement,
                clipToCoast = "no"
              ),silent=TRUE)
              
              if(!inherits(made.polygons, "try-error")) break
              iter = iter + 1
            }
            return(made.polygons)
          }else{
            return(err)
          }
        }
      )

    # If an error is returned try changing parameters 
    iter=0
    lowerAlpha = max(initialAlpha-alphaDecrement,0)
    while(inherits(made.polygons,c("error","try-error")) && lowerAlpha!=0 && iter < maxIter){
      # Try making polygons with lower initial alpha value and other buffer radius
      made.polygons <- tryCatch({
        made.polygons <-AlphaHullRangeModeller::getDynamicAlphaHull(
          x,
          coordHeaders = coordHeaders,
          fraction = fraction,
          partCount = partCount,
          buff = tryCatch(other_buffers[iter+1], error=function(err) return(tail(other_buffers,1))),
          initialAlpha = lowerAlpha,
          alphaIncrement = 0.5,
          clipToCoast = "no"
        )
        sf::st_buffer(made.polygons[[1]],dist=0)
      }, 
      error = function(err) return(err))# If this fails return error
      
      lowerAlpha = max(lowerAlpha-alphaDecrement,0)
      iter = iter + 1
    }

    # If another error is returned try making buffered points
    if(inherits(made.polygons, c("error","try-error"))) {
      print(made.polygons[[1]])
      made.polygons <- tryCatch({
        sf::st_buffer(sf::st_transform(ll,"+proj=eqearth +wktext"),
                      dist = tail(other_buffers,1)) %>%
          sf::st_union() %>%
          sf::st_transform(4326) %>% list()
      },
      error = function(err){ # If this fails then return NA
        return(list(NA))
      })
    }
    
  }
  else{ # If we have fewer points try making buffered points
    made.polygons <- tryCatch({
      sf::st_buffer(sf::st_transform(ll,"+proj=eqearth +wktext"),
                    dist = tail(other_buffers,1)) %>%
        sf::st_union() %>%
        sf::st_transform(4326) %>% list()
    },
    error=function(err){# If this fails return NA
      return(list(NA))
    })
  }

  if(inherits(made.polygons[[1]],c("sf","sfc")) && any(!sf::st_is_valid(made.polygons[[1]])))
    made.polygons[[1]] <- sf::st_make_valid(made.polygons[[1]])  

  # If we don't have NA then clip polygon to coastline (suppress warnings)
  if(inherits(made.polygons[[1]],c("sf","sfc")) && do.clipping){ 
    # Reproject to same projection as coastline
    made.polygons = suppressWarnings(
      suppressMessages(
        sf::st_transform(made.polygons[[1]],
                         crs=sf::st_crs(coastline))
        )
      )
    # Intersect
    made.polygons = tryCatch({
      sf_fun <- switch(clipToCoast, terrestrial=sf::st_intersection, aquatic=sf::st_difference)
      suppressWarnings(
        suppressMessages(
          sf_fun(made.polygons, sf::st_geometry(coastline))
          )
        )
      },error=function(err){
        stop(err)
      }
    )
    # transform back to original projection if necessary
    if(!is_lonlat) return(sf::st_transform(made.polygons, crs=proj))
    
    return(made.polygons)
  }

  # transform back to original projection if necessary
  if(!is_lonlat) return(sf::st_transform(made.polygons[[1]], crs=proj))
  
  return(made.polygons[[1]])
}

