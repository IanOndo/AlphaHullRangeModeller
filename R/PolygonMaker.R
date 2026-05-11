#' @title Generate species range polygons using adaptive alpha-hull methods
#'
#' @description Constructs species range polygons from geographic occurrence coordinates
#' using adaptive alpha-hull algorithms. The functions iteratively adjust alpha
#' parameters to generate polygons satisfying user-defined completeness and
#' fragmentation constraints, with multiple fallback strategies when alpha-hull
#' construction fails.
#'
#' Three implementations are provided:
#' \itemize{
#'   \item \code{PolygonMaker}: standard R implementation
#'   \item \code{PolygonMakerRcpp}: Rcpp-accelerated implementation
#'   \item \code{FastPolygonMaker}: optimized implementation with accelerated
#'     coastline clipping using spatial tiling
#' }
#'
#' If alpha-hull construction fails, the functions progressively:
#' \enumerate{
#'   \item decrease alpha values,
#'   \item modify buffer sizes,
#'   \item and finally fall back to buffered occurrence polygons.
#' }
#'
#' Special handling is implemented for collinear points and small occurrence
#' datasets.
#'
#' @name PolygonMaker
#'
#' @param x A two-column \code{data.frame} containing occurrence coordinates.
#' @param coordHeaders Optional character vector of length 2 specifying the
#'   longitude and latitude column names.
#' @param fraction Numeric value between 0 and 1 specifying the minimum
#'   fraction of occurrence points that must be enclosed within the polygon.
#' @param partCount Integer specifying the maximum number of disjunct polygon
#'   components allowed.
#' @param buffer Numeric buffering distance (in metres) applied around points.
#' @param initialAlpha Numeric initial alpha value used during alpha-hull
#'   construction.
#' @param alphaIncrement Numeric increment applied to alpha during iterative
#'   optimisation.
#' @param alphaDecrement Numeric decrement applied to alpha when polygon
#'   generation fails.
#' @param maxIter Integer specifying the maximum number of retry attempts.
#' @param other_buffers Numeric vector of alternative buffer distances (in
#'   metres) used during fallback attempts.
#' @param clipToCoast Character string specifying coastline clipping mode:
#'   \itemize{
#'     \item \code{"terrestrial"}: retain terrestrial portion only
#'     \item \code{"aquatic"}: retain aquatic portion only
#'     \item \code{"no"}: disable coastline clipping
#'   }
#' @param coastline Optional coastline polygons provided as \code{sf} or
#'   \code{sfc} objects.
#' @param proj Character string specifying the coordinate reference system.
#' @param tiles Optional \code{sf} tiling grid used by
#'   \code{FastPolygonMaker()} to accelerate coastline intersection operations.
#'
#' @details
#' These functions are designed for automated species range reconstruction from
#' occurrence records, particularly for biodiversity and conservation
#' applications.
#'
#' Polygon generation relies on iterative adaptive alpha-hull construction
#' implemented in:
#' \itemize{
#'   \item \code{getDynamicAlphaHull()}
#'   \item \code{getDynamicRcppAlphaHull()}
#'   \item \code{getDynamicFastAlphaHull()}
#' }
#'
#' The workflow:
#' \enumerate{
#'   \item Detects coordinate columns
#'   \item Converts coordinates to \code{sf} geometries
#'   \item Attempts alpha-hull construction
#'   \item Iteratively adjusts alpha and buffer parameters if needed
#'   \item Falls back to buffered occurrence polygons when necessary
#'   \item Optionally clips polygons to coastlines
#'   \item Returns valid \code{sf} geometries
#' }
#'
#' The \code{FastPolygonMaker()} implementation improves clipping performance
#' by spatially partitioning coastlines into tiles prior to intersection,
#' substantially reducing geometric filtering costs for large global datasets.
#'
#' @return An \code{sf} polygon geometry representing the estimated species
#'   range. Returns \code{NA} if polygon construction fails.
#'
#' @references
#' Rivers, M. C., Bachman, S. P., Meagher, T. R., Nic Lughadha, E. M.,
#' & Brummitt, N. A. (2010). Subpopulations, locations and fragmentation:
#' Applying IUCN Red List criteria to herbarium specimen data.
#' \emph{Biodiversity and Conservation}, 19, 2071--2085.
#'
#' @importFrom sf st_as_sf st_transform st_geometry st_buffer st_union
#'   st_make_valid st_is_valid st_filter st_intersection st_difference
#'   st_combine st_make_grid st_set_crs st_coordinates st_crs st_is_longlat
#' @importFrom dplyr rename mutate row_number filter
#'
#' @examples
#' \dontrun{
#' # Standard implementation
#' poly <- PolygonMaker(
#'   x = occ_data,
#'   buffer = 200000
#' )
#'
#' # Rcpp implementation
#' poly <- PolygonMakerRcpp(
#'   x = occ_data,
#'   buffer = 200000
#' )
#'
#' # Fast implementation with coastline clipping
#' poly <- FastPolygonMaker(
#'   x = occ_data,
#'   buffer = 200000,
#'   clipToCoast = "terrestrial",
#'   coastline = world_coast
#' )
#' }

#' @rdname PolygonMaker
#' @export
PolygonMaker <- function(x,
                         coordHeaders = NULL,
                         fraction=0.95,
                         partCount=10,
                         buffer,
                         initialAlpha=2,
                         alphaIncrement=1,
                         alphaDecrement=1,
                         maxIter=2,
                         other_buffers=rep(200000, maxIter),
                         clipToCoast="terrestrial",
                         coastline=NULL,
                         proj='+proj=longlat +datum=WGS84'){

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
    made.polygons <- tryCatch(getDynamicAlphaHull(
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
              
              made.polygons <- try(getDynamicAlphaHull(
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
        made.polygons <-getDynamicAlphaHull(
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

  if(!inherits(made.polygons,"list")){
    made.polygons<- list(made.polygons)
  }
  
  if(inherits(made.polygons[[1]],c("sf","sfc","sfg")) && any(!sf::st_is_valid(made.polygons[[1]])))
    made.polygons[[1]] <- sf::st_make_valid(made.polygons[[1]])  

  # If we don't have NA then clip polygon to coastline (suppress warnings)
  if(inherits(made.polygons[[1]],c("sf","sfc","sfg")) && do.clipping){ 
    # Reproject to same projection as coastline
    made.polygons = suppressWarnings(
      suppressMessages(
        sf::st_transform(made.polygons[[1]],
                         crs=sf::st_crs(coastline))
        )
      )
    # Intersect
    made.polygons = tryCatch({
      sf_fun <- switch(clipToCoast, 
                       terrestrial=sf::st_intersection, 
                       aquatic=function(x, y) sf::st_difference(x, sf::st_union(sf::st_combine(y))))
      coast_intersect <- sf::st_filter(coastline,
                                       made.polygons,
                                       .predicate = sf::st_intersects)
      suppressWarnings(
        suppressMessages(
          sf_fun(made.polygons, sf::st_geometry(coast_intersect))
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
#' @rdname PolygonMaker
#' @export
PolygonMakerRcpp <- function(x,
                         coordHeaders = NULL,
                         fraction=0.95,
                         partCount=10,
                         buffer,
                         initialAlpha=2,
                         alphaIncrement=1,
                         alphaDecrement=1,
                         maxIter=2,
                         other_buffers=rep(200000, maxIter),
                         clipToCoast="terrestrial",
                         coastline=NULL,
                         proj='+proj=longlat +datum=WGS84'){
  
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
    made.polygons <- tryCatch(getDynamicRcppAlphaHull(
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
          x.new <- sf::jitter(x) %>%
            sf::st_coordinates() %>%
            as.data.frame()
          
          made.polygons <- try(getDynamicRcppAlphaHull(
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
        made.polygons <-getDynamicRcppAlphaHull(
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
  
  if(inherits(made.polygons[[1]],c("sf","sfc","sfg")) && any(!sf::st_is_valid(made.polygons[[1]])))
    made.polygons[[1]] <- sf::st_make_valid(made.polygons[[1]])  
  
  # If we don't have NA then clip polygon to coastline (suppress warnings)
  if(inherits(made.polygons[[1]],c("sf","sfc","sfg")) && do.clipping){ 
    # Reproject to same projection as coastline
    made.polygons = suppressWarnings(
      suppressMessages(
        sf::st_transform(made.polygons[[1]],
                         crs=sf::st_crs(coastline))
      )
    )
    # Intersect
    made.polygons = tryCatch({
      sf_fun <- switch(clipToCoast, 
                       terrestrial=sf::st_intersection, 
                       aquatic=function(x, y) sf::st_difference(x, sf::st_union(sf::st_combine(y))))
      coast_intersect <- sf::st_filter(coastline,
                                       made.polygons,
                                       .predicate = sf::st_intersects)
      suppressWarnings(
        suppressMessages(
          sf_fun(made.polygons, sf::st_geometry(coast_intersect))
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

#' @rdname PolygonMaker
#' @export
FastPolygonMaker <- function(x,
                         coordHeaders = NULL,
                         fraction=0.95,
                         partCount=10,
                         buffer,
                         initialAlpha=2,
                         alphaIncrement=1,
                         alphaDecrement=1,
                         maxIter=2,
                         other_buffers=rep(200000, maxIter),
                         clipToCoast="terrestrial",
                         coastline=NULL,
                         proj='+proj=longlat +datum=WGS84',
                         tiles =  NULL){
  
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
    made.polygons <- tryCatch(getDynamicFastAlphaHull(
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
          
          made.polygons <- try(getDynamicFastAlphaHull(
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
        made.polygons <-getDynamicFastAlphaHull(
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
  
  if(!inherits(made.polygons,"list")){
    made.polygons<- list(made.polygons)
  }
  
  if(inherits(made.polygons[[1]],c("sf","sfc","sfg")) && any(!sf::st_is_valid(made.polygons[[1]])))
    made.polygons[[1]] <- sf::st_make_valid(made.polygons[[1]])  
  
  # If we don't have NA then clip polygon to coastline (suppress warnings)
  if(inherits(made.polygons[[1]],c("sf","sfc","sfg")) && do.clipping){ 
    # Reproject to same projection as coastline
    made.polygons = suppressWarnings(
      suppressMessages(
        sf::st_transform(made.polygons[[1]],
                         crs=sf::st_crs(coastline))
      )
    )
    # Intersect
    made.polygons = tryCatch({
      sf_fun <- switch(clipToCoast, 
                       terrestrial=sf::st_intersection, 
                       aquatic=function(x, y) sf::st_difference(x, sf::st_union(sf::st_combine(y))))
      
      if(is.null(tiles)){
        tiles <- sf::st_make_grid(coastline, n = c(72, 36)) |>
          sf::st_as_sf() |>
          sf::st_set_crs(sf::st_crs(coastline)) |>
          dplyr::rename(geometry=x) |>
          dplyr::mutate(tile_ID = dplyr::row_number(), .before=geometry)
      }

      xt <- sf::st_join(coastline, 
                            tiles, 
                            join = sf::st_intersects,
                            left = TRUE)
      yt <- tiles[made.polygons, op=sf::st_intersects]
      
      coast_intersect <- xt |>
        dplyr::filter(tile_ID %in% yt$tile_ID) #|>
        #dplyr::select(1,geometry) |>
        #dplyr::distinct()
      
      # coast_intersect <- sf::st_filter(coastline,
      #                                  made.polygons,
      #                                  .predicate = sf::st_intersects)

      # coast_simp <- sf::st_simplify(coastline, dTolerance = dTolerance * 0.008333333333333)
      # idx <- lengths(sf::st_intersects(coast_simp, made.polygons)) > 0
      # coast_intersect <- coastline[idx, ]
      
      suppressWarnings(
        suppressMessages(
          sf_fun(made.polygons, sf::st_geometry(coast_intersect))
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