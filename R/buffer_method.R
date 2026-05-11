#' @title Compute one-tenth of the maximum inter-point distance
#'
#' @description Computes a buffer distance equal to one tenth of the maximum geographic
#' distance between occurrence points. This value is primarily intended for
#' adaptive alpha-hull buffering and range polygon construction.
#'
#' The “one tenth maximum distance” heuristic follows the approach proposed by
#' Rivers et al. (2010) for estimating species range parameters from occurrence data
#' in conservation assessments.
#'
#' For very large spatial extents, occurrence points can optionally be split
#' into spatial clusters before distance estimation in order to reduce the
#' influence of extreme outliers or disjunct distributions.
#'
#' @param x A \code{data.frame} containing longitude and latitude coordinates.
#'   Coordinate columns are automatically detected using column names matching
#'   patterns such as \code{lon}, \code{lat}, \code{x}, or \code{y}.
#' @param default_buffer Numeric. Default buffer distance (in metres) returned
#'   when distances cannot be computed or when all points are identical.
#' @param maxDist Numeric. Maximum allowed inter-point distance (in metres)
#'   before triggering spatial clustering.
#' @param k Integer. Number of clusters used when splitting spatially distant
#'   occurrence points.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Detects longitude and latitude columns
#'   \item Computes pairwise geographic distances using
#'     \code{geosphere::distGeo()}
#'   \item Extracts the maximum inter-point distance
#'   \item Returns one tenth of this maximum distance
#' }
#'
#' If the maximum distance exceeds \code{maxDist}, occurrence points are first
#' partitioned into \code{k} spatial clusters using equal-size clustering
#' (\code{kmeansEqual()}). Distances are then computed independently within
#' clusters, and the minimum cluster-specific maximum distance is retained.
#'
#' This strategy helps stabilise buffer estimation for species with highly
#' fragmented or geographically disjunct distributions.
#'
#' @return Numeric value corresponding to one tenth of the maximum inter-point
#' geographic distance (in metres).
#'
#' @references
#' Rivers, M. C., Bachman, S. P., Meagher, T. R., Nic Lughadha, E. M.,
#' & Brummitt, N. A. (2010). Subpopulations, locations and fragmentation:
#' Applying IUCN Red List criteria to herbarium specimen data.
#' \emph{Biodiversity and Conservation}, 19, 2071--2085.
#'
#' @importFrom geosphere distGeo
#'
#' @examples
#' \dontrun{
#' occ <- data.frame(
#'   lon = c(10, 12, 15),
#'   lat = c(45, 46, 48)
#' )
#'
#' get_OneTenth_distmax(occ)
#'
#' # Using custom defaults
#' get_OneTenth_distmax(
#'   occ,
#'   default_buffer = 100000,
#'   maxDist = 5000000
#' )
#' }
#'
#' @export
get_OneTenth_distmax <- function(x, default_buffer=200000, maxDist=10000000, k=3){#, crs=sf::st_crs("+proj=longlat +datum=WGS84")

  if(!inherits(x,"data.frame"))
    stop("Argument x must be a data.frame")
  if(default_buffer<=0)
    stop("Argument 'default_buffer' must be > 0")
  #---------------------------------------------------------------------------
  #= 1. get coords (with a regular expression to look for longitude, latitude)
  #---------------------------------------------------------------------------
  x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|^[Xx]$",x = names(x),value=TRUE)[1]
  y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|^[Yy]$",x = names(x),value=TRUE)[1]
  longLatNames <- c(x_lon,y_lat)
  ll 		<- x[, longLatNames] #sp::SpatialPoints(x[, longLatNames], crs)
  #sf_points <- sf::st_as_sf(ll,coords=c(1,2),crs=4326)
  #-----------------------------------------------------
  #= 2. compute the 1/10th maximum inter-point distance
  #-----------------------------------------------------
  dMat 		<- tryCatch(geosphere::distGeo(ll), error=function(err) return(geosphere::distGeo(ll[,1],ll[,2])))
  # dMat <- sf::st_distance(sf::st_transform(sf_points,crs="+proj=eqearth")) %>%
  #   units::drop_units()
  if(length(dMat)==0L) return(default_buffer)
  dMax 		<- max(dMat, na.rm = TRUE)[1]
  if(dMax>maxDist){
    llsplit <- kmeansEqual(ll,k=k, verbose=FALSE, plot=FALSE, iter_max=20)
    if(llsplit$converged!=TRUE) warning("The clustering algorithm did not converged !")
    dMat <- lapply(split(llsplit$Data, llsplit$Data$assigned), 
                   function(dt) tryCatch(geosphere::distGeo(dt[,1:2]), error=function(err) return(geosphere::distGeo(dt[,1],dt[,2]))))
    dMat 		<- dMat[lengths(dMat)>0L]
    if(length(dMat)==0L) return(default_buffer)
    dMax 		<- min(sapply(dMat,max,na.rm = TRUE))
  }
  onetenth	<- dMax / 10.#units::drop_units(dMax / 10.)
  if(onetenth==0){
    onetenth=default_buffer
  }
  return(onetenth)
}

