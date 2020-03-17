#' Calculate the "1/10th max" buffer method
#'
#' Computes the tenth of the maximum inter-point distance
#'
#' @param x A two-column data.frame of occurrence records coordinates of the species.
#' @param default_buffer A numeric specifying the buffering distance in meters around each point in case the buffer method equals 0. Default is 2000 meters.
#' @export
get_OneTenth_distmax <- function(x, default_buffer=2000, crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")){

  if(!inherits(x,"data.frame"))
    stop("Argument x must be a data.frame")
  if(default_buffer<=0)
    stop("Argument 'default_buffer' must be > 0")
  #---------------------------------------------------------------------------
  #= 1. get coords (with a regular expression to look for longitude, latitude)
  #---------------------------------------------------------------------------
  x_lon 	<- grep(pattern = "[Ll][Oo][Nn]",x = names(x),value=TRUE)[1]
  y_lat 	<- grep(pattern = "[Ll][Aa][Tt]",x = names(x),value=TRUE)[1]
  longLatNames <- c(x_lon,y_lat)
  ll 		<- sp::SpatialPoints(x[, longLatNames], crs)

  #-----------------------------------------------------
  #= 2. compute the 1/10th maximum inter-point distance
  #-----------------------------------------------------
  dMat 		<- tryCatch(geosphere::distGeo(ll), error=function(err) return(geosphere::distGeo(ll@coords[1,],ll@coords[2,])))
  if(length(dMat)==0L) return(default_buffer)
  dMax 		<- max(dMat, na.rm = TRUE)[1]
  onetenth	<- dMax / 10.
  if(onetenth==0){

    onetenth=default_buffer
  }
  return(onetenth)
}
