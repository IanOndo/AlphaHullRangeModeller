#' Convert an arc into line segments given the center of the arc, the radius, the vector, and the angle (radians)
#' 
#' Modified from https://github.com/babichmorrowc/hull2spatial/blob/master/R/hull2line.R
#' See also https://babichmorrowc.github.io/post/2019-03-18-alpha-hull/
#'
#' @param center The coordinates of the center
#' @param r The radius of the arc
#' @param vector The vector of the arc
#' @param theta The angle of the arc in radians
#' @param npoints Number of points along the arc
#' @return a sf LINESTRING object that approximates the given arc
#' @export
#' @examples
#' arc2line(center = c(0.5, -0.1), r = 0.2, vector = c(0.1, 1), theta = 0.04)

arc2line <- function(center, r, vector, theta, npoints = 100) {
  # Get the angles at the extremes of the arcs
  angles <- alphahull::anglesArc(vector, theta)
  # Generate sequence of angles along the arc to determine the points
  seqang <- seq(angles[1], angles[2], length = npoints)
  # Generate x coordinates for points along the arc
  x <- center[1] + r * cos(seqang)
  # Generate y coordinates for points along the arc
  y <- center[2] + r * sin(seqang)
  coords.xy <- cbind(x,y)
  line <- sf::st_linestring(x = coords.xy, dim="XY")
  return(line)
}

#' Convert alpha hulls into sf MULTILINESTRING object
#'
#' @param hull Alpha hull
#' @return a sf MULTILINESTRING object
#' @export
#' @examples
#' library(alphahull)
#' set.seed(123)
#' x <- matrix(runif(100), nc = 2)
#' ahull_02 <- ahull(x, alpha = 0.2)
#' ahull_line_02 <- ahull2lines(ahull_02)

# Function to convert alpha hulls into MULTILINESTRING object
ahull2lines <- function(hull){
  if (!inherits(x, "ahull")) {
    stop("x needs to be an ahull class object")
  }
  arclist <- hull$arcs
  lines <- list()
  for (i in 1:nrow(arclist)) {
    # Extract the attributes of arc i
    center_i <- arclist[i, 1:2]
    radius_i <- arclist[i, 3]
    vector_i <- arclist[i, 4:5]
    theta_i <- arclist[i, 6]
    # Convert arc i into a LINESTRING object
    line_i <- arc2line(center = center_i, r = radius_i, vector = vector_i, theta = theta_i)
    list_length <- length(lines)
    if(list_length > 0){
      # If a line has already been added to the list of lines
      # Define last_line_coords as the coordinates of the last line added to the list before the ith line
      last_line_coords <- sf::st_coordinates(lines[[list_length]])[,1:2] # remove L1 id column if present
    }
    if(i == 1){
      # Add the first line to the list of lines
      lines[[i]] <- line_i
    } else if(isTRUE(all.equal(sf::st_coordinates(line_i)[1,], last_line_coords[nrow(last_line_coords),]))){
      # If the first coordinate in the ith line is equal to the last coordinate in the previous line
      # then those lines should be connected
      # Row bind the coordinates for the ith line to the coordinates of the previous line in the list
      coords = rbind(last_line_coords,
                     sf::st_coordinates(line_i)[2:nrow(sf::st_coordinates(line_i)),])
      lines[[list_length]]<- sf::st_linestring()
    } else {
      # If the first coordinate in the ith line does not match the last coordinate in the previous line
      # then the ith line represents a new line
      # Add the ith line to the list as a new element
      lines[[length(lines) + 1]] <- line_i
    }
  }

  # Convert the list of lines to a MULTILINESTRING object
  sf_lines <- sf::st_sfc(lines) %>% sf::st_cast("MULTILINESTRING")
  return(sf_lines)
}


#' Convert alpha hulls into sf POLYGON object
#'
#' @param hull Alpha hull
#' @return a sf POLYGON object
#' @export
#' @examples
#' library(alphahull)
#' set.seed(123)
#' x <- matrix(runif(100), nc = 2)
#' ahull_02 <- ahull(x, alpha = 0.2)
#' ahull_line_02 <- ah2sf(ahull_02)

# Function to convert alpha hulls into POLYGON object
ah2sf <- function(hull){
  if (!inherits(x, "ahull")) {
    stop("x needs to be an ahull class object")
  }
  # convert to MULTILINESTRING first
  sf_multilines <- ahull2lines(hull)
  
  # convert to POLYGON
  sf_poly <- sf_multilines[!sf::st_is_empty(sf_multilines)] %>% 
    sf::cast("POLYGON")
  
  return(sf_poly)
}

