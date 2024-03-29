% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDynamicAlphaHull.R
\name{getDynamicAlphaHull}
\alias{getDynamicAlphaHull}
\title{Generate polygon based on alpha hulls}
\usage{
getDynamicAlphaHull(
  x,
  fraction = 0.95,
  partCount = 3,
  buff = 10000,
  initialAlpha = 3,
  coordHeaders = c("Longitude", "Latitude"),
  clipToCoast = "terrestrial",
  alphaIncrement = 1,
  verbose = FALSE,
  alphaCap = 400
)
}
\arguments{
\item{x}{dataframe of coordinates in decimal degrees, with a minimum of 3
rows.}

\item{fraction}{the minimum fraction of occurrences that must be included
in polygon.}

\item{partCount}{the maximum number of disjunct polygons that are allowed.}

\item{buff}{buffering distance in meters}

\item{initialAlpha}{the starting value for alpha}

\item{coordHeaders}{the column names for the longitude and latitude
columns, respectively.  If x has two columns, these are assumed to be
longitude and latitude, and \code{coordHeaders} is ignored.}

\item{clipToCoast}{Either "no" (no clipping), "terrestrial" (only
terrestrial part of range is kept) or "aquatic" (only non-terrestrial part
is clipped). See Details.}

\item{alphaIncrement}{the amount to increase alpha with each iteration}

\item{verbose}{prints the alpha value to the console, intended for
debugging.}

\item{alphaCap}{Max alpha value before function aborts and returns a
minimum convex hull.}
}
\value{
a list with 2 elements: \item{hull}{ a sf polygon object }
\item{alpha}{ the alpha value that was found to satisfy the criteria.  If a
convex hull was returned, this will list MCH.  }
}
\description{
Generates an apha hull polygon, where the alpha parameter is determined by
the spatial distribution of the coordinates.
}
\details{
From a set of coordinates, this function will create an alpha hull with
\code{alpha = initialAlpha}, and will then increase \code{alpha} by
\code{alphaIncrement} until both the \code{fraction} and \code{partCount}
conditions are met.

If the conditions cannot be satisfied, then a minimum convex hull is
returned.

If \code{clipToCoast} is set to "terrestrial" or "aquatic", the resulting
polygon is clipped to the coastline, using a basemap from naturalearth. 
The first time this function is run, this basemap will be downloaded. 
Subsequent calls will use the downloaded map.
}
\examples{

data(crotalus)

# create a polygon range for Crotalus atrox
x <- crotalus[which(crotalus$genSp == 'Crotalus_atrox'),]
x <- x[sample(1:nrow(x), 50),]

range <- getDynamicAlphaHull(x, coordHeaders=c('decimallongitude','decimallatitude'), 
	clipToCoast = 'no')

plot(range[[1]], col=transparentColor('dark green', 0.5), border = NA)
points(x[,c('decimallongitude','decimallatitude')], cex = 0.5, pch = 3)

# to add a basic coastline, you can use the internal map
# world <- rangeBuilder:::loadWorldMap()
# plot(world, add = TRUE, lwd = 0.5)

}
\seealso{
Alpha hulls are created with \code{\link{ahull}}.
}
\author{
Pascal Title
}
