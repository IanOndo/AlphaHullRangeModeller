% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hull2sf.R
\name{arc2line}
\alias{arc2line}
\title{Convert an arc into line segments given the center of the arc, the radius, the vector, and the angle (radians)}
\usage{
arc2line(center, r, vector, theta, npoints = 100)
}
\arguments{
\item{center}{The coordinates of the center}

\item{r}{The radius of the arc}

\item{vector}{The vector of the arc}

\item{theta}{The angle of the arc in radians}

\item{npoints}{Number of points along the arc}
}
\value{
a sf LINESTRING object that approximates the given arc
}
\description{
Modified from https://github.com/babichmorrowc/hull2spatial/blob/master/R/hull2line.R
See also https://babichmorrowc.github.io/post/2019-03-18-alpha-hull/
}
\examples{
arc2line(center = c(0.5, -0.1), r = 0.2, vector = c(0.1, 1), theta = 0.04)
}
