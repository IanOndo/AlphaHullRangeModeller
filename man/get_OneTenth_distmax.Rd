% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buffer_method.R
\name{get_OneTenth_distmax}
\alias{get_OneTenth_distmax}
\title{Calculate the "1/10th max" buffer method}
\usage{
get_OneTenth_distmax(x, default_buffer = 2e+05, maxDist = 1e+07, k = 3)
}
\arguments{
\item{x}{A two-column data.frame of occurrence records coordinates of the species.}

\item{default_buffer}{A numeric specifying the buffering distance in meters around each point in case the buffer method equals 0. Default is 200000 meters.}

\item{maxDist}{A numeric value defining the maximum interpoint distance before geographically splitting the occurrence records in k groups. Defaults is 10,000km.}

\item{k}{number of groups used to split the records. Default is 3.}
}
\description{
Computes the tenth of the maximum inter-point distance
}
