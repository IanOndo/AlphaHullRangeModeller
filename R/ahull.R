#' Fast alpha hull computation using C++ triangulation
#'
#' Computes the alpha hull of a planar point set using a fast C++ backend for
#' the Delaunay–Voronoi mesh construction. This function behaves like
#' \code{\link[alphahull]{ahull}} but replaces the internal call to
#' \code{\link[alphahull]{delvor}} with a faster implementation.
#'
#' The alpha hull is a generalization of the convex hull that captures the
#' shape of a point cloud depending on the parameter \code{alpha}. Smaller
#' values of \code{alpha} allow concave features to appear, while large values
#' approach the convex hull.
#'
#' @param x A vector of x coordinates, a two-column matrix of coordinates, or
#' an object coercible via \code{\link[grDevices]{xy.coords}}.
#' @param y Optional vector of y coordinates if \code{x} contains only x values.
#' @param alpha A positive numeric value controlling the level of detail of the
#' alpha hull.
#'
#' @return
#' An object of class \code{"ahull"} identical in structure to the output of
#' \code{\link[alphahull]{ahull}}. The object contains:
#' \describe{
#'   \item{arcs}{Matrix describing the circular arcs that form the boundary
#'   of the alpha hull.}
#'   \item{xahull}{Coordinates of points used to define arc endpoints.}
#' }
#'
#' @details
#' This function accelerates alpha hull computation by replacing the Delaunay–
#' Voronoi construction step with a C++ implementation. All subsequent steps,
#' including arc construction and intersection handling, rely on the original
#' algorithms implemented in \pkg{alphahull}.
#'
#' The function is therefore compatible with existing utilities such as
#' plotting methods and conversion routines.
#'
#' @seealso
#' \code{\link[alphahull]{ahull}}, \code{\link{delvor_fast}},
#' \code{\link[alphahull]{ashape}}
#'
#' @examples
#' pts <- matrix(runif(200), ncol = 2)
#'
#' ah <- ahull_fast(pts, alpha = 1)
#'
#' plot(ah)
#'
#' @export
ahull_fast <- function(x, y = NULL, alpha) {
  ashape.obj <- ashape_fast(x, y, alpha)
  compl <- complement_fast(ashape.obj$delvor.obj, alpha = alpha)
  
  pm.x <- (compl[, "x1"] + compl[, "x2"]) * 0.5
  pm.y <- (compl[, "y1"] + compl[, "y2"]) * 0.5
  dm <- sqrt((compl[, "x1"] - compl[, "x2"])^2 + (compl[, "y1"] - compl[, "y2"])^2) * 0.5
  
  ashape.edges <- matrix(ashape.obj$edges[, c("ind1", "ind2")], ncol = 2, byrow = FALSE)
  noforget <- ashape.obj$alpha.extremes
  ind2 <- integer()
  j <- 0
  nshape <- length(ashape.edges) * 0.5
  arcs <- matrix(0, nrow = nshape, ncol = 6)
  indp <- matrix(0, nrow = nshape, ncol = 2)
  cutp <- ashape.obj$x
  
  if (nshape > 0) {
    for (i in 1:nshape) {
      ind <- which(ashape.edges[i, 1] == compl[, "ind1"] &
                     ashape.edges[i, 2] == compl[, "ind2"])
      if (length(ind) > 0) {
        if (!((1 <= sum((compl[ind, "ind"] == 1))) &
              (sum((compl[ind, "ind"] == 1)) < length(ind)))) {
          which_min <- which(compl[ind, "r"] == min(compl[ind[compl[ind, "r"] > 0], "r"]))
          j <- j + 1
          arcs[j, ] <- c(compl[ind[which_min], 1], compl[ind[which_min], 2], compl[ind[which_min], 3],
                         compl[ind[which_min], "v.x"], compl[ind[which_min], "v.y"], compl[ind[which_min], "theta"])
          vaux <- compl[ind[which_min], c("x1", "y1")] - cbind(pm.x, pm.y)[ind[which_min], ]
          theta.aux <- alphahull::rotation(compl[ind[which_min], c("v.x", "v.y")], compl[ind[which_min], "theta"])
          a2 <- sum(vaux * theta.aux)
          if (a2 > 0) {
            indp[j, ] <- compl[ind[which_min], c("ind1", "ind2")]
          } else {
            indp[j, ] <- compl[ind[which_min], c("ind2", "ind1")]
          }
        }
        ind2 <- c(ind2, ind)
      }
    }
  }
  
  arcs.old <- arcs[arcs[, 3] > 0, , drop = FALSE]
  colnames(arcs.old) <- c("c1", "c2", "r", "v.x", "v.y", "theta")
  arcs <- arcs.old
  indp <- indp[indp[, 1] != 0 & indp[, 2] != 0, , drop = FALSE]
  n.arc <- dim(arcs)[1]
  
  watch <- 1
  j <- 1
  
  if (n.arc > 0) {
    while (watch <= n.arc) {
      while (j <= n.arc) {
        if (j != watch) {
          intersection <- alphahull::inter(arcs[watch, 1], arcs[watch, 2], arcs[watch, 3],
                                           arcs[j, 1], arcs[j, 2], arcs[j, 3])
          if (intersection$n.cut == 2) {
            v.arc <- c(arcs[watch, "v.x"], arcs[watch, "v.y"])
            if (v.arc[2] >= 0) {
              ang.OX <- acos(v.arc[1])
            } else {
              ang.OX <- 2 * pi - acos(v.arc[1])
            }
            v.int <- intersection$v1
            v.int.rot <- alphahull::rotation(v.int, ang.OX)
            
            if (v.int.rot[2] >= 0) {
              ang.v.int.rot.OX <- acos(v.int.rot[1])
              angles <- c(-arcs[watch, "theta"], arcs[watch, "theta"],
                          ang.v.int.rot.OX - intersection$theta1,
                          ang.v.int.rot.OX + intersection$theta1)
              names(angles) <- c("theta1", "theta2", "beta1", "beta2")
              order <- names(sort(angles))
            } else {
              ang.v.int.rot.OX <- acos(v.int.rot[1])
              angles <- c(-arcs[watch, "theta"], arcs[watch, "theta"],
                          -ang.v.int.rot.OX - intersection$theta1,
                          -ang.v.int.rot.OX + intersection$theta1)
              names(angles) <- c("theta1", "theta2", "beta1", "beta2")
              order <- names(sort(angles))
            }
            
            if (sum(match(indp[watch, ], indp[j, ], nomatch = 0)) > 0) {
              if (indp[watch, 1] == indp[j, 2]) {
                if (all(order == c("beta1", "beta2", "theta1", "theta2"))) {
                  case <- 1
                } else if (all(order == c("theta1", "beta1", "beta2", "theta2"))) {
                  case <- 2
                } else if (all(order == c("beta1", "theta1", "beta2", "theta2"))) {
                  ang.control <- (angles["theta1"] - angles["beta1"]) / 2
                  case <- if (abs(ang.control) < 1e-05) 2 else 1
                }
                
                if (case == 2) {
                  ang.middle2 <- (angles["theta2"] - angles["beta2"]) / 2
                  v.new2 <- alphahull::rotation(c(1, 0), -arcs[watch, "theta"] + ang.middle2 - ang.OX)
                  cutp <- rbind(cutp, arcs[watch, 1:2] + arcs[watch, 3] * alphahull::rotation(v.new2, ang.middle2))
                  inn <- dim(cutp)[1]
                  arcs[watch, 4:6] <- c(v.new2, ang.middle2)
                  indp[watch, 1] <- inn
                  pmaux <- (cutp[inn, ] + cutp[indp[j, 1], ]) * 0.5
                  dmaux <- pmaux - cutp[indp[j, 1], ]
                  ndmaux <- sqrt(sum(dmaux^2))
                  vaux <- pmaux - arcs[j, 1:2]
                  nvaux <- sqrt(sum(vaux^2))
                  th <- atan(ndmaux / nvaux)
                  arcs[j, 4:6] <- c(vaux / nvaux, th)
                  indp[j, 2] <- inn
                }
              } else if (indp[watch, 2] == indp[j, 1]) {
                if (all(order == c("theta1", "theta2", "beta1", "beta2"))) {
                  case <- 1
                } else if (all(order == c("theta1", "beta1", "beta2", "theta2"))) {
                  case <- 2
                } else if (all(order == c("theta1", "beta1", "theta2", "beta2"))) {
                  ang.control <- (angles["theta2"] - angles["beta2"]) / 2
                  case <- if (abs(ang.control) < 1e-05) 2 else 1
                }
                
                if (case == 2) {
                  ang.middle <- (angles["beta1"] - angles["theta1"]) / 2
                  v.new <- alphahull::rotation(c(1, 0), arcs[watch, "theta"] - ang.middle - ang.OX)
                  cutp <- rbind(cutp, arcs[watch, 1:2] + arcs[watch, 3] * alphahull::rotation(v.new, -ang.middle))
                  inn <- dim(cutp)[1]
                  arcs[watch, 4:6] <- c(v.new, ang.middle)
                  indp[watch, 2] <- inn
                  pmaux <- (cutp[inn, ] + cutp[indp[j, 2], ]) * 0.5
                  dmaux <- pmaux - cutp[indp[j, 2], ]
                  ndmaux <- sqrt(sum(dmaux^2))
                  vaux <- pmaux - arcs[j, 1:2]
                  nvaux <- sqrt(sum(vaux^2))
                  th <- atan(ndmaux / nvaux)
                  arcs[j, 4:6] <- c(vaux / nvaux, th)
                  indp[j, 1] <- inn
                }
              }
            } else if (all(order == c("theta1", "beta1", "beta2", "theta2"))) {
              v.arcj <- c(arcs[j, "v.x"], arcs[j, "v.y"])
              if (v.arcj[2] >= 0) {
                ang.OXj <- acos(v.arcj[1])
              } else {
                ang.OXj <- 2 * pi - acos(v.arcj[1])
              }
              
              v.intj <- intersection$v2
              v.int.rotj <- alphahull::rotation(v.intj, ang.OXj)
              
              if (v.int.rotj[2] >= 0) {
                ang.v.int.rot.OXj <- acos(v.int.rotj[1])
                anglesj <- c(-arcs[j, "theta"], arcs[j, "theta"],
                             ang.v.int.rot.OXj - intersection$theta2,
                             ang.v.int.rot.OXj + intersection$theta2)
                names(anglesj) <- c("theta1", "theta2", "beta1", "beta2")
                orderj <- names(sort(anglesj))
              } else {
                ang.v.int.rot.OXj <- acos(v.int.rotj[1])
                anglesj <- c(-arcs[j, "theta"], arcs[j, "theta"],
                             -ang.v.int.rot.OXj - intersection$theta2,
                             -ang.v.int.rot.OXj + intersection$theta2)
                names(anglesj) <- c("theta1", "theta2", "beta1", "beta2")
                orderj <- names(sort(anglesj))
              }
              
              if (all(orderj == c("theta1", "beta1", "beta2", "theta2"))) {
                ang.middle <- (angles["beta1"] - angles["theta1"]) / 2
                v.new <- alphahull::rotation(c(1, 0), arcs[watch, "theta"] - ang.middle - ang.OX)
                ang.middle2 <- (angles["theta2"] - angles["beta2"]) / 2
                v.new2 <- alphahull::rotation(c(1, 0), -arcs[watch, "theta"] + ang.middle2 - ang.OX)
                
                arcs <- rbind(arcs, c(arcs[watch, 1], arcs[watch, 2], arcs[watch, 3], v.new2[1], v.new2[2], ang.middle2))
                arcs[watch, ] <- c(arcs[watch, 1], arcs[watch, 2], arcs[watch, 3], v.new[1], v.new[2], ang.middle)
                n.arc <- n.arc + 1
                
                np1 <- arcs[watch, 1:2] + arcs[watch, 3] * alphahull::rotation(v.new, -ang.middle)
                np2 <- arcs[watch, 1:2] + arcs[watch, 3] * alphahull::rotation(v.new2, ang.middle2)
                
                indold <- indp[watch, 2]
                inn1 <- dim(cutp)[1] + 1
                inn2 <- dim(cutp)[1] + 2
                indp[watch, 2] <- inn1
                indp <- rbind(indp, c(inn2, indold))
                cutp <- rbind(cutp, np1, np2)
                
                indold <- indp[j, 1]
                indp[j, 1] <- inn1
                indp <- rbind(indp, c(indold, inn2))
                
                pmaux <- (cutp[inn1, ] + cutp[indp[j, 2], ]) * 0.5
                dmaux <- pmaux - cutp[indp[j, 2], ]
                ndmaux <- sqrt(sum(dmaux^2))
                vaux <- pmaux - arcs[j, 1:2]
                nvaux <- sqrt(sum(vaux^2))
                th <- atan(ndmaux / nvaux)
                arcs[j, 4:6] <- c(vaux / nvaux, th)
                
                pmaux <- (cutp[inn2, ] + cutp[indold, ]) * 0.5
                dmaux <- pmaux - cutp[indold, ]
                ndmaux <- sqrt(sum(dmaux^2))
                vaux <- pmaux - arcs[j, 1:2]
                nvaux <- sqrt(sum(vaux^2))
                th <- atan(ndmaux / nvaux)
                arcs <- rbind(arcs, c(arcs[j, 1:3], vaux / nvaux, th))
                n.arc <- n.arc + 1
              }
            }
          }
          case <- 0
        }
        j <- j + 1
      }
      watch <- watch + 1
      j <- 1
    }
    
    ord.old <- 1:dim(indp)[1]
    ord.new <- numeric()
    
    while (length(ord.new) < length(ord.old)) {
      if (length(ord.new) == 0) {
        ord.new <- 1
      } else {
        ord.new <- c(ord.new, ord.old[-ord.new][1])
      }
      coinc <- match(indp[ord.new[length(ord.new)], 2], indp[-ord.new, 1])
      while (!is.na(coinc)) {
        ord.new <- c(ord.new, ord.old[-ord.new][coinc])
        coinc <- match(indp[ord.new[length(ord.new)], 2], indp[-ord.new, 1])
      }
    }
    
    indp <- indp[ord.new, , drop = FALSE]
    ahull.arcs <- cbind(arcs[ord.new, , drop = FALSE], indp)
    colnames(ahull.arcs) <- c("c1", "c2", "r", "v.x", "v.y", "theta", "end1", "end2")
    
    lengthah <- alphahull::lengthahull(arcs)
    
    addp <- noforget[is.na(match(noforget, indp))]
    num <- length(addp)
    if (num > 0) {
      mat.noforget <- cbind(
        matrix(ashape.obj$x[addp, 1:2], ncol = 2, byrow = FALSE),
        rep(0, num), rep(0, num), rep(0, num), rep(0, num),
        addp, addp
      )
      ahull.arcs <- rbind(ahull.arcs, mat.noforget)
    }
  } else {
    num <- length(noforget)
    ahull.arcs <- cbind(
      ashape.obj$x[noforget, ],
      rep(0, num), rep(0, num), rep(0, num), rep(0, num),
      noforget, noforget
    )
    colnames(ahull.arcs) <- c("c1", "c2", "r", "v.x", "v.y", "theta", "end1", "end2")
    lengthah <- 0
  }
  
  ahull.obj <- list(
    arcs = ahull.arcs,
    xahull = cutp,
    length = lengthah,
    complement = compl,
    alpha = alpha,
    ashape.obj = ashape.obj
  )
  class(ahull.obj) <- "ahull"
  invisible(ahull.obj)
}
