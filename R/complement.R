#' Fast alpha-shape complement computation
#'
#' Computes the complement structure associated with an alpha shape using a
#' fast implementation designed to work with Delaunay--Voronoi meshes returned
#' by \code{\link{delvor_fast}} or other compatible \code{"delvor"} objects.
#'
#' This function mirrors the role of \code{\link[alphahull]{complement}} in the
#' \pkg{alphahull} workflow. It derives the circular arcs and half-plane
#' components that contribute to the alpha hull boundary from the Voronoi
#' representation of the point set.
#'
#' The implementation keeps the original logic of \pkg{alphahull} but reduces
#' overhead by vectorising several geometric operations, including midpoint,
#' radius, and between-point tests.
#'
#' @param x A \code{"delvor"} object, a two-column matrix of planar
#'   coordinates, or any object accepted by \code{\link[grDevices]{xy.coords}}.
#' @param y Optional vector of y coordinates if \code{x} contains only x
#'   coordinates.
#' @param alpha A non-negative numeric value controlling the level of detail of
#'   the alpha shape and alpha hull.
#'
#' @return A numeric matrix describing the complement components. Columns are:
#' \describe{
#'   \item{\code{c1}, \code{c2}}{Coordinates of the center of the circle or
#'   coefficients defining a supporting line for hull cases.}
#'   \item{\code{r}}{Radius of the complement circle, or a code used for hull
#'   half-plane cases.}
#'   \item{\code{ind1}, \code{ind2}}{Indices of the two data points associated
#'   with the Delaunay edge.}
#'   \item{\code{x1}, \code{y1}, \code{x2}, \code{y2}}{Coordinates of the two
#'   data points defining the Delaunay edge.}
#'   \item{\code{mx1}, \code{my1}, \code{mx2}, \code{my2}}{Coordinates of the
#'   two Voronoi endpoints associated with the edge.}
#'   \item{\code{bp1}, \code{bp2}}{Boundary indicators for the two Voronoi
#'   endpoints.}
#'   \item{\code{ind}}{Indicator identifying which side of the Voronoi edge the
#'   complement component is associated with.}
#'   \item{\code{v.x}, \code{v.y}}{Unit direction vector associated with the
#'   component.}
#'   \item{\code{theta}}{Half-angle subtended by the corresponding arc.}
#' }
#'
#' @details
#' The complement is an intermediate geometric structure used in the
#' construction of alpha hulls. It is derived from the dual Delaunay--Voronoi
#' representation and contains:
#' \itemize{
#'   \item circular components centered on Voronoi vertices,
#'   \item truncated components defined by intersections with the alpha radius,
#'   \item hull-specific half-plane components for unbounded Voronoi edges.
#' }
#'
#' This function is intended for use with \code{\link{ahull_fast3}} or other
#' workflows that rely on a fast and internally consistent mesh convention.
#' Because \code{\link{delvor_fast}} does not reproduce the raw row ordering and
#' hull conventions of \code{\link[alphahull]{delvor}} exactly, the output of
#' \code{complement_fast()} is designed to be consistent with
#' \code{\link{delvor_fast}} rather than identical row-for-row to
#' \code{\link[alphahull]{complement}}.
#'
#' @seealso
#' \code{\link{delvor_fast}}, \code{\link{ashape_fast}},
#' \code{\link[alphahull]{complement}}, \code{\link[alphahull]{ahull}}
#'
#' @examples
#' set.seed(1)
#' pts <- cbind(runif(100), runif(100))
#'
#' dv <- delvor_fast(pts)
#' comp <- complement_fast(dv, alpha = 0.2)
#'
#' dim(comp)
#' head(comp)
#'
#' @importFrom grDevices xy.coords
#' @export
complement_fast <- function(x, y = NULL, alpha) {
  if (alpha < 0) {
    stop("Parameter alpha must be greater or equal to zero")
  }
  
  if (!inherits(x, "delvor")) {
    dd.obj <- delvor_fast(x, y)
  } else {
    dd.obj <- x
  }
  
  mesh <- dd.obj$mesh
  
  dm1 <- sqrt((mesh[, "x1"] - mesh[, "mx1"])^2 + (mesh[, "y1"] - mesh[, "my1"])^2)
  dm2 <- sqrt((mesh[, "x1"] - mesh[, "mx2"])^2 + (mesh[, "y1"] - mesh[, "my2"])^2)
  
  pm.x <- (mesh[, "x1"] + mesh[, "x2"]) * 0.5
  pm.y <- (mesh[, "y1"] + mesh[, "y2"]) * 0.5
  dm   <- sqrt((mesh[, "x1"] - mesh[, "x2"])^2 + (mesh[, "y1"] - mesh[, "y2"])^2) * 0.5
  
  d.pm.m1 <- sqrt((mesh[, "mx1"] - pm.x)^2 + (mesh[, "my1"] - pm.y)^2)
  d.pm.m2 <- sqrt((mesh[, "mx2"] - pm.x)^2 + (mesh[, "my2"] - pm.y)^2)
  
  theta.m1 <- atan(dm / d.pm.m1)
  theta.m2 <- atan(dm / d.pm.m2)
  
  v.x.m1 <- ifelse(d.pm.m1 != 0, (pm.x - mesh[, "mx1"]) / d.pm.m1, 0)
  v.y.m1 <- ifelse(d.pm.m1 != 0, (pm.y - mesh[, "my1"]) / d.pm.m1, 0)
  v.x.m2 <- ifelse(d.pm.m2 != 0, (pm.x - mesh[, "mx2"]) / d.pm.m2, 0)
  v.y.m2 <- ifelse(d.pm.m2 != 0, (pm.y - mesh[, "my2"]) / d.pm.m2, 0)
  
  thetam.m1 <- acos(pmin(1, pmax(-1, v.x.m1)))
  thetam.m2 <- acos(pmin(1, pmax(-1, v.x.m2)))
  
  n.edges <- nrow(mesh)
  
  ## exact same betw semantics as original rank()[3] == 2
  betw <- numeric(n.edges)
  vert <- mesh[, "mx1"] == mesh[, "mx2"]
  
  if (sum(vert) > 0) {
    a <- mesh[vert, "my1"]
    b <- mesh[vert, "my2"]
    p <- pm.y[vert]
    betw[vert] <- ((p > pmin(a, b)) & (p < pmax(a, b))) | ((a == b) & (p == a))
  }
  
  if (sum(!vert) > 0) {
    a <- mesh[!vert, "mx1"]
    b <- mesh[!vert, "mx2"]
    p <- pm.x[!vert]
    betw[!vert] <- ((p > pmin(a, b)) & (p < pmax(a, b))) | ((a == b) & (p == a))
  }
  
  aux <- alpha^2 - dm^2
  
  make_block <- function(idx, center_id, ex, ey, rr, vx, vy, th) {
    if (!any(idx)) return(NULL)
    
    m <- mesh[idx, , drop = FALSE]
    cbind(
      ex, ey, rr,
      m,
      center_id,
      vx, vy, th
    )
  }
  
  ## ---- side 1 ----
  case1.1 <- mesh[, "bp1"] == 0 & dm1 > alpha
  comp1.1 <- make_block(
    case1.1, 1,
    mesh[case1.1, "mx1"],
    mesh[case1.1, "my1"],
    dm1[case1.1],
    v.x.m1[case1.1],
    v.y.m1[case1.1],
    theta.m1[case1.1]
  )
  
  case1.2 <- (mesh[, "bp1"] == 0 & dm1 > alpha & aux > 0 & betw == 1) |
    (mesh[, "bp1"] == 1 & aux > 0 & betw == 1)
  if (any(case1.2)) {
    a1 <- sqrt(aux[case1.2])
    e.x <- pm.x[case1.2] - a1 * v.x.m1[case1.2]
    e.y <- pm.y[case1.2] - a1 * v.y.m1[case1.2]
    theta <- atan(dm[case1.2] / a1)
    comp1.2 <- make_block(case1.2, 1, e.x, e.y, alpha, v.x.m1[case1.2], v.y.m1[case1.2], theta)
  } else {
    comp1.2 <- NULL
  }
  
  case1.3 <- (mesh[, "bp1"] == 0 & dm1 > alpha & aux > 0 & betw == 0) |
    (mesh[, "bp1"] == 1 & aux > 0 & betw == 0)
  if (any(case1.3)) {
    a1 <- sqrt(aux[case1.3])
    e.x <- pm.x[case1.3] - a1 * v.x.m1[case1.3]
    e.y <- pm.y[case1.3] - a1 * v.y.m1[case1.3]
    theta <- atan(dm[case1.3] / a1)
    keep <- alpha - a1 < dm2[case1.3] - d.pm.m2[case1.3]
    
    if (any(keep)) {
      sub_idx <- which(case1.3)[keep]
      comp1.3.1 <- cbind(
        e.x[keep], e.y[keep], alpha,
        mesh[sub_idx, , drop = FALSE],
        1,
        v.x.m1[sub_idx],
        v.y.m1[sub_idx],
        theta[keep]
      )
    } else {
      comp1.3.1 <- NULL
    }
  } else {
    comp1.3.1 <- NULL
  }
  
  ## ---- side 2 ----
  case2.1 <- mesh[, "bp2"] == 0 & dm2 > alpha
  comp2.1 <- make_block(
    case2.1, 2,
    mesh[case2.1, "mx2"],
    mesh[case2.1, "my2"],
    dm2[case2.1],
    v.x.m2[case2.1],
    v.y.m2[case2.1],
    theta.m2[case2.1]
  )
  
  case2.2 <- (mesh[, "bp2"] == 0 & dm2 > alpha & aux > 0 & betw == 1) |
    (mesh[, "bp2"] == 1 & aux > 0 & betw == 1)
  if (any(case2.2)) {
    a1 <- sqrt(aux[case2.2])
    e.x <- pm.x[case2.2] - a1 * v.x.m2[case2.2]
    e.y <- pm.y[case2.2] - a1 * v.y.m2[case2.2]
    theta <- atan(dm[case2.2] / a1)
    comp2.2 <- make_block(case2.2, 2, e.x, e.y, alpha, v.x.m2[case2.2], v.y.m2[case2.2], theta)
  } else {
    comp2.2 <- NULL
  }
  
  case2.3 <- (mesh[, "bp2"] == 0 & dm2 > alpha & aux > 0 & betw == 0) |
    (mesh[, "bp2"] == 1 & aux > 0 & betw == 0)
  if (any(case2.3)) {
    a1 <- sqrt(aux[case2.3])
    e.x <- pm.x[case2.3] - a1 * v.x.m2[case2.3]
    e.y <- pm.y[case2.3] - a1 * v.y.m2[case2.3]
    theta <- atan(dm[case2.3] / a1)
    keep <- alpha - a1 < dm1[case2.3] - d.pm.m1[case2.3]
    
    if (any(keep)) {
      sub_idx <- which(case2.3)[keep]
      comp2.3.1 <- cbind(
        e.x[keep], e.y[keep], alpha,
        mesh[sub_idx, , drop = FALSE],
        2,
        v.x.m2[sub_idx],
        v.y.m2[sub_idx],
        theta[keep]
      )
    } else {
      comp2.3.1 <- NULL
    }
  } else {
    comp2.3.1 <- NULL
  }
  
  ## ---- hull-specific block ----
  case.h <- mesh[, "bp2"] == 1
  if (any(case.h)) {
    idx <- which(case.h)
    comp.h <- matrix(NA_real_, nrow = length(idx), ncol = 19)
    
    for (k in seq_along(idx)) {
      i <- idx[k]
      
      if (mesh[i, "x1"] == mesh[i, "x2"]) {
        a <- mesh[i, "x1"]
        sig <- if (mesh[i, "mx2"] <= a) -4 else -3
        comp.h[k, ] <- c(a, 0, sig, mesh[i, ], 2, 0, 0, 0)
      } else {
        b <- (mesh[i, "y2"] - mesh[i, "y1"]) / (mesh[i, "x2"] - mesh[i, "x1"])
        a <- mesh[i, "y1"] - mesh[i, "x1"] * b
        sig <- if (mesh[i, "my2"] <= a + b * mesh[i, "mx2"]) -2 else -1
        comp.h[k, ] <- c(a, b, sig, mesh[i, ], 2, 0, 0, 0)
      }
    }
  } else {
    comp.h <- NULL
  }
  
  comp <- rbind(comp1.1, comp1.2, comp1.3.1, comp2.1, comp2.2, comp2.3.1, comp.h)
  colnames(comp) <- c("c1", "c2", "r", colnames(mesh), "ind", "v.x", "v.y", "theta")
  rownames(comp) <- NULL
  
  invisible(comp)
}