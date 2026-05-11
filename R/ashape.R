#' Fast alpha-shape computation with vectorised edge filtering
#'
#' A drop-in replacement for \code{alphahull::ashape()} that preserves the
#' original algorithm but removes the main row-wise bottlenecks in the edge
#' filtering stage.
#'
#' @param x A \code{delvor} object, a two-column matrix of coordinates, or any
#'   input accepted by \code{grDevices::xy.coords()}.
#' @param y Optional y coordinates.
#' @param alpha Non-negative alpha value.
#'
#' @return An object of class \code{"ashape"}.
#' @export
ashape_fast <- function(x, y = NULL, alpha) {
  if (alpha < 0) {
    stop("Parameter alpha must be greater or equal to zero")
  }
  
  if (!inherits(x, "delvor")) {
    dd.obj <- delvor_fast(x, y)
  } else {
    dd.obj <- x
  }
  
  xy.data <- dd.obj$x
  mesh <- dd.obj$mesh
  
  dm1 <- sqrt((mesh[, "x1"] - mesh[, "mx1"])^2 + (mesh[, "y1"] - mesh[, "my1"])^2)
  dm2 <- sqrt((mesh[, "x1"] - mesh[, "mx2"])^2 + (mesh[, "y1"] - mesh[, "my2"])^2)
  
  dm1[mesh[, "bp1"] == 1] <- Inf
  dm2[mesh[, "bp2"] == 1] <- Inf
  
  n <- nrow(xy.data)
  ind <- seq_len(n)
  ind.on <- chull(xy.data)
  ind.in <- ind[-ind.on]
  
  if (length(ind.in) > 0) {
    aux_alpha <- rbind(
      cbind(mesh[, c("ind1", "ind2")], dm1),
      cbind(mesh[, c("ind1", "ind2")], dm2),
      cbind(mesh[, c("ind2", "ind1")], dm1),
      cbind(mesh[, c("ind2", "ind1")], dm2)
    )
    fc <- factor(aux_alpha[, 1])
    aux_max <- tapply(aux_alpha[, 3], fc, max)
    alpha.max <- cbind(aux_max, as.numeric(levels(fc)))
    alpha.ext <- c(
      ind.on,
      ind.in[na.omit(match(alpha.max[alpha < alpha.max[, 1], 2], ind.in))]
    )
  } else {
    alpha.ext <- ind.on
  }
  
  i1 <- match(mesh[, 1], alpha.ext)
  i2 <- match(mesh[, 2], alpha.ext)
  is.edge <- which(i1 & i2)
  
  aux <- mesh[is.edge, , drop = FALSE]
  pm.x <- (aux[, "x1"] + aux[, "x2"]) * 0.5
  pm.y <- (aux[, "y1"] + aux[, "y2"]) * 0.5
  dm   <- sqrt((aux[, "x1"] - aux[, "x2"])^2 + (aux[, "y1"] - aux[, "y2"])^2) * 0.5
  
  betw <- rep(NA_real_, nrow(aux))
  vertical <- aux[, "mx1"] == aux[, "mx2"]
  
  # Exact rank(...) semantics:
  # TRUE if strictly between the two endpoints,
  # or if all three values are exactly equal.
  if (any(vertical)) {
    idx <- which(vertical)
    a <- aux[idx, "my1"]
    b <- aux[idx, "my2"]
    p <- pm.y[idx]
    
    ok <- ((p > pmin(a, b)) & (p < pmax(a, b))) | ((a == b) & (p == a))
    betw[idx[ok]] <- 1
  }
  
  if (any(!vertical)) {
    idx <- which(!vertical)
    a <- aux[idx, "mx1"]
    b <- aux[idx, "mx2"]
    p <- pm.x[idx]
    
    ok <- ((p > pmin(a, b)) & (p < pmax(a, b))) | ((a == b) & (p == a))
    betw[idx[ok]] <- 1
  }
  
  z3 <- dm * betw
  l.min <- pmin(dm1[is.edge], dm2[is.edge], z3, na.rm = TRUE)
  l.max <- pmax(dm1[is.edge], dm2[is.edge], z3, na.rm = TRUE)
  
  in.ashape <- (l.min <= alpha & alpha <= l.max)
  
  edges <- aux[in.ashape, , drop = FALSE]
  edges <- as.matrix(edges)
  colnames(edges) <- colnames(aux)
  rownames(edges) <- NULL
  
  ashape.obj <- list(
    edges = edges,
    length = sum(2 * dm[in.ashape]),
    alpha = alpha,
    alpha.extremes = alpha.ext,
    delvor.obj = dd.obj,
    x = xy.data
  )
  
  class(ashape.obj) <- "ashape"
  invisible(ashape.obj)
}