## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/Boov/
## https://github.com/stla/PolygonSoup/
## https://github.com/stla/cgalMeshes/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## adapted by
## Daniel Wollschlaeger
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title Smooth shape
#' @description Smooths the overall shape of the mesh by using the mean
#'   curvature flow. The mesh must be triangle.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#'   The mesh must be triangle or be able to be made triangle.
#' @param indices The indices of the faces to be smoothed. If missing, the whole mesh
#'   is smoothed.
#' @param nIter Positive \code{integer}: Number of iterations.
#' @param time Positive number: A time step that corresponds to the speed by
#'   which the surface is smoothed (the larger the faster); typical values
#'   lie between \code{1e-6} and \code{1}.
#' @param triangulate Boolean: Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @returns A \code{CGALmesh} object.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' mesh       <- makeMesh(mesh=dataTruncIcosahedron, triangulate=TRUE)
#' mesh_rgl   <- toRGL(mesh)
#' mesh_s     <- smoothShape(mesh, nIter=2, time=0.01)
#' mesh_s_rgl <- toRGL(mesh_s)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_rgl)
#' next3d()
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_s_rgl)
#'
#' @export
smoothShape <- function(x, indices, nIter = 1, time = 0.0001, triangulate = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isPositiveNumber(time))
  stopifnot(isBoolean(triangulate))
  if(missing(indices)) {
    indices <- integer(0L)
  } else {
    stopifnot(isAtomicVector(indices))
    stopifnot(is.numeric(indices))
    integers <- isTRUE(all.equal(indices, floor(indices)))
    if(!integers) {
      stop("The indices must be positive integers.")
    }
    positive <- all(indices >= 1L)
    if(!positive) {
      stop("The indices must be positive integers.")
    }
    indices <- unique(as.integer(indices)) - 1L
  }
  meshCPP <- fromR(x)
  meshOut <- smoothShape_cpp(meshCPP, indices, as.integer(nIter), time, triangulate)
  fromCPP(meshOut)
}
