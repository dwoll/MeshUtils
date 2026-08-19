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

#' @title Scale-space surface reconstruction
#' @description Reconstruction of a surface from a cloud of 3D points.
#'
#' @param x Numeric matrix which stores the points, one point per row.
#' @param scaleIterations Number of iterations used to increase the scale.
#' @param neighbors Number of neighbors used to smooth the point cloud.
#' @param samples Number of samples used to smooth the point cloud.
#' @param separateShells Boolean, whether to separate the shells.
#' @param forceManifold Boolean, whether to force a manifold output mesh.
#' @param borderAngle Bound on the angle in degrees used to detect border edges.
#' @param repairSoup Boolean. Attempt to fix polygon soup?
#'
#' @returns A \code{CGALmesh} object or a \code{\link[rgl]{mesh3d}} object from package \strong{rgl}.
#'
#' @details See \href{https://doc.cgal.org/latest/Scale_space_reconstruction_3/index.html}{Scale-space Surface Reconstruction}.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @seealso \code{\link[MeshUtils]{reconstructAFS}},
#'    \code{\link[MeshUtils]{reconstructPoisson}}, \code{\link[MeshUtils]{alphaWrap}}
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' mesh     <- dataHeart1
#' mesh_rgl <- toRGL(mesh)
#' mesh_sss <- reconstructSSS(mesh[["vertices"]],
#'                            scaleIterations=1,
#'                            forceManifold  =TRUE,
#'                            neighbors      =6)
#'
#' mesh_sss_rgl <- toRGL(mesh_sss)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' wire3d(mesh_rgl)
#' next3d()
#' wire3d(mesh_sss_rgl)
#'
#' @export
reconstructSSS <- function(
  x,
  scaleIterations=1,
  neighbors      =12,
  samples        =300,
  separateShells =FALSE,
  forceManifold  =TRUE,
  borderAngle    =45,
  repairSoup     =TRUE) {
  if(!is.matrix(x) || !is.numeric(x)) {
    stop("The `x` argument must be a numeric matrix.", call. = TRUE)
  }
  if(ncol(x) != 3L) {
    stop("The `x` matrix must have three columns.", call. = TRUE)
  }
  if(nrow(x) <= 3L) {
    stop("Insufficient number of points.", call. = TRUE)
  }
  if(anyNA(x)) {
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  storage.mode(x) <- "double"
  stopifnot(isStrictPositiveInteger(scaleIterations))
  stopifnot(isPositiveInteger(neighbors), neighbors >= 2)
  stopifnot(isBoolean(separateShells))
  stopifnot(isBoolean(forceManifold))
  stopifnot(isNonNegativeNumber(borderAngle))
  stopifnot(isBoolean(repairSoup))
  mesh_cpp <- reconstructSSS_cpp(
    t(x),
    as.integer(scaleIterations),
    as.integer(neighbors),
    as.integer(samples),
    separateShells,
    forceManifold,
    as.double(borderAngle),
    repairSoup)
  fromCPP(mesh_cpp)
}
