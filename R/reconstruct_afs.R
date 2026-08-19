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

#' @title Advancing front surface reconstruction
#' @description Reconstruction of a surface from a cloud of 3D points.
#'
#' @param x Numeric matrix which stores the points, one point per row.
#' @param jetSmoothing If not \code{NULL}, must be an integer >= 2.
#'   Then, the point cloud is smoothed before the reconstruction, using
#'   this integer as the number of neighbors for the smoothing. Note that this
#'   smoothing preprocessing relocates the points and then should not be used
#'   if the points have been sampled without noise on the surface.
#' @param repairSoup Boolean. Attempt to fix polygon soup?
#'
#' @returns A \code{CGALmesh} object.
#'
#' @details See \href{https://doc.cgal.org/latest/Advancing_front_surface_reconstruction/index.html#Chapter_Advancing_Front_Surface_Reconstruction}{Advancing Front Surface Reconstruction}.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @seealso \code{\link[MeshUtils]{reconstructSSS}},
#'    \code{\link[MeshUtils]{reconstructPoisson}}, \code{\link[MeshUtils]{alphaWrap}}
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' # no smoothing
#' mesh          <- makeMesh(mesh=dataHopfTorus)
#' mesh_afs1     <- reconstructAFS(mesh[["vertices"]])
#' mesh_afs1_rgl <- toRGL(mesh_afs1)
#'
#' # jet smoothing
#' mesh_afs2     <- reconstructAFS(mesh[["vertices"]],
#'                                 jetSmoothing=30)
#' mesh_afs2_rgl <- toRGL(mesh_afs2)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(20, -40, zoom=0.85)
#' shade3d(mesh_afs1_rgl, color="gold")
#' next3d()
#' shade3d(mesh_afs2_rgl, color="gold")
#'
#' @export
reconstructAFS <- function(x, jetSmoothing=NULL, repairSoup=TRUE) {
  if(!is.matrix(x) || !is.numeric(x)) {
    stop("The `x` argument must be a numeric matrix.", call. = TRUE)
  }
  if(ncol(x) != 3L) {
    stop("The `x` matrix must have three columns.", call. = TRUE)
  }
  if(nrow(x) <= 3L) {
    stop("Insufficient number of points in `x`.", call. = TRUE)
  }
  if(anyNA(x)) {
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  storage.mode(x) <- "double"
  if(!is.null(jetSmoothing)) {
    stopifnot(isPositiveInteger(jetSmoothing), jetSmoothing >= 2L)
  } else {
    jetSmoothing <- 0L
  }
  stopifnot(isBoolean(repairSoup))
  mesh_cpp <- reconstructAFS_cpp(t(x), as.integer(jetSmoothing), repairSoup)
  fromCPP(mesh_cpp)
}
