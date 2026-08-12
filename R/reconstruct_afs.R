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
#' @param jetSmoothing If not \code{NULL}, must be an integer higher than two.
#'   Then, the point cloud is smoothed before the reconstruction, using
#'   this integer as the number of neighbors for the smoothing. Note that this
#'   smoothing preprocessing relocates the points and then should not be used
#'   if the points have been sampled without noise on the surface.
#' @param out Character to indicate output mesh format.
#'
#' @returns A \code{CGALmesh} object or a \code{mesh3d} object from package \strong{rgl}.
#'
#' @details See \href{https://doc.cgal.org/latest/Advancing_front_surface_reconstruction/index.html#Chapter_Advancing_Front_Surface_Reconstruction}{Advancing Front Surface Reconstruction}.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @seealso [reconstructSSS()], [reconstructPoisson()], [alphaWrap()]
#'
#' @examples
#' library(MeshUtils)
#' # no smoothing
#' mesh    <- makeMesh(mesh=HopfTorus)
#' mesh_r1 <- reconstructAFS(mesh, out="rgl")
#'
#' # jet smoothing
#' mesh_r2 <- reconstructAFS(mesh, jetSmoothing=30, out="rgl")
#' # plot
#' library(rgl)
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(20, -40, zoom=0.85)
#' shade3d(mesh_r1, color="gold")
#' next3d()
#' view3d(20, -40, zoom=0.85)
#' shade3d(mesh_r2, color="gold")
#'
#' @export
#' @importFrom Rvcg vcgUpdateNormals
#' @importFrom rgl tmesh3d
reconstructAFS <- function(x, jetSmoothing=NULL, out=c("CGALmesh", "rgl")) {
  out <- match.arg(out)
  if(!is.matrix(x) || !is.numeric(x)){
    stop("The `x` argument must be a numeric matrix.", call. = TRUE)
  }
  if(ncol(x) != 3L) {
    stop("The `x` matrix must have three columns.", call. = TRUE)
  }
  if(nrow(x) <= 3L) {
    stop("Insufficient number of x.", call. = TRUE)
  }
  if(anyNA(x)){
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  storage.mode(x) <- "double"
  if(!is.null(jetSmoothing)) {
    stopifnot(isPositiveInteger(jetSmoothing), jetSmoothing >= 2L)
  } else {
    jetSmoothing <- 0L
  }
  mesh_r   <- reconstructAFS_cpp(t(x), as.integer(jetSmoothing))
  mesh_rwn <- vcgUpdateNormals(
    tmesh3d(mesh_r[["vertices"]],
            mesh_r[["faces"]],
            normals    =NULL,
            homogeneous=FALSE))
  mesh_out <- if(out == "rgl") {
    mesh_rwn
  } else {
    makeMesh(mesh=mesh_rwn)
  }
  mesh_out
}
