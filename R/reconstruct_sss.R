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
#' @param out Character to indicate output mesh format.
#'
#' @returns A \code{CGALmesh} object or a \code{mesh3d} object from package \strong{rgl}.
#'
#' @details See \href{https://doc.cgal.org/latest/Scale_space_reconstruction_3/index.html}{Scale-space Surface Reconstruction}.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @seealso [reconstructAFS()], [reconstructPoisson()], [alphaWrap()]
#'
#' @examples
#' library(MeshUtils)
#' mesh   <- makeMesh(mesh=HopfTorus)
#' mesh_r <- reconstructSSS(
#'   mesh,
#'   scaleIterations=4,
#'   forceManifold  =TRUE,
#'   neighbors      =30,
#'   out            ="rgl")
#'
#' library(rgl)
#' open3d(windowRect=50 + c(0, 0, 512, 512))
#' view3d(20, -40, zoom=0.85)
#' shade3d(mesh_r, color="tomato")
#'
#' @export
#' @importFrom Rvcg vcgUpdateNormals
#' @importFrom rgl tmesh3d
reconstructSSS <- function(
  x,
  scaleIterations=1,
  neighbors      =12,
  samples        =300,
  separateShells =FALSE,
  forceManifold  =TRUE,
  borderAngle    =45,
  out            =c("CGALmesh", "rgl")
) {
  out <- match.arg(out)
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
  mesh_r <- reconstructSSS_cpp(
    t(x),
    as.integer(scaleIterations),
    as.integer(neighbors),
    as.integer(samples),
    separateShells,
    forceManifold,
    as.double(borderAngle))
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
