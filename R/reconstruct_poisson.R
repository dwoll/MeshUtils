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

#' @title Poisson surface reconstruction
#' @description Poisson reconstruction of a surface, from a cloud of 3D points.
#'
#' @param x Numeric matrix which stores the points, one point per row.
#' @param normals Numeric matrix which stores the normals, one normal per row
#'   (it must have the same size as the \code{x} matrix). If you do not
#'   have normals, set \code{normals=NULL} (the default), and some normals will
#'   be computed with the help of \code{\link[Rvcg]{vcgUpdateNormals}}, or
#'   use the \code{\link[MeshUtils]{getSomeNormals}} function.
#' @param spacing Size parameter. Smaller values increase the precision of the
#'   output mesh at the cost of higher computation time. Set to \code{NULL}
#'   (the default) for a reasonable automatic value: an average spacing whose
#'   value will be displayed in a message and that you can also get in the
#'   \code{"spacing"} attribute of the output.
#' @param smAngle Bound for the minimum facet angle in degrees.
#' @param smRadius Relative bound for the radius of the surface Delaunay balls.
#' @param smDistance Relative bound for the center-center distances.
#' @param clean Boolean. Attempt to fix polygon soup?
#'
#' @returns A \code{CGALmesh} object.
#'
#' @details See \href{https://doc.cgal.org/latest/Poisson_surface_reconstruction_3/index.html}{Poisson Surface Reconstruction}.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @seealso \code{\link[MeshUtils]{reconstructAFS}},
#'    \code{\link[MeshUtils]{reconstructSSS}},
#'    \code{\link[MeshUtils]{alphaWrap}}

#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' mesh     <- makeMesh(mesh=dataHopfTorus)
#' mesh_rgl <- toRGL(mesh)
#' mesh_psr <- reconstructPoisson(mesh[["vertices"]],
#'                                smAngle=10,
#'                                smRadius=3,
#'                                smDistance=0.3)
#' mesh_psr_rgl <- toRGL(mesh_psr)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' wire3d(mesh_rgl)
#' next3d()
#' wire3d(mesh_psr_rgl)
#'
#' @export
#' @importFrom Rvcg vcgUpdateNormals
reconstructPoisson <- function(
  x,
  normals   = NULL,
  spacing   = NULL,
  smAngle   = 20,
  smRadius  = 30,
  smDistance= 0.375,
  clean     =TRUE) {
  if(!is.matrix(x) || !is.numeric(x)) {
    stop("The `x` argument must be a numeric matrix.", call. = TRUE)
  }
  if(anyNA(x)) {
    stop("Missing values in the `x` matrix are not allowed.", call. = TRUE)
  }
  dimension <- ncol(x)
  if(dimension != 3L) {
    stop("The `x` matrix must have three columns.", call. = TRUE)
  }
  storage.mode(x) <- "double"
  if(is.null(normals)) {
    normals <- vcgUpdateNormals(x, silent = TRUE)[["normals"]][-4L, ]
  } else if(is.function(normals) && inherits(normals, "CGALnormalsFunc")) {
    normals <- normals(x)
  } else {
    stop("Invalid argument `normals`: it must be `NULL` or a function ",
         "returned by the `getSomeNormals` function.")
  }
  if(nrow(x) <= dimension) {
    stop("Insufficient number of points.", call. = TRUE)
  }
  # if(any(is.na(points)) || (!is.null(normals) && any(is.na(normals)))) {
  #   stop("Points or normals with missing values are not allowed.", call. = TRUE)
  # }
  if(is.null(spacing)) {
    spacing <- -1
  } else {
    stopifnot(isPositiveNumber(spacing))
  }
  stopifnot(isPositiveNumber(smAngle))
  stopifnot(isPositiveNumber(smRadius))
  stopifnot(isPositiveNumber(smDistance))
  mesh_cpp <- reconstructPoisson_cpp(
    t(x), normals, spacing, smAngle, smRadius, smDistance, clean)
  fromCPP(mesh_cpp)
}
