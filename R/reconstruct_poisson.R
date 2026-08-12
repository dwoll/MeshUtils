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
#' @param out Character to indicate output mesh format.
#'
#' @returns A \code{CGALmesh} object or a \code{\link[rgl]{mesh3d}} object from package \strong{rgl}.
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
#' # Hopf torus
#' mesh   <- makeMesh(mesh=HopfTorus)
#' mesh_r <- reconstructPoisson(mesh[["vertices"]], spacing=0.2, out="rgl")
#' shade3d(mesh_r, color="darkorange")
#' wire3d(mesh_r,  color="black")
#'
#' @export
#' @importFrom rgl tmesh3d
#' @importFrom Rvcg vcgUpdateNormals
reconstructPoisson <- function(
  x,
  normals   = NULL,
  spacing   = NULL,
  smAngle   = 20,
  smRadius  = 30,
  smDistance= 0.375,
  out       =c("CGALmesh", "rgl")
) {
  out <- match.arg(out)
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
    stop(
      "Invalid argument `normals`: it must be `NULL` or a function returned ",
      "by the `getSomeNormals` function."
    )
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
  mesh_r <- reconstructPoisson_cpp(
    t(x), normals, spacing, smAngle, smRadius, smDistance)
  mesh_rwn <- vcgUpdateNormals(
    tmesh3d(mesh_r[["vertices"]],
            mesh_r[["faces"]],
            normals    =NULL,
            homogeneous=FALSE))
  if(spacing == -1) {
    message(sprintf(
      "Poisson reconstruction using average spacing: %s.",
      formatC(mesh_r[["spacing"]])))
    attr(mesh_rwn, "spacing") <- mesh_r[["spacing"]]
  }
  mesh_out <- if(out == "rgl") {
    mesh_rwn
  } else {
    makeMesh(mesh=mesh_rwn)
  }
  mesh_out
}
