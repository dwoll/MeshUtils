## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/Boov/
## https://github.com/stla/PolygonSoup/
## https://github.com/stla/cgalMeshes/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title Poisson surface reconstruction
#' @description Poisson reconstruction of a surface, from a cloud of 3D points.
#'
#' @param points numeric matrix which stores the points, one point per row
#' @param normals numeric matrix which stores the normals, one normal per row
#'   (it must have the same size as the \code{points} matrix); if you don't
#'   have normals, set \code{normals=NULL} (the default) and some normals will
#'   be computed with the help of \code{\link[Rvcg]{vcgUpdateNormals}}, or
#'   use the \code{\link{getSomeNormals}} function
#' @param spacing size parameter; smaller values increase the precision of the
#'   output mesh at the cost of higher computation time; set to \code{NULL}
#'   (the default) for a reasonable automatic value: an average spacing whose
#'   value will be displayed in a message and that you can also get in the
#'   \code{"spacing"} attribute of the output
#' @param sm_angle bound for the minimum facet angle in degrees
#' @param sm_radius relative bound for the radius of the surface Delaunay balls
#' @param sm_distance relative bound for the center-center distances
#'
#' @return A triangle mesh, of class \code{mesh3d} (ready for plotting
#'   with \strong{rgl}).
#'
#' @details See \href{https://doc.cgal.org/latest/Poisson_surface_reconstruction_3/index.html}{Poisson Surface Reconstruction}.
#'
#' @export
#' @importFrom rgl tmesh3d
#' @importFrom Rvcg vcgUpdateNormals
#'
#'
#' @examples
#' library(SurfaceReconstruction)
#' library(rgl)
#'
#' # Solid Möbius strip
#' Psr_mesh <- PoissonReconstruction(SolidMobiusStrip)
#' shade3d(Psr_mesh, color= "yellow")
#' wire3d(Psr_mesh, color = "black")
#'
#' # Hopf torus
#' Psr_mesh <- PoissonReconstruction(HopfTorus, spacing = 0.2)
#' shade3d(Psr_mesh, color= "darkorange")
#' wire3d(Psr_mesh, color = "black")
reconstructPoisson <- function(
  points,
  normals = NULL,
  spacing = NULL,
  smAngle = 20,
  smRadius = 30,
  smDistance = 0.375
) {
  if(!is.matrix(points) || !is.numeric(points)) {
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  if(anyNA(points)) {
    stop("Missing values in the `points` matrix are not allowed.", call. = TRUE)
  }
  dimension <- ncol(points)
  if(dimension != 3L) {
    stop("The `points` matrix must have three columns.", call. = TRUE)
  }
  storage.mode(points) <- "double"
  if(is.null(normals)) {
    normals <- vcgUpdateNormals(points, silent = TRUE)[["normals"]][-4L, ]
  } else if(is.function(normals) && inherits(normals, "CGALnormalsFunc")) {
    normals <- normals(points)
  } else {
    stop(
      "Invalid argument `normals`: it must be `NULL` or a function returned ",
      "by the `getSomeNormals` function."
    )
  }
  if(nrow(points) <= dimension) {
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
  psr <- reconstructPoisson_cpp(
    t(points), normals, spacing, smAngle, smRadius, smDistance
  )
  out <- vcgUpdateNormals(
    tmesh3d(psr[["vertices"]],
            psr[["faces"]],
            normals = NULL,
            homogeneous = FALSE)
  )
  if(spacing == -1) {
    message(sprintf(
      "Poisson reconstruction using average spacing: %s.",
      formatC(psr[["spacing"]])
    ))
    attr(out, "spacing") <- psr[["spacing"]]
  }
  out
}
