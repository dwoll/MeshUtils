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

#' @title 3D alpha wrapping
#' @description Reconstruction of a surface from a cloud of 3D points by
#'   alpha wrapping.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}},
#'    or a numeric matrix which stores the points, one point per row.
#' @param alphaRel Relative alpha parameter. The actual alpha parameter (see
#'   details) is defined as the length of the diagonal of the bounding box of
#'   the point cloud divided by the relative alpha parameter. Increase for
#'   more detailed mesh.
#' @param offsetRel Relative offset. The actual offset parameter (see details)
#'   is defined as the length of the diagonal of the bounding box of the
#'   point cloud divided by the relative offset parameter. Increase for
#'   output that is closer to input mesh.
#' @param normals Boolean: Whether to return vertex normals.
#'
#' @returns A \code{CGALmesh} object.
#'
#' @details See \url{https://doc.cgal.org/latest/Alpha_wrap_3/} for details.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @seealso \code{\link[MeshUtils]{reconstructAFS}},
#'    \code{\link[MeshUtils]{reconstructSSS}},
#'    \code{\link[MeshUtils]{reconstructPoisson}}
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' mesh            <- dataHeart1
#' mesh_rgl        <- toRGL(mesh)
#' mesh_alwrap     <- alphaWrap(mesh[["vertices"]],
#'                              alphaRel =5,
#'                              offsetRel=300)
#' mesh_alwrap_rgl <- toRGL(mesh_alwrap)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' wire3d(mesh_rgl)
#' next3d()
#' wire3d(mesh_alwrap_rgl)
#'
#' @export
#' @importFrom rgl tmesh3d
alphaWrap <- function(x, alphaRel, offsetRel, normals=FALSE) {
  stopifnot(isPositiveNumber(alphaRel))
  stopifnot(isPositiveNumber(offsetRel))
  stopifnot(isBoolean(normals))

  if(!inherits(x, "CGALmesh") && !is.matrix(x)) {
      stop("The `x` argument must be either of class 'CGALmesh'",
           " (i.e., the output of the `makeMesh()` function),",
           " or a numeric matrix with 3 columns.")
  }

  meshOut <- if(is.matrix(x)) {
    if(!is.numeric(x)) {
      stop("The `x` matrix must be numeric.", call. = TRUE)
    }
    if(ncol(x) != 3L) {
      stop("The `x` matrix must have three columns.", call. = TRUE)
    }
    if(nrow(x) <= 3L) {
      stop("Insufficient number of points.", call. = TRUE)
    }
    storage.mode(x) <- "double"
    if(anyNA(x)) {
      stop("Points in `x` with missing values are not allowed.", call. = TRUE)
    }
    alphaWrapPoints_cpp(t(x), alphaRel, offsetRel, normals)
  } else {
    meshCPP <- fromR(x)
    alphaWrapMesh_cpp(meshCPP, alphaRel, offsetRel, normals)
  }
  fromCPP(meshOut)
}
