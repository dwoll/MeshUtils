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

#' @title Hausdorff distance between two meshes
#' @description Hausdorff distance between two meshes. Either
#'   approximate distance, or distance estimate with a given error bound.
#' @param mesh1 A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param mesh2 A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param triangulate1 Boolean: Whether to triangulate the faces of \code{mesh1}.
#'   Ignored if faces are already triangle.
#' @param triangulate2 Boolean: Whether to triangulate the faces of \code{mesh2}.
#'   Ignored if faces are already triangle.
#' @param symmetric Boolean: Whether to consider the symmetric Hausdorff
#'   distance.
#' @param errorBound A positive number, a bound on the error of the
#'   estimate. If missing, the approximate distance is returned.
#' @returns A number. For the apprixmate distance, the algorithm uses
#'   some simulations and thus the result can vary.
#' @details See \url{https://doc.cgal.org/latest/Polygon_mesh_processing/index.html#PMPDistance} for details.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' ## approximate symmetric Hausdorff distance
#' getHausdorffDistance(dataHeart1, dataHeart2, symmetric=TRUE)
#'
#' ## estimate with error bound
#' getHausdorffDistance(dataHeart1, dataHeart2, symmetric=TRUE,
#'                      errorBound=0.001)
#'
#' @export
getHausdorffDistance <- function(mesh1, mesh2, symmetric = TRUE, errorBound,
	triangulate1 = FALSE, triangulate2 = FALSE) {
  stopifnot(inherits(mesh1, "CGALmesh"))
  stopifnot(inherits(mesh2, "CGALmesh"))
  stopifnot(isBoolean(symmetric))
  stopifnot(isBoolean(triangulate1))
  stopifnot(isBoolean(triangulate1))
  meshCPP1 <- fromR(mesh1)
  meshCPP2 <- fromR(mesh2)
  if(!missing(errorBound)) {
    stopifnot(isPositiveNumber(errorBound))
    getHausdorffEst_cpp(meshCPP1, meshCPP2, symmetric, errorBound,
                        triangulate1, triangulate2)
  } else {
    getHausdorffApprox_cpp(meshCPP1, meshCPP2, symmetric,
    	                     triangulate1, triangulate2)
  }
}
