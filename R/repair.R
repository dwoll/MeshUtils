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

#' @title Remove self intersections
#' @description Remove self intersections.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param method \code{character}. One of \code{"auto"} (for auto-refine) and
#'   \code{"auto_snap"} (auto-refine with iterative snap). See details.
#' @param normals Boolean. Whether to return vertex normals.
#' @param verbose Boolean. Whether to print out messages about mesh processing.
#' @returns \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#' @details See \url{https://www.cgal.org/2025/06/13/autorefine-and-snap/}.
#'     If faces are not triangle, the mesh is triangulated.
#'
#' @examples
#' library(MeshUtils)
#' mesh     <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' mesh_nsi <- removeSelfIntersections(mesh)
#' getVolume(mesh_nsi)
#'
#' @export
removeSelfIntersections <- function(
  x, method=c("auto", "auto_snap"), normals = FALSE, verbose = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(normals))
  stopifnot(isBoolean(verbose))
  method_choices <- c("auto", "auto_snap")
  method    <- match.arg(method, choices=method_choices)
  methodInt <- match(method, method_choices)
  meshCPP   <- fromR(x)
  mesh      <- removeSelfIntersections_cpp(meshCPP, methodInt, normals, verbose)
  fromCPP(mesh)
}

#' @title Fill boundary holes
#' @description Fill boundary holes.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param fairHole Boolean. Use CGAL \code{triangulate_refine_and_fair_hole()} (\code{TRUE})
#'     or \code{triangulate_and_refine_hole()} (\code{FALSE})?
#' @param maxNumHoles \code{integer}. Maximum number of holes to be filled. May be 0.
#' @param normals Boolean. Whether to return vertex normals.
#' @param verbose Boolean. Whether to print out messages about mesh processing.
#' @returns \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#' @details See \url{https://www.cgal.org/2025/06/13/autorefine-and-snap/}.
#'     If faces are not already triangle, the mesh is first triangulated.
#'
#' @examples
#' library(MeshUtils)
#' mesh      <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' mesh_fill <- fillBoundaryHoles(mesh)
#'
#' @export
fillBoundaryHoles <- function(
  x, fairHole = TRUE, maxNumHoles=10L, normals = FALSE, verbose = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(fairHole))
  stopifnot(isStrictPositiveInteger(maxNumHoles))
  stopifnot(isBoolean(normals))
  stopifnot(isBoolean(verbose))
  meshCPP <- fromR(x)
  mesh    <- fillBoundaryHoles_cpp(meshCPP, fairHole, maxNumHoles, normals, verbose)
  fromCPP(mesh)
}
