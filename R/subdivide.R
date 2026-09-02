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

#' @title Catmull-Clark subdivision and deformation
#' @description Performs the Catmull-Clark subdivision and deformation.
#'   Includes triangulation if mesh is not already triangle.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#'   The mesh must be triangle or be able to be made triangle.
#' @param nIter Positive \code{integer}: Number of iterations.
#' @param normals Boolean: Whether to return vertex normals.
#' @returns A \code{CGALmesh} object.
#'
#' @details See \url{https://doc.cgal.org/latest/Subdivision_method_3/} for details.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' mesh        <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' mesh_rgl    <- toRGL(mesh)
#' mesh_sd     <- subdivideCatmullClark(mesh, nIter=2)
#' mesh_sd_rgl <- toRGL(mesh_sd)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_rgl)
#' next3d()
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_sd_rgl)
#'
#' @export
subdivideCatmullClark <- function(x, nIter = 1, normals = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			     " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isBoolean(normals))
  meshCPP <- fromR(x)
  meshOut <- subdivideCatmullClark_cpp(meshCPP, as.integer(nIter), normals)
  fromCPP(meshOut)
}

#' @title Doo-Sabin subdivision and deformation
#' @description Performs the Doo-Sabin subdivision and deformation.
#'   Includes triangulation if mesh is not already triangle.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#'   The mesh must be triangle or be able to be made triangle.
#' @param nIter Positive \code{integer}: Number of iterations.
#' @param normals Boolean: Whether to return vertex normals.
#' @returns A \code{CGALmesh} object.
#' @details See \url{https://doc.cgal.org/latest/Subdivision_method_3/} for details.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' mesh        <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' mesh_rgl    <- toRGL(mesh)
#' mesh_ds     <- subdivideDooSabin(mesh, nIter=2)
#' mesh_ds_rgl <- toRGL(mesh_ds)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_rgl)
#' next3d()
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_ds_rgl)
#'
#' @export
subdivideDooSabin <- function(x, nIter = 1, normals = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			     " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isBoolean(normals))
  meshCPP <- fromR(x)
  meshOut <- subdivideDooSabin_cpp(meshCPP, as.integer(nIter), normals)
  fromCPP(meshOut)
}

#' @title Sqrt3 subdivision and deformation
#' @description Performs the 'Sqrt3' subdivision and deformation.
#'   Includes triangulation if mesh is not already triangle.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#'   The mesh must be triangle or be able to be made triangle.
#' @param nIter Positive \code{integer}: Number of iterations.
#' @param normals Boolean: Whether to return vertex normals.
#' @returns A \code{CGALmesh} object.
#' @details See \url{https://doc.cgal.org/latest/Subdivision_method_3/} for details.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' mesh        <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' mesh_rgl    <- toRGL(mesh)
#' mesh_s3     <- subdivideSqrt3(mesh, nIter=2)
#' mesh_s3_rgl <- toRGL(mesh_s3)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_rgl)
#' next3d()
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_s3_rgl)
#'
#' @export
subdivideSqrt3 <- function(x, nIter = 1, normals = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isBoolean(normals))
  meshCPP <- fromR(x)
  meshOut <- subdivideSqrt3_cpp(meshCPP, as.integer(nIter), normals)
  fromCPP(meshOut)
}
