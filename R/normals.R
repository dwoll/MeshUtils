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

#' @title Normals for a point cloud
#' @description Returns a function which computes some normals for a 3D point
#'   cloud.
#'
#' @param x Integer, number of neighbors used to compute the normals.
#' @param method One of \code{"pca"} or \code{"jet"}.
#'
#' @returns A function which takes just one argument: a numeric matrix with
#'   three columns, each row represents a point, and the function returns a
#'   matrix of the same size as the input matrix, with each row giving one
#'   unit normal for the point.
#'
#' @note The \code{getSomeNormals} function is intended to be used in the
#'   \code{\link{reconstructPoisson}} function. If you want to use it for
#'   another purpose, be careful because the function it returns does not
#'   check the matrix it takes as argument.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#' fun  <- getSomeNormals(6)
#' mesh <- makeMesh(HopfTorus$vertices,
#'                  HopfTorus$faces)
#' mesh_rgl <- reconstructPoisson(mesh, fun, out="rgl")
#' open3d()
#' wire3d(mesh_rgl)
#'
#' @export
getSomeNormals <- function(x, method = c("PCA", "Jet")) {
  method <- match.arg(tolower(method), choices=c("pca", "jet"))
  x <- as.integer(x)
  if(x < 2L) {
    stop("There must be at least two neighbors.", call. = TRUE)
  }
  fun <- if(method == "pca") {
    function(points) { pca_normals_cpp(t(points), x) }
  } else{
    function(points) { jet_normals_cpp(t(points), x) }
  }
  class(fun) <- "CGALnormalsFunc"
  fun
}
