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
#'   \code{\link[MeshUtils]{reconstructPoisson}} function. If you want to use it for
#'   another purpose, be careful because the function it returns does not
#'   check the matrix it takes as argument.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#' mesh     <- dataHeart1
#' mesh_rgl <- toRGL(mesh)
#' fun      <- getSomeNormals(6)
#' mesh_psr <- reconstructPoisson(mesh[["vertices"]],
#'                                normals=fun,
#'                                smAngle=10,
#'                                smRadius=1.5,
#'                                smDistance=0.3)
#'
#' mesh_psr_rgl <- toRGL(mesh_psr)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' wire3d(mesh_rgl)
#' next3d()
#' wire3d(mesh_psr_rgl)
#'
#' @export
getSomeNormals <- function(x, method = c("PCA", "Jet")) {
  method <- match.arg(tolower(method), choices=c("pca", "jet"))
  x <- as.integer(x)
  if(x < 2L) {
    stop("There must be at least two neighbors.", call. = TRUE)
  }
  fun <- if(method == "pca") {
    function(points) {
      if(!is.matrix(points) || !is.numeric(points)) {
        stop("The `points` argument must be a numeric matrix.", call. = TRUE)
      }
      if(ncol(points) != 3L) {
        stop("The `points` matrix must have three columns.", call. = TRUE)
      }
      if(nrow(points) <= 3L) {
        stop("Insufficient number of points.", call. = TRUE)
      }
      storage.mode(points) <- "double"
      pca_normals_cpp(t(points), x)
    }
  } else{
    function(points) {
      if(!is.matrix(points) || !is.numeric(points)) {
        stop("The `points` argument must be a numeric matrix.", call. = TRUE)
      }
      if(ncol(points) != 3L) {
        stop("The `points` matrix must have three columns.", call. = TRUE)
      }
      if(nrow(points) <= 3L) {
        stop("Insufficient number of points.", call. = TRUE)
      }
      storage.mode(points) <- "double"
      jet_normals_cpp(t(points), x)
    }
  }
  class(fun) <- "CGALnormalsFunc"
  fun
}
