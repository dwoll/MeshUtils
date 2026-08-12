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
#' @param x Numeric matrix which stores the points, one point per row.
#' @param alphaRel Relative alpha parameter; the actual alpha parameter (see
#'   details) is defined as the length of the diagonal of the bounding box of
#'   the point cloud divided by the relative alpha parameter.
#' @param offsetRel Relative offset; the actual offset parameter (see details)
#'   is defined as the length of the diagonal of the bounding box of the
#'   point cloud divided by the relative offset parameter.
#' @param out Character to indicate output mesh format.
#'
#' @returns A \code{CGALmesh} object.
#'
#' @details See \href{https://doc.cgal.org/latest/Alpha_wrap_3/index.html}{3D Alpha Wrapping}
#'   for details. The smallest alpha parameter, the smallest triangles in the
#'   output mesh. The offset is the distance from the input points to the output
#'   mesh.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#' # take the diplodocus mesh
#' off  <- system.file("extdata", "diplodocus.off", package="MeshUtils")
#' vf   <- readMeshFile(off)
#' mesh <- makeMesh(mesh=vf, normals=TRUE)
#'
#' # reconstruct the mesh from the points
#' mesh_alwrap_rgl <- alphaWrap(mesh_raw[["vertices"]],
#'                              alphaRel =70,
#'                              offsetRel=3000,
#'                              out      ="rgl")
#'
#' # plot
#' mesh_rgl <- toRGL(mesh)
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(20, 0, zoom=0.85)
#' shade3d(mesh_rgl, color="forestgreen")
#' next3d()
#' view3d(20, 0, zoom=0.85)
#' shade3d(mesh_alwrap_rgl, color="forestgreen")
#'
#' @export
#' @importFrom Rvcg vcgUpdateNormals
#' @importFrom rgl tmesh3d
alphaWrap <- function(x, alphaRel, offsetRel, out=c("CGALmesh", "rgl")) {
  out <- match.arg(out)
  stopifnot(isPositiveNumber(alphaRel))
  stopifnot(isPositiveNumber(offsetRel))
  if(!is.matrix(x) || !is.numeric(x)){
    stop("The `x` argument must be a numeric matrix.", call. = TRUE)
  }
  if(ncol(x) != 3L) {
    stop("The `points` matrix must have three columns.", call. = TRUE)
  }
  if(nrow(x) <= 3L) {
    stop("Insufficient number of points.", call. = TRUE)
  }
  storage.mode(x) <- "double"
  if(anyNA(x)){
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  alwrap   <- alphaWrap_cpp(t(x), alphaRel, offsetRel)
  mesh_rwn <- vcgUpdateNormals(
    tmesh3d(alwrap[["vertices"]],
            alwrap[["faces"]],
            normals    =NULL,
            homogeneous=FALSE))
  mesh_out <- if(out == "rgl") {
    mesh_rwn
  } else {
    makeMesh(mesh=mesh_rwn)
  }
  mesh_out
}
