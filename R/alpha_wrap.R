## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/Boov/
## https://github.com/stla/PolygonSoup/
## https://github.com/stla/cgalMeshes/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title 3D alpha wrapping
#' @description Reconstruction of a surface from a cloud of 3D points by
#'   alpha wrapping.
#'
#' @param points numeric matrix which stores the points, one point per row
#' @param alpha_rel relative alpha parameter; the actual alpha parameter (see
#'   details) is defined as the length of the diagonal of the bounding box of
#'   the points cloud divided by the relative alpha parameter
#' @param offset_rel relative offset; the actual offset parameter (see details)
#'   is defined as the length of the diagonal of the bounding box of the
#'   points cloud divided by the relative offset parameter
#'
#' @return A \code{CGALmesh} object.
#'
#' @details See \href{https://doc.cgal.org/latest/Alpha_wrap_3/index.html}{3D Alpha Wrapping}
#'   for details. The smallest alpha parameter, the smallest triangles in the
#'   output mesh. The offset is the distance from the input points to the output
#'   mesh.
#'
#' @export
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#' # take the diplodocus mesh
#' off <- system.file("extdata", "diplodocus.off", package="MeshUtils")
#' mesh_raw <- readMeshFile(off)
#' diplodocusMesh <- makeMesh(mesh_raw$vertices, mesh_raw$faces, normals=TRUE)
#' # reconstruct the mesh from its vertices
#' pts <- mesh_raw$vertices
#' \donttest{wrapMesh <- alphaWrap(pts, 70, 3000)
#' # plot
#' diplodocus_rgl <- toRGL(diplodocusMesh)
#' wrap_rgl       <- toRGL(wrapMesh)
#' open3d(windowRect = 50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(20, 0, zoom=0.85)
#' shade3d(diplodocus_rgl, color="forestgreen")
#' next3d()
#' view3d(20, 0, zoom=0.85)
#' shade3d(wrap_rgl, color="forestgreen")}
alphaWrap <- function(points, alphaRel, offsetRel){
  stopifnot(isPositiveNumber(alphaRel))
  stopifnot(isPositiveNumber(offsetRel))
  if(!is.matrix(points) || !is.numeric(points)){
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  if(ncol(points) != 3L){
    stop("The `points` matrix must have three columns.", call. = TRUE)
  }
  if(nrow(points) <= 3L){
    stop("Insufficient number of points.", call. = TRUE)
  }
  storage.mode(points) <- "double"
  if(anyNA(points)){
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  aw  <- alphaWrap_cpp(t(points), alphaRel, offsetRel)
  out <- vcgUpdateNormals(
    tmesh3d(
      aw[["vertices"]], aw[["faces"]], normals = NULL,
      homogeneous = FALSE
    )
  )
  out
}
