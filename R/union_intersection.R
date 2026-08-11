## ----------------------------------------------------------------------- //
## Code adapted from packages
## https://github.com/stla/Boov/
## https://github.com/stla/PolygonSoup/
## https://github.com/stla/cgalMeshes/
## developed and copyright by
## Stéphane Laurent <laurent_step@outlook.fr>
## License: GPL-3
## ----------------------------------------------------------------------- //

#' @title Meshes intersection
#' @description Computes the intersection of the given meshes.
#'
#' @param meshes a list of meshes, each being either a
#'   \strong{rgl} mesh, or as a list with (at least) two fields:
#'   \code{vertices} and \code{faces}; the \code{vertices} matrix can
#'   be a numeric matrix or a matrix of \code{bigq} rational numbers
#'   (from the \strong{gmp} package)
#' @param clean Boolean, whether to clean the meshes (merging
#'   duplicated vertices, duplicated faces, removing isolated vertices);
#'   set to \code{FALSE} if you are sure your meshes are clean, to
#'   gain some speed
#' @param normals Boolean, whether to return the per-vertex normals of the
#'   output mesh
#'
#' @return A triangle mesh given as a list with fields \code{vertices},
#'   \code{faces}, \code{edges}, \code{exteriorEdges}, \code{gmpvertices}
#'   if using \strong{gmp} meshes, and \code{normals} if \code{normals=TRUE}.
#'
#' @export
#'
#' @examples
#' library(Boov)
#' library(rgl)
#'
#' # mesh one: truncated icosahedron; we triangulate it for plotting
#' library(PolygonSoup)
#' mesh1 <- makeMesh(
#'   mesh = truncatedIcosahedron,
#'   triangulate = TRUE, normals = FALSE
#' )
#'
#' # mesh two: a cube
#' mesh2 <- translate3d( # (from the rgl package)
#'   cube3d(), 2, 0, 0
#' )
#'
#' # compute the intersection
#' inter <- MeshesIntersection(list(mesh1, mesh2))
#'
#' # plot
#' rglmesh1 <- toRGL(mesh1)
#' rglinter <- toRGL(inter)
#' open3d(windowRect = c(50, 50, 562, 562))
#' shade3d(rglmesh1, color = "yellow", alpha = 0.2)
#' shade3d(mesh2, color = "cyan", alpha = 0.2)
#' shade3d(rglinter, color = "red")
#' plotEdges(
#'   vertices = inter[["vertices"]], edges = inter[["exteriorEdges"]],
#'   edgesAsTubes = FALSE, lwd = 3, verticesAsSpheres = FALSE
#' )
#'
#' # other example, with 'gmp' rational numbers ####
#' library(Boov)
#' library(gmp)
#' library(rgl)
#'
#' cube <- cube3d()
#'
#' rglmesh1 <- cube
#' mesh1 <-
#'   list(vertices = t(cube[["vb"]][-4L, ]), faces = t(cube[["ib"]]))
#' mesh1[["vertices"]] <- as.bigq(mesh1[["vertices"]])
#'
#' rotMatrix <- t(cbind( # pi/3 around a great diagonal
#'   as.bigq(c(2, -1, 2), c(3, 3, 3)),
#'   as.bigq(c(2, 2, -1), c(3, 3, 3)),
#'   as.bigq(c(-1, 2, 2), c(3, 3, 3))
#' ))
#' mesh2 <-
#'   list(vertices = t(cube[["vb"]][-4L, ]), faces = t(cube[["ib"]]))
#' mesh2[["vertices"]] <- as.bigq(mesh2[["vertices"]]) %*% rotMatrix
#' rglmesh2 <- rotate3d(cube, pi/3, 1, 1, 1)
#'
#' inter <- MeshesIntersection(list(mesh1, mesh2))
#' # perfect vertices:
#' inter[["gmpVertices"]]
#' rglinter <- toRGL(inter)
#'
#' open3d(windowRect = c(50, 50, 562, 562), zoom = 0.9)
#' bg3d("#363940")
#' shade3d(rglmesh1, color = "yellow", alpha = 0.2)
#' shade3d(rglmesh2, color = "orange", alpha = 0.2)
#' shade3d(rglinter, color = "hotpink")
#' plotEdges(
#'   inter[["vertices"]], inter[["exteriorEdges"]],
#'   only = inter[["exteriorVertices"]],
#'   color = "firebrick",
#'   tubesRadius = 0.05, spheresRadius = 0.07
#' )
meshIntersection <- function(meshes, clean=TRUE, normals=FALSE) {
  stopifnot(is.list(meshes))
  stopifnot(length(meshes) >= 2L)
  checkMeshes <- lapply(meshes, function(mesh) {
    if(inherits(mesh, "mesh3d")) {
      vft  <- getVFT(mesh, beforeCheck = TRUE)
      mesh <- vft[["rmesh"]]
    }
    checkMesh(mesh[["vertices"]], mesh[["faces"]], aslist = TRUE)
  })
  areTriangle <- vapply(checkMeshes, `[[`, logical(1L), "isTriangle")
  triangulate <- !areTriangle
  meshes      <- lapply(checkMeshes, `[`, c("vertices", "faces"))
  inter       <- Intersection_EK(meshes, clean, normals, triangulate)
  fromCPP(inter)
}

#' @title Meshes difference
#' @description Computes the difference between two meshes.
#'
#' @param mesh1,mesh2 two meshes, each being given as either a
#'   \strong{rgl} mesh, or a list with (at least) two fields:
#'   \code{vertices} and \code{faces}; the \code{vertices} matrix can be
#'   a numeric matrix or a matrix of \code{bigq} rational numbers (from the
#'   \strong{gmp} package)
#' @param clean Boolean, whether to clean the meshes (merging duplicated
#'   vertices, duplicated faces, removing isolated vertices); set to
#'   \code{FALSE} if you know your meshes are clean
#' @param normals Boolean, whether to return the per-vertex normals of the
#'   output mesh
#'
#' @return A triangle mesh given as a list with fields \code{vertices},
#'   \code{faces}, \code{edges}, \code{exteriorEdges}, \code{gmpvertices}
#'   if using \strong{gmp} meshes, and \code{normals} if \code{normals=TRUE}.
#'
#' @export
#'
#' @examples
#' library(Boov)
#' library(rgl)
#'
#' # mesh one: a cube
#' mesh1 <- cube3d() # (from the rgl package)
#'
#' # mesh two: another cube
#' mesh2 <- translate3d( # (from the rgl package)
#'   cube3d(), 1, 1, 0
#' )
#'
#' # compute the difference
#' differ <- MeshesDifference(mesh1, mesh2)
#'
#' # plot
#' rgldiffer <- toRGL(differ)
#' open3d(windowRect = c(50, 50, 562, 562))
#' shade3d(mesh1, color = "yellow", alpha = 0.2)
#' shade3d(mesh2, color = "cyan", alpha = 0.2)
#' shade3d(rgldiffer, color = "red")
#' plotEdges(
#'   vertices = differ[["vertices"]], edges = differ[["exteriorEdges"]],
#'   edgesAsTubes = TRUE, verticesAsSpheres = TRUE
#' )
meshDifference <- function(mesh1, mesh2, clean=TRUE, normals=FALSE) {
  stopifnot(is.list(mesh1), is.list(mesh2))

  if(inherits(mesh1, "mesh3d")) {
    vft   <- getVFT(mesh1, beforeCheck = TRUE)
    mesh1 <- vft[["rmesh"]]
  }
  checkMesh1   <- checkMesh(mesh1[["vertices"]], mesh1[["faces"]], aslist = TRUE)
  triangulate1 <- !checkMesh1[["isTriangle"]]

  if(inherits(mesh2, "mesh3d")) {
    vft   <- getVFT(mesh2, beforeCheck = TRUE)
    mesh2 <- vft[["rmesh"]]
  }
  checkMesh2   <- checkMesh(mesh2[["vertices"]], mesh2[["faces"]], aslist = TRUE)
  triangulate2 <- !checkMesh2[["isTriangle"]]

  mesh1  <- checkMesh1[c("vertices", "faces")]
  mesh2  <- checkMesh2[c("vertices", "faces")]
  differ <- Difference_EK(mesh1, mesh2, clean, normals, triangulate1, triangulate2)
  fromCPP(differ)
}

#' @title Meshes union
#' @description Computes the union of the given meshes.
#'
#' @param meshes a list of meshes, each being either a
#'   \strong{rgl} mesh, or as a list with (at least) two fields:
#'   \code{vertices} and \code{faces}; the \code{vertices} matrix can
#'   be a numeric matrix or a matrix of \code{bigq} rational numbers
#'   (from the \strong{gmp} package)
#' @param clean Boolean, whether to clean the meshes (merging
#'   duplicated vertices, duplicated faces, removing isolated vertices);
#'   set to \code{FALSE} if you are sure your meshes are clean, to
#'   gain some speed
#' @param normals Boolean, whether to return the per-vertex normals of the
#'   output mesh
#'
#' @return A triangle mesh given as a list with fields \code{vertices},
#'   \code{faces}, \code{edges}, \code{exteriorEdges}, \code{gmpvertices}
#'   if using \strong{gmp} meshes, and \code{normals} if \code{normals=TRUE}.
#'
#' @export
#'
#' @examples
#' library(Boov)
#' library(rgl)
#'
#' # mesh one: a cube
#' mesh1 <- cube3d() # (from the rgl package)
#'
#' # mesh two: another cube
#' mesh2 <- translate3d( # (from the rgl package)
#'   cube3d(), 1, 1, 1
#' )
#'
#' # compute the union
#' umesh <- MeshesUnion(list(mesh1, mesh2))
#'
#' # plot
#' rglumesh <- toRGL(umesh)
#' open3d(windowRect = c(50, 50, 562, 562))
#' shade3d(rglumesh, color = "red")
#' plotEdges(
#'   vertices = umesh[["vertices"]], edges = umesh[["exteriorEdges"]],
#'   edgesAsTubes = TRUE, verticesAsSpheres = TRUE
#' )
meshUnion <- function(meshes, clean = TRUE, normals = FALSE) {
  stopifnot(is.list(meshes))
  stopifnot(length(meshes) >= 2L)

  checkMeshes <- lapply(meshes, function(mesh) {
    if(inherits(mesh, "mesh3d")) {
      vft  <- getVFT(mesh, beforeCheck = TRUE)
      mesh <- vft[["rmesh"]]
    }
    checkMesh(mesh[["vertices"]], mesh[["faces"]], aslist = TRUE)
  })

  areTriangle <- vapply(checkMeshes, `[[`, logical(1L), "isTriangle")
  triangulate <- !areTriangle
  meshes      <- lapply(checkMeshes, `[`, c("vertices", "faces"))
  umesh       <- Union_EK(meshes, clean, normals, triangulate)
  fromCPP(umesh)
}
