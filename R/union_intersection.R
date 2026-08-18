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

#' @title Meshes intersection
#' @description Computes the intersection of the given meshes.
#'
#' @param x A list of meshes, each being either a \code{mesh3d} object
#'   from package \strong{rgl}, or as a list with (at least) two fields:
#'   \code{vertices} and \code{faces}, such as a \code{CGALmesh} object.
#' @param repairSoup Boolean, whether to clean the meshes (merging
#'   duplicated vertices, duplicated faces, removing isolated vertices).
#'   Set to \code{FALSE} if you are sure your meshes are clean, to
#'   gain some speed.
#' @param normals Boolean, whether to return the per-vertex normals of the
#'   output mesh.
#'
#' @returns A \code{CGALmesh} object.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' # mesh one: truncated icosahedron; we triangulate it for plotting
#' mesh1 <- makeMesh(mesh       =dataTruncIcosahedron,
#'                   triangulate=TRUE,
#'                   normals    = FALSE)
#'
#' # mesh two: a cube
#' mesh2_rgl <- translate3d(cube3d(), 2, 0, 0)
#' mesh2     <- makeMesh(mesh       =mesh2_rgl,
#'                       triangulate=TRUE,
#'                       normals    =FALSE)
#'
#' # compute the intersection
#' mesh_i <- boolIntersection(list(mesh1, mesh2))
#'
#' # plot
#' mesh1_rgl  <- toRGL(mesh1)
#' mesh_i_rgl <- toRGL(mesh_i)
#' open3d(windowRect=c(50, 50, 562, 562))
#' shade3d(mesh1_rgl,  color="yellow", alpha=0.2)
#' shade3d(mesh2_rgl,  color="cyan",   alpha=0.2)
#' shade3d(mesh_i_rgl, color="red")
#' plotEdges(vertices         =mesh_i[["vertices"]],
#'           edges            =mesh_i[["exteriorEdges"]],
#'           edgesAsTubes     =FALSE,
#'           lwd              =3,
#'           verticesAsSpheres=FALSE)
#'
#' @export
boolIntersection <- function(x, repairSoup=TRUE, normals=FALSE) {
  stopifnot(is.list(x))
  stopifnot(length(x) >= 2L)
  stopifnot(isBoolean(repairSoup))
  stopifnot(isBoolean(normals))
  checkMeshes <- lapply(x, function(mesh) {
    if(inherits(mesh, "mesh3d")) {
      vft  <- getVFT(mesh, beforeCheck = TRUE)
      mesh <- vft[["rmesh"]]
    }
    checkMesh(mesh[["vertices"]], mesh[["faces"]], aslist = TRUE)
  })
  areTriangle <- vapply(checkMeshes, `[[`, logical(1L), "isTriangle")
  triangulate <- !areTriangle
  meshes      <- lapply(checkMeshes, `[`, c("vertices", "faces"))
  inter       <- intersectionEK_cpp(meshes, triangulate, repairSoup, normals)
  fromCPP(inter)
}

#' @title Mesh difference
#' @description Computes the difference between two meshes.
#'
#' @param mesh1 A mesh, either being given a \code{mesh3d} object
#'   from package \strong{rgl}, or by a list with (at least) two fields:
#'   \code{vertices} and \code{faces}, such as a \code{CGALmesh} object.
#' @param mesh2 A mesh, either being given a \code{mesh3d} object
#'   from package \strong{rgl}, or by a list with (at least) two fields:
#'   \code{vertices} and \code{faces}, such as a \code{CGALmesh} object.
#' @param repairSoup Boolean, whether to clean the meshes (merging duplicated
#'   vertices, duplicated faces, removing isolated vertices). Set to
#'   \code{FALSE} if you know your meshes are clean.
#' @param normals Boolean, whether to return the per-vertex normals of the
#'   output mesh.
#'
#' @returns A \code{CGALmesh} object.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' # mesh one: a cube
#' mesh1_rgl <- cube3d() # (from the rgl package)
#'
#' # mesh two: another cube
#' mesh2_rgl <- translate3d(cube3d(), 1, 1, 0)
#'
#' # compute the difference
#' mesh_d <- boolDifference(mesh1_rgl, mesh2_rgl)
#'
#' # plot
#' mesh_d_rgl <- toRGL(mesh_d)
#' open3d(windowRect=c(50, 50, 562, 562))
#' shade3d(mesh1_rgl,  color="yellow", alpha=0.2)
#' shade3d(mesh2_rgl,  color="cyan",   alpha=0.2)
#' shade3d(mesh_d_rgl, color="red")
#' plotEdges(vertices         =mesh_d[["vertices"]],
#'           edges            =mesh_d[["exteriorEdges"]],
#'           edgesAsTubes     =TRUE,
#'           verticesAsSpheres=TRUE)
#'
#' @export
boolDifference <- function(mesh1, mesh2, repairSoup=TRUE, normals=FALSE) {
  stopifnot(is.list(mesh1), is.list(mesh2))
  stopifnot(isBoolean(repairSoup))
  stopifnot(isBoolean(normals))

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
  differ <- differenceEK_cpp(mesh1, mesh2, triangulate1, triangulate2, repairSoup, normals)
  fromCPP(differ)
}

#' @title Meshes union
#' @description Computes the union of the given meshes.
#'
#' @param x A list of meshes, each being either a \code{mesh3d} object
#'   from package \strong{rgl}, or as a list with (at least) two fields:
#'   \code{vertices} and \code{faces}, such as a \code{CGALmesh} object.
#' @param repairSoup Boolean, whether to clean the meshes (merging
#'   duplicated vertices, duplicated faces, removing isolated vertices).
#'   Set to \code{FALSE} if you are sure your meshes are clean, to
#'   gain some speed.
#' @param normals Boolean, whether to return the per-vertex normals of the
#'   output mesh.
#'
#' @returns A \code{CGALmesh} object.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' # mesh one: a cube
#' mesh1_rgl <- cube3d() # (from the rgl package)
#'
#' # mesh two: another cube
#' mesh2_rgl <- translate3d(cube3d(), 1, 1, 1)
#'
#' # compute the union
#' mesh_u <- boolUnion(list(mesh1_rgl, mesh2_rgl))
#'
#' # plot
#' mesh_u_rgl <- toRGL(mesh_u)
#' open3d(windowRect=c(50, 50, 562, 562))
#' shade3d(mesh_u_rgl, color="red")
#' plotEdges(vertices         =mesh_u[["vertices"]],
#'           edges            =mesh_u[["exteriorEdges"]],
#'           edgesAsTubes     =TRUE,
#'           verticesAsSpheres=TRUE)
#'
#' @export
boolUnion <- function(x, repairSoup = TRUE, normals = FALSE) {
  stopifnot(is.list(x))
  stopifnot(length(x) >= 2L)
  stopifnot(isBoolean(repairSoup))
  stopifnot(isBoolean(normals))

  checkMeshes <- lapply(x, function(mesh) {
    if(inherits(mesh, "mesh3d")) {
      vft  <- getVFT(mesh, beforeCheck = TRUE)
      mesh <- vft[["rmesh"]]
    }
    checkMesh(mesh[["vertices"]], mesh[["faces"]], aslist = TRUE)
  })

  areTriangle <- vapply(checkMeshes, `[[`, logical(1L), "isTriangle")
  triangulate <- !areTriangle
  meshes      <- lapply(checkMeshes, `[`, c("vertices", "faces"))
  umesh       <- unionEK_cpp(meshes, triangulate, repairSoup, normals)
  fromCPP(umesh)
}
