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

#' @title Read a mesh file
#' @description Read mesh vertices and faces from a file.
#'
#' @param x Path to the mesh file; supported formats are \code{stl},
#'   \code{ply}, \code{obj} and \code{off}.
#' @param method Either \code{"soup"} when the file is a polygon soup, or \code{"mesh"} when
#'   the file is a valid mesh.
#' @param binary Boolean: For \code{method="mesh"}: Whether input file is binary.
#' @param normals Boolean: Whether to return vertex normals for \code{method="mesh"}.
#'
#' @returns For \code{method="soup"}: A list with two components: \code{vertices},
#'   a numeric matrix with three
#'   columns, and \code{faces}, either a list of integer vectors or, in the
#'   case if all faces have the same number of sides, an integer matrix.
#'   For \code{method="mesh"}: A list of class \code{CGALmesh} giving the vertices,
#'   the edges, the faces of the mesh, the exterior edges, the exterior vertices and
#'   optionally the normals.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#' ply  <- system.file("extdata", "dataHeart3.ply", package="MeshUtils")
#' vf   <- readMeshFile(ply, method="soup")
#' mesh <- makeMesh(mesh=vf, normals=TRUE)
#'
#' mesh_rgl <- toRGL(mesh)
#' open3d(windowRect=c(50, 50, 562, 562))
#' view3d(0, 0, zoom=0.8)
#' shade3d(mesh_rgl, color="palevioletred")
#'
#' @export
readMeshFile <- function(x, method=c("soup", "mesh"), binary=FALSE, normals=FALSE) {
  stopifnot(isString(x))
  stopifnot(isBoolean(binary))
  method <- match.arg(method)
  if(!file.exists(x)) {
    stop("File not found.")
  }
  if(method == "soup") {
    mesh   <- readFileSoup_cpp(x, binary)
    faces  <- mesh[["faces"]]
    usizes <- length(unique(lengths(faces)))
    if(usizes == 1L) {
      mesh[["faces"]] <- do.call(rbind, faces)
    }
    mesh
  } else {
    stopifnot(isBoolean(normals))
    mesh_cpp <- readFileMesh_cpp(x, binary, normals)
    fromCPP(mesh_cpp)
  }
}

#' @title Export mesh to a file
#' @description Export a mesh to a file.
#'
#' @param x A mesh given either as a list containing (at least) the fields
#'   \code{vertices} and \code{faces}, otherwise a \strong{rgl} mesh
#'   (i.e. a \code{\link[rgl]{mesh3d}} object).
#' @param filename Name of the file to be written, with extension \code{stl},
#'   \code{ply}, \code{obj} or \code{off}
#' @param precision Positive integer, number of decimal digits for the vertices.
#' @param binary Boolean, whether to write a binary file or an ASCII file.
#'
#' @returns No value, just generates the file.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @export
writeMeshFile <- function(x, filename, precision = 17L, binary = FALSE) {
  stopifnot(isString(filename))
  stopifnot(isPositiveInteger(precision))
  stopifnot(isBoolean(binary))
  if(inherits(x, "mesh3d")){
    vft <- getVFT(x, beforeCheck = TRUE)
    x   <- vft[["rmesh"]]
  }
  vertices    <- x[["vertices"]]
  faces       <- x[["faces"]]
  checkedMesh <- checkMesh(vertices, faces, aslist=TRUE)
  vertices    <- checkedMesh[["vertices"]]
  faces       <- checkedMesh[["faces"]]
  writeFile_cpp(filename,
                binary,
                as.integer(precision),
                vertices,
                faces)
  invisible(NULL)
}
