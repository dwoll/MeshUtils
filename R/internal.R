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

isFalsy <- function(x){
  isFALSE(x) || is.null(x) || is.na(x)
}

makeTriangle <- function(vertices, indices) {
  vertices[indices, ]
}

isAtomicVector <- function(x) {
  is.atomic(x) && is.vector(x)
}

isPositiveNumber <- function(x) {
	is.numeric(x) && (length(x) == 1L) && (x > 0) && !is.na(x)
}

isNonNegativeNumber <- function(x){
	is.numeric(x) && (length(x) == 1L) && (x >= 0) && !is.na(x)
}

isPositiveInteger <- function(x){
	is.numeric(x) && (length(x) == 1L) && !is.na(x) && (floor(x) == x) && (x >= 0)
}

isStrictPositiveInteger <- function(x){
	isPositiveInteger(x) && (x > 0)
}

isBoolean <- function(x){
	is.logical(x) && (length(x) == 1L) && !is.na(x)
}

isString <- function(x){
  is.character(x) && (length(x) == 1L) && !is.na(x)
}

getVFT <- function(x, beforeCheck = FALSE) {
  transposed <- !beforeCheck
  i0 <- as.integer(transposed)
  if(inherits(x, "mesh3d")) {
    triangles <- x[["it"]]
    if(!is.null(triangles)) {
      triangles <- lapply(seq_len(ncol(triangles)), function(i) { triangles[, i] - i0 })
    }
    quads <- x[["ib"]]
    isTriangle <- is.null(quads)
    if(!isTriangle) {
      quads <- lapply(seq_len(ncol(quads)), function(i) { quads[, i] - i0 })
    }
    faces <- c(triangles, quads)
    vertices <- x[["vb"]][-4L, ]
    if(!transposed) {
      vertices <- t(vertices)
    }
    rmesh <- list("vertices" = vertices, "faces" = faces)
  } else if(inherits(x, "CGALmesh")) {
    isTriangle <- attr(x, "toRGL") == 3L
    vertices <- x[["vertices"]]
    if(transposed) {
      vertices <- t(vertices)
    }
    faces <- x[["faces"]]
    if(is.matrix(faces)) {
      faces <- lapply(seq_len(nrow(faces)), function(i) { faces[i, ] - i0 })
    } else if(!beforeCheck) {
      faces <- lapply(faces, function(face) { face - 1L })
    }
    rmesh <- list("vertices" = vertices, "faces" = faces)
  } else if(is.list(x)) {
    rmesh <- checkMesh(x[["vertices"]], x[["faces"]], aslist = TRUE)
    isTriangle <- rmesh[["isTriangle"]]
    if(beforeCheck) {
      rmesh[["vertices"]] <- t(rmesh[["vertices"]])
      rmesh[["faces"]] <- lapply(rmesh[["faces"]], function(face) { face + 1L })
    }
  } else {
    stop("Invalid `x` argument.", call. = FALSE)
  }
  list("rmesh" = rmesh, "isTriangle" = isTriangle)
}

## convert R mesh to format required for CPP
fromR <- function(x) {
  vertices <- t(x[["vertices"]])
  faces    <- lapply(seq_len(nrow(x[["faces"]])),
                     function(i) { x[["faces"]][i, ] - 1L })
  list("vertices"=vertices, "faces"=faces)
}

## convert CPP mesh to format required in R
## TODO check whether t(mesh[["faces"]]) is really necessary - may be list of 3-vectors anyway
#' @importFrom utils hasName
#' @noRd
fromCPP <- function(x) {
  x[["vertices"]]         <- t(x[["vertices"]])
  x[["faces"]]            <- t(x[["faces"]])
  edgesDF                 <- x[["edges"]]
  x[["edgesDF"]]          <- edgesDF
  x[["edges"]]            <- as.matrix(edgesDF[, c("i1", "i2")])
  exteriorEdges           <- as.matrix(subset(edgesDF, exterior)[, c("i1", "i2")])
  x[["exteriorEdges"]]    <- exteriorEdges
  x[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
  if(hasName(x, "normals")) {
    x[["normals"]] <- t(x[["normals"]])
  }

  attr(x, "toRGL") <- 3L
  class(x) <- "CGALmesh"
  x
}
