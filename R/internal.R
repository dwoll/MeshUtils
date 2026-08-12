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

#' @noRd
isFalsy <- function(x){
  isFALSE(x) || is.null(x) || is.na(x)
}

#' @noRd
makeTriangle <- function(vertices, indices) {
  vertices[indices, ]
}

#' @noRd
isAtomicVector <- function(x) {
  is.atomic(x) && is.vector(x)
}

#' @noRd
isPositiveNumber <- function(x) {
	is.numeric(x) && (length(x) == 1L) && (x > 0) && !is.na(x)
}

#' @noRd
isNonNegativeNumber <- function(x){
	is.numeric(x) && (length(x) == 1L) && (x >= 0) && !is.na(x)
}

#' @noRd
isPositiveInteger <- function(x){
	is.numeric(x) && (length(x) == 1L) && !is.na(x) && (floor(x) == x) && (x >= 0)
}

#' @noRd
isStrictPositiveInteger <- function(x){
	isPositiveInteger(x) && (x > 0)
}

#' @noRd
isBoolean <- function(x){
	is.logical(x) && (length(x) == 1L) && !is.na(x)
}

#' @noRd
isString <- function(x){
  is.character(x) && (length(x) == 1L) && !is.na(x)
}

#' @noRd
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
#' @noRd
fromR <- function(x) {
  vertices <- t(x[["vertices"]])
  faces    <- lapply(seq_len(nrow(x[["faces"]])),
                     function(i) { x[["faces"]][i, ] - 1L })
  list("vertices"=vertices, "faces"=faces)
}

## convert CPP mesh to format required in R
#' @importFrom utils hasName
#' @noRd
fromCPP <- function(x) {
  x[["vertices"]] <- t(x[["vertices"]])
  if(hasName(x, "edges")) {
    edgesDF                 <- x[["edges"]]
    x[["edgesDF"]]          <- edgesDF
    x[["edges"]]            <- as.matrix(edgesDF[, c("i1", "i2")])
    exteriorEdges           <- as.matrix(subset(edgesDF, exterior)[, c("i1", "i2")])
    x[["exteriorEdges"]]    <- exteriorEdges
    x[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
  }

  if(hasName(x, "normals")) {
    x[["normals"]] <- t(x[["normals"]])
  }

  attr(x, "toRGL") <- FALSE
  if(is.matrix(x[["faces"]])) {
    x[["faces"]] <- t(x[["faces"]])
    if(ncol(x[["faces"]]) %in% c(3L, 4L)) {
      attr(x, "toRGL") <- ncol(x[["faces"]])
    }
  } else {
    n_verts   <- lengths(x[["faces"]])
    n_verts_u <- unique(n_verts)
    if(length(n_verts_u) == 1L) {
      x[["faces"]]     <- do.call(rbind, x[["faces"]])
      attr(x, "toRGL") <- n_verts_u
    } else {
      if(all(n_verts_u %in% c(3L, 4L))) {
        attr(x, "toRGL") <- 34L
      }
    }
  }

  class(x) <- "CGALmesh"
  x
}

#' @noRd
checkMesh <- function(vertices, faces, aslist) {
  if(!is.matrix(vertices) || (ncol(vertices) != 3L)) {
    stop("The `vertices` argument must be a matrix with three columns.")
  }

  stopifnot(is.numeric(vertices))
  storage.mode(vertices) <- "double"

  if(anyNA(vertices)) {
    stop("Found missing values in `vertices`.")
  }

  homogeneousFaces <- FALSE
  isTriangle       <- FALSE
  toRGL            <- FALSE
  if(is.matrix(faces)) {
    if(ncol(faces) < 3L) {
      stop("Faces must be given by at least three indices.")
    }
    storage.mode(faces) <- "integer"
    if(anyNA(faces)) {
      stop("Found missing values in `faces`.")
    }
    if(any(faces < 1L)) {
      stop("Faces cannot contain indices lower than 1.")
    }
    if(any(faces > nrow(vertices))) {
      stop("Faces cannot contain indices higher than the number of vertices.")
    }

    homogeneousFaces <- ncol(faces)
    if(homogeneousFaces %in% c(3L, 4L)) {
      isTriangle <- homogeneousFaces == 3L
      toRGL <- homogeneousFaces
    }

    if(aslist) {
      faces <- lapply(seq_len(nrow(faces)), function(i) faces[i, ] - 1L)
    } else {
      faces <- t(faces - 1L)
    }
  } else if(is.list(faces)) {
    check <- all(vapply(faces, isAtomicVector, logical(1L)))
    if(!check) {
      stop("The `faces` argument must be a list of integer vectors.")
    }

    check <- any(vapply(faces, anyNA, logical(1L)))
    if(check) {
      stop("Found missing values in `faces`.")
    }

    faces <- lapply(faces, function(x) as.integer(x) - 1L)
    sizes <- lengths(faces)
    if(any(sizes < 3L)) {
      stop("Faces must be given by at least three indices.")
    }

    check <- any(vapply(faces, function(f) {
              any(f < 0L) || any(f >= nrow(vertices))
            }, logical(1L)))
    if(check) {
      stop(
          "Faces cannot contain indices lower than 1 or higher than the ",
          "number of vertices."
      )
    }
    usizes <- length(unique(sizes))
    if(usizes == 1L) {
      homogeneousFaces <- sizes[1L]
      isTriangle <- homogeneousFaces == 3L
      if(homogeneousFaces %in% c(3L, 4L)) {
        toRGL <- homogeneousFaces
      }
    } else if((usizes == 2L) && all(sizes %in% c(3L, 4L))) {
      toRGL <- 34L
    }
  } else {
    stop("The `faces` argument must be a list or a matrix.")
  }
  list("vertices"        =t(vertices),
       "faces"           =faces,
       "homogeneousFaces"=homogeneousFaces,
       "isTriangle"      =isTriangle,
       "toRGL"           =toRGL)
}
