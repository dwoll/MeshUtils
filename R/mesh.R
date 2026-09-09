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

#' @exportS3Method print CGALmesh
print.CGALmesh <- function(x, ...) {
	rgl <- attr(x, "toRGL")
	nv  <- nrow(x[["vertices"]])
	nf  <- if(is.list(x[["faces"]])) { length(x[["faces"]]) } else { nrow(x[["faces"]]) }
	msg <- sprintf("CGALmesh with %d vertices and %d faces.\n", nv, nf)
	cat(msg)
	elr <- formatC(range(x[["edgesDF"]][["length"]]))
	msg <- sprintf("The edge lengths vary from %s to %s.\n", elr[1L], elr[2L])
	cat(msg)
	is  <- if(rgl == 3L) { " is " } else { " is not " }
	msg <- paste0("This mesh", is, "triangle.\n")
	cat(msg)
	can <- if(isFALSE(rgl)) { " cannot " } else { " can " }
	msg <- paste0(
			"This mesh", can, "be converted to a 'rgl' mesh (see `?toRGL`).\n"
	)
	cat(msg)
	normals <- !is.null(x[["normals"]])
	has     <- if(normals) { " has " } else { " does not have " }
	msg     <- paste0("This mesh", has, "vertex normals.\n")
	cat(msg)
	invisible(NULL)
}

#' @title Make a 3D mesh
#' @description Make a 3D surface mesh from an input file,
#'   from an existing \strong{rgl} mesh object,
#'   or from given vertices and faces. The mesh is optionally cleaned:
#'   duplicated vertices or faces are merged, and isolated vertices are removed.
#'   The returned faces are coherently oriented, normals are computed if requested, and
#'   triangulation is performed if requested.
#'
#' @param x One of three options: 1) A numeric matrix with 3 columns providing the
#'   coordinates of the vertices of the mesh. 2) Either a list containing the components
#'   \code{vertices} and \code{faces}, or a \code{\link[rgl]{mesh3d}} object from
#'   package \strong{rgl}. 3) A filename to read a
#'   mesh file as in \code{\link[rgl]{readMeshFile}}.
#' @param faces When \code{x} is a numeric matrix with the vertices: Either an integer
#'   matrix (each row provides the vertex indices
#'   of the corresponding face) or a list of integer vectors, each one
#'   providing the vertex indices of the corresponding face.
#' @param triangulate Boolean. Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @param repairSoup Boolean. Whether to clean the mesh (merging duplicated
#'   vertices, duplicated faces, removing isolated vertices).
#' @param removeIntersections Boolean. Whether to attempt to remove self intersections.
#' @param removeMethod \code{character}. One of \code{"auto"} (for auto-refine)
#'   and \code{"auto_snap"} (auto-refine with iterative snap).
#' @param fillHoles Boolean. Whether to attempt to fill boundary holes.
#' @param fairHole Boolean. Use CGAL \code{triangulate_refine_and_fair_hole()}
#'   (\code{TRUE}) or \code{triangulate_and_refine_hole()} (\code{FALSE})?
#' @param maxNumHoles \code{integer}. Maximum number of holes to be filled. May be 0.
#' @param normals Boolean. Whether to compute vertex normals.
#' @param verbose Boolean. Whether to print out messages about mesh processing.
#'
#' @returns A list of class \code{CGALmesh} giving the vertices, the edges, the faces
#'   of the mesh, the exterior edges, the exterior vertices and optionally the
#'   vertex normals.
#'
#' @seealso See \code{\link[MeshUtils]{plotEdges}} for more details about the edges
#'   returned by this function.
#'   See \code{\link[MeshUtils]{makeMeshValid}} for a similar function that assumes
#'   that the input defines a valid mesh, and does not perform mesh repair to gain
#'   some speed.
#'   See \code{\link[MeshUtils]{toRGL}} for conversion to class
#'   \code{\link[rgl]{mesh3d}} from package \strong{rgl}.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' ## make mesh based on file
#' f_mesh    <- system.file("extdata", "dataHeart3.ply", package="MeshUtils")
#' mesh1     <- makeMesh(f_mesh, verbose=TRUE)
#' mesh1_rgl <- toRGL(mesh1)
#'
#' open3d(windowRect=c(50, 50, 562, 562))
#' wire3d(mesh1_rgl)
#'
#' ## make mesh from vertices and faces
#' ## tetrahedron with ill-oriented faces
#' vertices <- rbind(
#'   c(-1, -1, -1),
#'   c(1, 1, -1),
#'   c(1, -1, 1),
#'   c(-1, 1, 1))
#'
#' faces <- rbind(
#'   c(1, 2, 3),
#'   c(3, 4, 2),
#'   c(4, 2, 1),
#'   c(4, 3, 1))
#'
#' mesh2 <- makeMesh(vertices, faces, verbose=TRUE)
#' mesh2_rgl <- toRGL(mesh2)
#' open3d(windowRect=c(50, 50, 562, 562))
#' wire3d(mesh2_rgl)
#'
#' ## illustration of the `triangulate` option
#' ## the faces of the truncated icosahedron are hexagonal or pentagonal:
#' dataTruncIcosahedron[["faces"]]
#' # so we triangulate them:
#' mesh3     <- makeMesh(dataTruncIcosahedron, triangulate=TRUE)
#' mesh3_rgl <- toRGL(mesh3)
#' open3d(windowRect=c(50, 50, 562, 562))
#' wire3d(mesh3_rgl)
#'
#' @export
makeMesh <- function(x,
                     faces,
                     triangulate        =FALSE,
                     repairSoup         =TRUE,
                     removeIntersections=FALSE,
                     removeMethod       =c("auto", "auto_snap"),
                     fillHoles          =FALSE,
                     fairHole           =FALSE,
                     maxNumHoles        =5L,
                     normals            =FALSE,
                     verbose            =FALSE) {
    stopifnot(isBoolean(triangulate))
    stopifnot(isBoolean(repairSoup))
    stopifnot(isBoolean(removeIntersections))
    method_choices  <- c("auto", "auto_snap")
    removeMethod    <- match.arg(removeMethod, choices=method_choices)
    removeMethodInt <- match(removeMethod, method_choices)
    stopifnot(isBoolean(fillHoles))
    stopifnot(isBoolean(fairHole))
    stopifnot(isPositiveInteger(maxNumHoles))
    stopifnot(isBoolean(normals))
    stopifnot(isBoolean(verbose))
    mesh_cpp <- if(is.character(x)) { # filename
        stopifnot(length(x) == 1L, file.exists(x))
        makeMeshFF_cpp(x,
                       triangulate,
                       repairSoup,
                       removeIntersections,
                       removeMethodInt,
                       fillHoles,
                       fairHole,
                       maxNumHoles,
                       normals,
                       verbose)
    } else {
	      if(is.matrix(x)) {
		        vertices <- x
		        stopifnot(!missing(faces))
	      } else if(inherits(x, "mesh3d")) {
		    	  vft      <- getVFT(x, beforeCheck = TRUE)
		    	  mesh     <- vft[["rmesh"]]
		    	  vertices <- mesh[["vertices"]]
		    	  faces    <- mesh[["faces"]]
	      } else if(is.list(x)) {
		    	  stopifnot(hasName(x, "vertices"), hasName(x, "faces"))
		        vertices <- x[["vertices"]]
		        faces    <- x[["faces"]]
	      } else {
		        stop("`x` needs to be a matrix with vertex coords,\n or a `mesh3d` object, or a filename.")
	      }
		    ## ensure 0-based indexing, transposed vertices
		    mesh_r <- checkMesh(vertices, faces, aslist = TRUE)
		    makeMesh_cpp(mesh_r,
                     triangulate,
                     repairSoup,
                     removeIntersections,
                     removeMethodInt,
                     fillHoles,
                     fairHole,
                     maxNumHoles,
                     normals,
                     verbose)
	  }
		fromCPP(mesh_cpp)
}

#' @title Make a 3D mesh assuming valid input
#' @description Make a 3D surface mesh from an input file,
#'   from an existing \strong{rgl} mesh object, or
#'   from given vertices and faces - assuming that the input defines a valid mesh.
#'   Omitted validity checks save some processing time.
#'   The returned faces are coherently oriented (if possible),
#'   normals are computed if requested, triangulation is performed if requested
#'
#' @param x One of three options: 1) A numeric matrix with three columns providing the
#'   coordinates of the vertices of the mesh. 2) Either a list containing the components
#'   \code{vertices} and \code{faces}, or a \strong{rgl} mesh (i.e. a
#'   \code{\link[rgl]{mesh3d}} object). 3) A filename to read a
#'   mesh file as in \code{\link[rgl]{readMeshFile}}.
#' @param faces Either an integer matrix (each row provides the vertex indices
#'   of the corresponding face) or a list of integer vectors, each one
#'   providing the vertex indices of the corresponding face.
#' @param soup Boolean. Whether to assume a polygon soup. If \code{FALSE}, assume
#'   correctly oriented faces - which is a bit faster.
#' @param triangulate Boolean. Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @param normals Boolean. Whether to compute vertex normals.
#' @param verbose Boolean. Whether to print out messages about mesh processing.
#'
#' @returns A list of class \code{CGALmesh} giving the vertices, the edges, the faces
#'   of the mesh, the exterior edges, the exterior vertices and optionally the normals.
#'
#' @seealso See \code{\link[MeshUtils]{plotEdges}} for more details about the edges
#'   returned by this function.
#'   See \code{\link[MeshUtils]{makeMesh}} for a similar function that does not assume
#'   that the input defines a valid mesh, and performs mesh repair.
#'   See \code{\link[MeshUtils]{toRGL}} for conversion to class
#'   \code{\link[rgl]{mesh3d}} from package \strong{rgl}.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' ## mesh from file
#' f_mesh    <- system.file("extdata", "corner.off", package="MeshUtils")
#' mesh1     <- makeMeshValid(f_mesh, soup=TRUE, triangulate=TRUE)
#' mesh1_rgl <- toRGL(mesh1)
#'
#' open3d(windowRect=c(50, 50, 562, 562))
#' wire3d(mesh1_rgl)
#'
#' ## mesh from vertices, faces
#' head(dataPentaPrism[["vertices"]])
#' head(dataPentaPrism[["faces"]])
#' mesh2     <- makeMeshValid(dataPentaPrism, soup=TRUE, triangulate=TRUE)
#' mesh2_rgl <- toRGL(mesh2)
#' open3d(windowRect=c(50, 50, 562, 562))
#' wire3d(mesh2_rgl)
#'
#' @export
makeMeshValid <- function(x,
                          faces,
                          soup       =FALSE,
                          triangulate=FALSE,
                          normals    =FALSE,
                          verbose    =FALSE) {
    stopifnot(isBoolean(soup))
    stopifnot(isBoolean(triangulate))
    stopifnot(isBoolean(normals))
    stopifnot(isBoolean(verbose))
    mesh_cpp <- if(is.character(x)) { # filename
        stopifnot(length(x) == 1L, file.exists(x))
        makeMeshValidFF_cpp(x,
                            soup,
                            triangulate,
                            normals,
                            verbose)
    } else {
        if(is.matrix(x)) {
          vertices <- x
          stopifnot(!missing(faces))
        } else if(inherits(x, "mesh3d")) {
   	      vft      <- getVFT(x, beforeCheck = TRUE)
   	      mesh     <- vft[["rmesh"]]
   	      vertices <- mesh[["vertices"]]
   	      faces    <- mesh[["faces"]]
        } else if(is.list(x)) {
   	      stopifnot(hasName(x, "vertices"), hasName(x, "faces"))
          vertices <- x[["vertices"]]
          faces    <- x[["faces"]]
        } else {
          stop("`x` needs to be a matrix with vertex coords,\n or a `mesh3d` object, or a filename.")
        }
        ## ensure 0-based indexing, transposed vertices
        mesh_r <- checkMeshValid(vertices, faces, aslist = TRUE)
        makeMeshValid_cpp(mesh_r,
                          soup,
                          triangulate,
                          normals,
                          verbose)
    }
    fromCPP(mesh_cpp)
}

#' @title Conversion to 'rgl' mesh
#' @description Converts a \code{CGALmesh} object (the output of the \code{\link{makeMesh}}
#'   function) to a \code{\link[rgl]{mesh3d}} object from package \strong{rgl}.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#'   In order to be convertible to a \code{\link[rgl]{mesh3d}} object from package
#'   \strong{rgl}, its faces must have at most four sides.
#' @param ... Arguments passed to \code{\link[rgl]{mesh3d}}.
#'
#' @returns A \code{\link[rgl]{mesh3d}} object from package \strong{rgl}.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#' mesh     <- makeMesh(dataTruncIcosahedron, triangulate=TRUE)
#' mesh_rgl <- toRGL(mesh, segments=t(mesh[["edges"]]))
#' open3d(windowRect=c(50, 50, 562, 562), zoom=0.9)
#' wire3d(mesh_rgl, color="darkred")
#'
#' @export
#' @importFrom rgl mesh3d
toRGL <- function(x, ...) {
	if(!inherits(x, "CGALmesh")) {
		stop("The `x` argument must be of class 'CGALmesh'",
				 " (i.e., the output of the `makeMesh()` function).")
	}
	rgl <- attr(x, "toRGL")
	if(isFALSE(rgl)) {
		stop("Impossible to convert this mesh to a 'rgl' mesh ",
				 "(the faces must have at most four sides).")
	}
	if(rgl == 3L) {
		mesh3d(x        =x[["vertices"]],
				   normals  =x[["normals"]],
				   triangles=t(x[["faces"]]),
				   ...)
	} else if(rgl == 4L) {
		mesh3d(x      =x[["vertices"]],
				   normals=x[["normals"]],
				   quads  =t(x[["faces"]]),
				   ...)
	} else {
		faces <- split(x[["faces"]], lengths(x[["faces"]]))
		mesh3d(x        =x[["vertices"]],
				   normals  =x[["normals"]],
				   triangles=do.call(cbind, faces[["3"]]),
				   quads    =do.call(cbind, faces[["4"]]),
				   ...)
	}
}

#' @title Plot some edges
#' @description Plot the given edges with functions from package \strong{rgl}.
#'
#' @param vertices A matrix with 3 columns giving the coordinates of the vertices.
#' @param edges \code{integer} matrix with two columns giving the edges by pairs of
#'   vertex indices.
#' @param color \code{character}. Color for the edges.
#' @param lwd Positive number. Line width. Ignored if \code{edgesAsTubes=TRUE}.
#' @param edgesAsTubes Boolean. Whether to draw the edges as tubes.
#' @param tubesRadius Positive number. Radius of the tubes when \code{edgesAsTubes=TRUE}.
#' @param verticesAsSpheres Boolean. Whether to draw the vertices as spheres.
#' @param only \code{integer} vector made of the indices of the vertices you want
#'   to plot (as spheres). If missing, all vertices are plotted as spheres.
#' @param spheresRadius The radius of the spheres when
#'   \code{verticesAsSpheres=TRUE}.
#' @param spheresColor \code{character}. Color of the spheres when
#'   \code{verticesAsSpheres=TRUE}.
#'
#' @returns No value.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' # triangulate and plot the pentagrammic prism mesh
#' mesh     <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' mesh_rgl <- toRGL(mesh)
#' open3d(windowRect=c(50, 50, 562, 562), zoom=0.9)
#' shade3d(mesh_rgl, color="navy")
#' # plot the exterior edges only, given in `mesh[["exteriorEdges"]]`
#' plotEdges(mesh[["vertices"]],
#'           mesh[["exteriorEdges"]],
#'           color        ="gold",
#'           tubesRadius  =0.02,
#'           spheresRadius=0.02)
#'
#' @export
#' @importFrom rgl cylinder3d shade3d lines3d spheres3d
plotEdges <- function(
		vertices,
		edges,
		color = "black",
		lwd = 2,
		edgesAsTubes = TRUE,
		tubesRadius = 0.03,
		verticesAsSpheres = TRUE,
		only,
		spheresRadius = 0.05,
		spheresColor = color) {
  for(i in seq_len(nrow(edges))) {
		edge <- edges[i, ]
		if(edgesAsTubes) {
			tube <- cylinder3d(
					vertices[edge, , drop=FALSE], radius = tubesRadius, sides = 90)
			shade3d(tube, color = color)
		} else {
			lines3d(vertices[edge, ], color = color, lwd = lwd)
		}
	}
	if(verticesAsSpheres) {
		if(!missing(only)) {
			vertices <- vertices[only, , drop=FALSE]
		}
		spheres3d(vertices, radius=spheresRadius, color=spheresColor)
	}
	invisible(NULL)
}

#' @title Is mesh a valid mesh?
#' @description Is the given mesh a valid 3D surface mesh?
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns \code{TRUE} or \code{FALSE}.
#' @export
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(dataPentaPrism)
#' isValid(mesh)
#'
#' @export
isValid <- function(x) {
    if(!inherits(x, "CGALmesh")) {
        stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
    }
    meshCPP <- fromR(x)
    isValid_cpp(meshCPP)
}

#' @title Does mesh have garbage?
#' @description Does the given 3D surface mesh have garbage?
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns \code{TRUE} or \code{FALSE}.
#' @export
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(dataPentaPrism)
#' hasGarbage(mesh)
#'
#' @export
hasGarbage <- function(x) {
    if(!inherits(x, "CGALmesh")) {
        stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
    }
    meshCPP <- fromR(x)
    hasGarbage_cpp(meshCPP)
}

#' @title Is mesh a triangle mesh?
#' @description Is the given 3D surface mesh a triangle mesh?
#'
#' @param x A \code{list} with components \code{vertices} and \code{faces},
#'   e.g., a \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns \code{TRUE} or \code{FALSE}.
#' @export
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(dataPentaPrism, triangulate=FALSE)
#' isTriangle(mesh)
#'
#' @export
isTriangle <- function(x) {
  checkedMesh <- checkMesh(x[["vertices"]],
                           x[["faces"]],
                           aslist = TRUE)

  checkedMesh[["isTriangle"]]
}

#' @title Is mesh a quad mesh?
#' @description Is the given 3D surface mesh a quad mesh?
#'
#' @param x A \code{list} with components \code{vertices} and \code{faces},
#'   e.g., a \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns \code{TRUE} or \code{FALSE}.
#' @export
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' f_mesh <- system.file("extdata", "corner.off", package="MeshUtils")
#' mesh   <- makeMesh(f_mesh)
#' isQuad(mesh)
#'
#' @export
isQuad <- function(x) {
  checkedMesh <- checkMesh(x[["vertices"]],
                           x[["faces"]],
                           aslist = TRUE)

  checkedMesh[["isQuad"]]
}

#' @title Does mesh bound a volume?
#' @description Does given 3D surface mesh bound a volume?
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns \code{TRUE} or \code{FALSE}.
#' @export
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' doesBoundVolume(mesh)
#'
#' @export
doesBoundVolume <- function(x) {
    if(!inherits(x, "CGALmesh")) {
        stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
    }
    meshCPP <- fromR(x)
    doesBoundVolume_cpp(meshCPP)
}

#' @title Does mesh self intersect?
#' @description Does the given 3D surface mesh self intersect?
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns \code{TRUE} or \code{FALSE}.
#' @seealso See \code{\link[MeshUtils]{removeSelfIntersections}} for a function to
#'   remove self-intersections.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' doesSelfIntersect(mesh)
#'
#' @export
doesSelfIntersect <- function(x) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  meshCPP <- fromR(x)
  doesSelfIntersect_cpp(meshCPP)
}

#' @title Is mesh closed?
#' @description Is the given 3D surface mesh closed?
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns \code{TRUE} or \code{FALSE}.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#' @seealso See \code{\link[MeshUtils]{fillBoundaryHoles}} for a function to fill holes.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' isClosed(mesh)
#'
#' @export
isClosed <- function(x) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  meshCPP <- fromR(x)
  isClosed_cpp(meshCPP)
}

#' @title Orient mesh to bound a volume
#' @description Orient a given 3D surface mesh to bound a volume
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param normals Boolean. Whether to return vertex normals.
#' @returns A \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh   <- makeMesh(dataPentaPrism)
#' mesh_o <- orientToBoundVolume(mesh)
#' getVolume(mesh_o)
#'
#' @export
orientToBoundVolume <- function(x, normals = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(normals))
  meshCPP <- fromR(x)
  mesh    <- orientToBoundVolume_cpp(meshCPP, normals)
  fromCPP(mesh)
}

#' @title Get mesh area
#' @description Get the surface area of a given 3D surface mesh.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns \code{numeric}: The mesh area.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' getArea(mesh)
#'
#' @export
getArea <- function(x) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  meshCPP <- fromR(x)
  getArea_cpp(meshCPP)
}

#' @title Get axis-parallel bounding box
#' @description Get axis-parallel bounding box
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param triangulate Boolean. Whether to triangulate the faces of the bounding box.
#' @param normals Boolean. Whether to return vertex normals.
#'
#' @returns A \code{CGALmesh} object.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' mesh     <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' mesh_rgl <- toRGL(mesh)
#' bb       <- getBoundingBox(mesh)
#' bb_rgl   <- toRGL(bb)
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' wire3d(mesh_rgl)
#' wire3d(bb_rgl)
#'
#' @export
#' @importFrom rgl translate3d scale3d cube3d
getBoundingBox <- function(x, triangulate = FALSE, normals = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(triangulate))
  stopifnot(isBoolean(normals))
  meshCPP <- fromR(x)
  outL    <- boundingBox_cpp(meshCPP)
  lcorner <- outL[["lcorner"]]
  ucorner <- outL[["ucorner"]]
  center  <- (lcorner + ucorner) / 2
  ax <- ucorner[1L] - lcorner[1L]
  ay <- ucorner[2L] - lcorner[2L]
  az <- ucorner[3L] - lcorner[3L]
  m_rgl <- translate3d(scale3d(cube3d(), ax/2, ay/2, az/2),
                       center[1L], center[2L], center[3L])

  makeMesh(m_rgl, repairSoup=FALSE, triangulate=triangulate, normals=normals)
}

#' @title Get optimal bounding box
#' @description Get the optimal (oriented) bounding box of a given 3D surface mesh.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param triangulate Boolean. Whether to triangulate the faces of the bounding box.
#' @param normals Boolean. Whether to return vertex normals.
#' @returns A \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#' @details See \url{https://doc.cgal.org/latest/Optimal_bounding_box/} for details.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#' mesh     <- dataHeart1
#' mesh_rgl <- toRGL(mesh)
#' obb      <- getBoundingBoxOptimal(mesh)
#' obb_rgl  <- toRGL(obb[["mesh"]])
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' wire3d(mesh_rgl)
#' wire3d(obb_rgl)
#'
#' @export
getBoundingBoxOptimal <- function(x, triangulate = FALSE, normals = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(triangulate))
  stopifnot(isBoolean(normals))
  meshCPP <- fromR(x)
  outL    <- optimalBoundingBox_cpp(meshCPP, triangulate, normals)
  outL[["mesh"]] <- fromCPP(outL[["mesh"]])
  outL
}

#' @title Get mesh centroid
#' @description Get the centroid of a given 3D surface mesh.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}. The mesh must be triangle.
#' @returns \code{numeric} 3-vector with the coordinates of the mesh centroid.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' getCentroid(mesh)
#'
#' @export
getCentroid <- function(x) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  meshCPP <- fromR(x)
  getCentroid_cpp(meshCPP)
}

#' @title Get distance from points to a mesh
#' @description Get the Euclidean distance of points to a given 3D surface mesh.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#'   The mesh must be triangle or must be able to be made triangle.
#' @param points \code{numeric} matrix with 3 columns with one point per row.
#' @returns \code{numeric} vector: The distance of each point in \code{points}
#'   to the mesh \code{x}.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh   <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' points <- matrix(2*runif(3*4), ncol=3)
#' getDistance(mesh, points)
#'
#' @export
#' @importFrom stats na.omit
getDistance <- function(x, points) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  if(!is.matrix(points) || !is.numeric(points)) {
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  storage.mode(points) <- "double"
  n_pts <- nrow(points)
  is_na <- vapply(seq_len(n_pts), function(i) {
    anyNA(points[i, ]) }, logical(1))

  dst         <- rep(NA_real_, n_pts)
  pts_nona    <- na.omit(points)
  meshCPP     <- fromR(x)
  dst[!is_na] <- getDistance_cpp(meshCPP, t(pts_nona))
  dst
}

#' @title Get mesh volume
#' @description Get the volume of a given 3D surface mesh.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns \code{numeric}: The mesh volume - if mesh bounds a volume.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' getVolume(mesh)
#'
#' @export
getVolume <- function(x) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  meshCPP <- fromR(x)
  getVolume_cpp(meshCPP)
}

#' @title Add normals to a mesh
#' @description Add normal vectors to a given 3D surface mesh.
#'   Currently, only vertex normals are supported.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns A \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#' mesh        <- makeMesh(dataPentaPrism, triangulate=TRUE)
#' mesh_vn     <- assignNormals(mesh)
#' mesh_vn_rgl <- toRGL(mesh_vn)
#' open3d(windowRect=c(50, 50, 562, 562))
#' wire3d(mesh_vn_rgl)
#'
#' @export
assignNormals <- function(x) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  meshCPP <- fromR(x)
  mesh    <- addVNormals_cpp(meshCPP)
  fromCPP(mesh)
}

#' @title Triangulate mesh
#' @description Triangulate a given 3D surface mesh.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param normals Boolean. Whether to return vertex normals.
#' @returns A \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#' ## quad mesh
#' f_mesh   <- system.file("extdata", "corner.off", package="MeshUtils")
#' mesh     <- makeMesh(f_mesh)
#' mesh_rgl <- toRGL(mesh)
#' isTriangle(mesh)
#'
#' mesh_tri     <- triangulateMesh(mesh)
#' mesh_tri_rgl <- toRGL(mesh_tri)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' wire3d(mesh_rgl)
#' next3d()
#' wire3d(mesh_tri_rgl)
#'
#' @export
triangulateMesh <- function(x, normals = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(normals))
  meshCPP <- fromR(x)
  mesh    <- triangulateMesh_cpp(meshCPP, normals)
  fromCPP(mesh)
}
