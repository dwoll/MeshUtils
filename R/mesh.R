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
	nf  <- if(is.list(x[["faces"]])) length(x[["faces"]]) else nrow(x[["faces"]])
	msg <- sprintf("CGALmesh with %d vertices and %d faces.\n", nv, nf)
	cat(msg)
	elr <- formatC(range(x[["edgesDF"]][["length"]]))
	msg <- sprintf("The edge lengths vary from %s to %s.\n", elr[1L], elr[2L])
	cat(msg)
	is  <- if(rgl == 3L) " is " else " is not "
	msg <- paste0("This mesh", is, "triangle.\n")
	cat(msg)
	can <- if(isFALSE(rgl)) " cannot " else " can "
	msg <- paste0(
			"This mesh", can, "be converted to a 'rgl' mesh (see `?toRGL`).\n"
	)
	cat(msg)
	normals <- !is.null(x[["normals"]])
	has     <- if(normals) " has " else " does not have "
	msg     <- paste0("This mesh", has, "vertex normals.\n")
	cat(msg)
	invisible(NULL)
}

#' @title Make a 3D mesh
#' @description Make a 3D mesh from given vertices and faces; the returned
#'   faces are coherently oriented, normals are computed if desired, and
#'   triangulation is performed if desired. The mesh is also cleaned:
#'   duplicated vertices or faces are merged, and isolated vertices are removed.
#'
#' @param vertices A numeric matrix with three columns.
#' @param faces Either an integer matrix (each row provides the vertex indices
#'   of the corresponding face) or a list of integer vectors, each one
#'   providing the vertex indices of the corresponding face.
#' @param mesh If not \code{NULL}, this argument takes precedence over \code{vertices}
#'   and \code{faces}, and must be either a list containing the components \code{vertices}
#'   and \code{faces} (objects as described above), otherwise a \strong{rgl} mesh
#'   (i.e. a \code{\link[rgl]{mesh3d}} object).
#' @param triangulate Boolean: Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @param repairSoup Boolean: Whether to do some mesh cleaning.
#' @param removeIntersections Boolean: Whether to attempt to remove self intersections.
#' @param method One of \code{"auto"} (for auto-refine) and \code{"auto_snap"} (auto-refine with iterative snap).
#' @param fillHoles Boolean: Whether to attempt to fill boundary holes.
#' @param fairHole Boolean: Use CGAL `triangulate_refine_and_fair_hole()` (\code{TRUE})
#'     or `triangulate_and_refine_hole()` (\code{FALSE})?
#' @param maxNumHoles \code{integer}: Maximum number of holes to be filled. May be 0.
#' @param normals Boolean: Whether to compute the normals.
#'
#' @returns A list of class \code{CGALmesh} giving the vertices, the edges, the faces
#'   of the mesh, the exterior edges, the exterior vertices and optionally the normals.
#'
#' @seealso See \code{\link[MeshUtils]{plotEdges}} for more details about the edges
#'   returned by this function. See \code{\link[MeshUtils]{toRGL}} for conversion to class
#'   \code{\link[rgl]{mesh3d}} from package \strong{rgl}.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' ## a tetrahedron with ill-oriented faces
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
#' ## plot the tetrahedron, hiding the back of the faces
#' ## then some faces do not appear, as their orientation is not correct
#' mesh1_rgl_a <- tmesh3d(vertices=t(vertices),
#'                        indices =t(faces))
#'
#' open3d(windowRect=c(50, 50, 562, 562))
#' shade3d(mesh1_rgl_a, color="green", back = "cull")
#'
#' ## now run the `makeMesh` function
#' mesh1 <- makeMesh(vertices, faces, normals=FALSE)
#' ## plot the tetrahedron, hiding the back of the faces
#' ## then all faces appear now
#' mesh1_rgl_b <- toRGL(mesh1)
#' open3d(windowRect=c(50, 50, 562, 562))
#' shade3d(mesh1_rgl_b, color="blue", back="cull")
#'
#' ## illustration of the `triangulate` option
#' ## the faces of the truncated icosahedron are hexagonal or pentagonal:
#' dataTruncIcosahedron[["faces"]]
#' # so we triangulate them:
#' mesh2     <- makeMesh(mesh=dataTruncIcosahedron, triangulate=TRUE)
#' mesh2_rgl <- toRGL(mesh2)
#' open3d(windowRect=c(50, 50, 562, 562), zoom=0.9)
#' shade3d(mesh2_rgl, color="orange")
#'
#' @export
makeMesh <- function(vertices,
                     faces,
                     mesh               =NULL,
                     triangulate        =FALSE,
                     repairSoup         =TRUE,
                     removeIntersections=FALSE,
                     removeMethod       =c("auto", "auto_snap"),
                     fillHoles          =FALSE,
                     fairHole           =FALSE,
                     maxNumHoles        =0L,
                     normals            =FALSE) {
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
	if(!is.null(mesh)) {
		if(inherits(mesh, "mesh3d")) {
			vft  <- getVFT(mesh, beforeCheck = TRUE)
			mesh <- vft[["rmesh"]]
		}
		vertices <- mesh[["vertices"]]
		faces    <- mesh[["faces"]]
	}
	checkedMesh <- checkMesh(vertices, faces, aslist = TRUE)
	vertices    <- checkedMesh[["vertices"]]
	faces       <- checkedMesh[["faces"]]
	isTriangle  <- checkedMesh[["isTriangle"]]
	if(triangulate && isTriangle) {
		message("Ignored option `triangulate` as mesh is already triangle")
		triangulate <- FALSE
	}

	mesh_r   <- list("vertices"=vertices, "faces"=faces)
	mesh_cpp <- makeMesh_cpp(mesh_r,
                           triangulate,
	                         repairSoup,
	                         removeIntersections,
													 removeMethodInt,
	                         fillHoles,
	                         fairHole,
													 maxNumHoles,
											     normals)
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
#' @returns A \strong{rgl} mesh object, i.e., a list of class \code{\link[rgl]{mesh3d}}.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#' mesh     <- makeMesh(mesh=dataTruncIcosahedron, triangulate=TRUE)
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
#' @description Plot the given edges with \strong{rgl}.
#'
#' @param vertices A three-columns matrix giving the coordinates of the vertices.
#' @param edges A two-columns integer matrix giving the edges by pairs of
#'   vertex indices.
#' @param color A color for the edges.
#' @param lwd Line width, a positive number, ignored if \code{edgesAsTubes=TRUE}.
#' @param edgesAsTubes Boolean: Whether to draw the edges as tubes.
#' @param tubesRadius The radius of the tubes when \code{edgesAsTubes=TRUE}.
#' @param verticesAsSpheres Boolean: Whether to draw the vertices as spheres.
#' @param only Integer vector made of the indices of the vertices you want
#'   to plot (as spheres), or \code{NULL} to plot all vertices.
#' @param spheresRadius The radius of the spheres when
#'   \code{verticesAsSpheres=TRUE}.
#' @param spheresColor The color of the spheres when
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
#' # triangulate the pentagrammic prism mesh
#' mesh <- makeMesh(
#'   mesh       =dataPentaPrism,
#'   triangulate=TRUE,
#'   normals    =FALSE)
#'
#' # plot the pentagrammic prism
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
		only = NULL,
		spheresRadius = 0.05,
		spheresColor = color) {
  for(i in seq_len(nrow(edges))) {
		edge <- edges[i, ]
		if(edgesAsTubes) {
			tube <- cylinder3d(
					vertices[edge, ], radius = tubesRadius, sides = 90)
			shade3d(tube, color = color)
		} else {
			lines3d(vertices[edge, ], color = color, lwd = lwd)
		}
	}
	if(verticesAsSpheres) {
		if(!is.null(only)) {
			vertices <- vertices[only, , drop=FALSE]
		}
		spheres3d(vertices, radius=spheresRadius, color=spheresColor)
	}
	invisible(NULL)
}

#' @title Is mesh a valid mesh?
#' @description Is the given mesh a valid surface mesh?
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns TRUE or FALSE.
#' @export
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=dataPentaPrism)
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
#' @description Does the given surface mesh have garbage?
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns TRUE or FALSE.
#' @export
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=dataPentaPrism)
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
#' @description Is the given surface mesh a triangle mesh?
#'
#' @param x A \code{list} with components \code{vertices} and \code{faces},
#'   e.g., a \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns TRUE or FALSE.
#' @export
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=dataPentaPrism, triangulate=FALSE)
#' isTriangle(mesh)
#'
#' @export
isTriangle <- function(x) {
  checkedMesh <- checkMesh(x[["vertices"]],
                           x[["faces"]],
                           aslist = TRUE)

  checkedMesh[["isTriangle"]]
}

#' @title Does mesh bound a volume?
#' @description Does mesh bound a volume?
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param triangulate Boolean: Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @returns TRUE or FALSE.
#' @export
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' doesBoundVolume(mesh)
#'
#' @export
doesBoundVolume <- function(x, triangulate = FALSE) {
    if(!inherits(x, "CGALmesh")) {
        stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
    }
    stopifnot(isBoolean(triangulate))
    meshCPP <- fromR(x)
    doesBoundVolume_cpp(meshCPP, triangulate)
}

#' @title Does mesh self intersect?
#' @description Does mesh self intersect?
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param triangulate Boolean: Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @returns TRUE or FALSE.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' doesSelfIntersect(mesh)
#'
#' @export
doesSelfIntersect <- function(x, triangulate = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(triangulate))
  meshCPP <- fromR(x)
  doesSelfIntersect_cpp(meshCPP, triangulate)
}

#' @title Is mesh closed?
#' @description Is mesh closed?
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns TRUE or FALSE.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
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
#' @description Orient mesh to bound a volume
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param triangulate Boolean: Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @returns \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh   <- makeMesh(mesh=dataPentaPrism, triangulate=FALSE)
#' mesh_o <- orientToBoundVolume(mesh, triangulate=TRUE)
#' getVolume(mesh_o)
#'
#' @export
orientToBoundVolume <- function(x, triangulate = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(triangulate))
  meshCPP <- fromR(x)
  mesh    <- orientToBoundVolume_cpp(meshCPP, triangulate)
  fromCPP(mesh)
}

#' @title Remove self intersections
#' @description Remove self intersections.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param triangulate Boolean: Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @param method One of \code{"auto"} (for auto-refine) and \code{"auto_snap"} (auto-refine with iterative snap).
#' @returns \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#' @details See \url{https://www.cgal.org/2025/06/13/autorefine-and-snap/}.
#'
#' @examples
#' library(MeshUtils)
#' mesh     <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' mesh_nsi <- removeSelfIntersections(mesh)
#' getVolume(mesh_nsi)
#'
#' @export
removeSelfIntersections <- function(x, triangulate = FALSE, method=c("auto", "auto_snap")) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(triangulate))
  method_choices <- c("auto", "auto_snap")
  method    <- match.arg(method, choices=method_choices)
  methodInt <- match(method, method_choices)
  meshCPP   <- fromR(x)
  mesh      <- removeSelfIntersections_cpp(meshCPP, triangulate, methodInt)
  fromCPP(mesh)
}

#' @title Fill boundary holes
#' @description Fill boundary holes.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param fairHole Boolean: Use CGAL `triangulate_refine_and_fair_hole()` (\code{TRUE})
#'     or `triangulate_and_refine_hole()` (\code{FALSE})?
#' @param maxNumHoles \code{integer}: Maximum number of holes to be filled.
#' @returns \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#' @details See \url{https://www.cgal.org/2025/06/13/autorefine-and-snap/}.
#'     If faces are not triangle, the mesh is triangulated.
#'
#' @examples
#' library(MeshUtils)
#' mesh      <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' mesh_fill <- fillBoundaryHoles(mesh)
#'
#' @export
fillBoundaryHoles <- function(x, fairHole = TRUE, maxNumHoles=10L) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(fairHole))
  stopifnot(isStrictPositiveInteger(maxNumHoles))
  meshCPP <- fromR(x)
  mesh    <- fillBoundaryHoles_cpp(meshCPP, fairHole, maxNumHoles)
  fromCPP(mesh)
}

#' @title Get mesh area
#' @description Get the surface area of a 3D mesh.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param triangulate Boolean: Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @returns \code{numeric}: The mesh area.
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' getArea(mesh)
#'
#' @export
getArea <- function(x, triangulate = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(triangulate))
  meshCPP <- fromR(x)
  getArea_cpp(meshCPP, triangulate)
}

#' @title Get mesh volume
#' @description Get the volume of a 3D mesh.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param triangulate Boolean: Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @returns \code{numeric}: The mesh volume - if mesh bounds a volume.
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' getVolume(mesh)
#'
#' @export
getVolume <- function(x, triangulate = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(triangulate))
  meshCPP <- fromR(x)
  getVolume_cpp(meshCPP, triangulate)
}

#' @title Get mesh centroid
#' @description Get mesh centroid.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}. The mesh must be triangle.
#' @param triangulate Boolean: Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @returns \code{numeric} 3-vector with the cartesian coordinates of the mesh centroid.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' getCentroid(mesh)
#'
#' @export
getCentroid <- function(x, triangulate = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isBoolean(triangulate))
  meshCPP <- fromR(x)
  getCentroid_cpp(meshCPP, triangulate)
}

#' @title Get optimal bounding box
#' @description Get oriented bounding box
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @returns A \code{CGALmesh} object.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#' mesh     <- dataHeart1
#' mesh_rgl <- toRGL(mesh)
#' obb      <- getOptimalBoundingBox(mesh)
#' obb_rgl  <- toRGL(obb[["mesh"]])
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' wire3d(mesh_rgl)
#' wire3d(obb_rgl)
#'
#' @export
getOptimalBoundingBox <- function(x) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  meshCPP <- fromR(x)
  outL    <- optimalBoundingBox_cpp(meshCPP)
  outL[["mesh"]] <- fromCPP(outL[["mesh"]])
  outL
}

#' @title Get axis-parallel bounding box
#' @description Get axis-parallel bounding box
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param out Character to indicate output mesh format.
#'
#' @returns A \code{CGALmesh} object or a \code{\link[rgl]{mesh3d}} object from package \strong{rgl}.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' mesh     <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' mesh_rgl <- toRGL(mesh)
#' bb_rgl   <- getBoundingBox(mesh, out="rgl")
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' wire3d(mesh_rgl)
#' wire3d(bb_rgl)
#'
#' @export
#' @importFrom rgl translate3d scale3d cube3d
getBoundingBox <- function(x, out=c("CGALmesh", "rgl")) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  out     <- match.arg(out)
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

  if(out == "rgl") {
    m_rgl
  } else {
    makeMesh(mesh=m_rgl, repairSoup=FALSE)
  }
}

#' @title Get distance from points to a mesh
#' @description Get the Euclidean distance of points
#' to a 3D mesh.
#'
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}. The mesh must be triangle.
#' @param points \code{numeric} matrix with 3 columns with one point per row.
#' @param triangulate Boolean: Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @returns \code{numeric} vector: The distance of each point in \code{points}
#'     to the mesh \code{x}.
#' @examples
#' library(MeshUtils)
#' mesh   <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' points <- matrix(2*runif(3*4), ncol=3)
#' getDistance(mesh, points)
#'
#' @export
#' @importFrom stats na.omit
getDistance <- function(x, points, triangulate = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  if(!is.matrix(points) || !is.numeric(points)) {
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  stopifnot(isBoolean(triangulate))
  if((ncol(points) != 3L) && !triangulate) {
    stop("The `points` matrix must have three columns. Need 'triangulate=TRUE'", call. = TRUE)
  }
  storage.mode(points) <- "double"
  n_pts <- nrow(points)
  is_na <- vapply(seq_len(n_pts), function(i) {
    anyNA(points[i, ]) }, logical(1))

  dst         <- rep(NA_real_, n_pts)
  pts_nona    <- na.omit(points)
  meshCPP     <- fromR(x)
  dst[!is_na] <- getDistance_cpp(meshCPP, t(pts_nona), triangulate)
  dst
}

#' @title Hausdorff distance between two meshes
#' @description Hausdorff distance between two meshes. Either
#'   approximate distance, or distance estimate with a given error bound.
#' @param mesh1 A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param mesh2 A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param triangulate1 Boolean: Whether to triangulate the faces of \code{mesh1}.
#'   Ignored if faces are already triangle.
#' @param triangulate2 Boolean: Whether to triangulate the faces of \code{mesh2}.
#'   Ignored if faces are already triangle.
#' @param symmetric Boolean: Whether to consider the symmetric Hausdorff
#'   distance.
#' @param errorBound A positive number, a bound on the error of the
#'   estimate. If missing, the approximate distance is returned.
#' @returns A number. For the apprixmate distance, the algorithm uses
#'   some simulations and thus the result can vary.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#' @examples
#' library(MeshUtils)
#' ## approximate symmetric Hausdorff distance
#' getHausdorffDistance(dataHeart1, dataHeart2, symmetric=TRUE)
#'
#' ## estimate with error bound
#' getHausdorffDistance(dataHeart1, dataHeart2, symmetric=TRUE,
#'                      errorBound=0.001)
#'
#' @export
getHausdorffDistance <- function(mesh1, mesh2, symmetric = TRUE, errorBound,
	triangulate1 = FALSE, triangulate2 = FALSE) {
  stopifnot(inherits(mesh1, "CGALmesh"))
  stopifnot(inherits(mesh2, "CGALmesh"))
  stopifnot(isBoolean(symmetric))
  stopifnot(isBoolean(triangulate1))
  stopifnot(isBoolean(triangulate1))
  meshCPP1 <- fromR(mesh1)
  meshCPP2 <- fromR(mesh2)
  if(!missing(errorBound)) {
    stopifnot(isPositiveNumber(errorBound))
    getHausdorffEst_cpp(meshCPP1, meshCPP2, symmetric, errorBound,
                        triangulate1, triangulate2)
  } else {
    getHausdorffApprox_cpp(meshCPP1, meshCPP2, symmetric,
    	                     triangulate1, triangulate2)
  }
}

#' @title Isotropic remeshing
#' @description Isotropic remeshing.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#' @param method \code{character}: Either \code{"uniform"} for uniform sizing field
#'    or \code{"adaptive"} for adapative sizing field.
#' @param targetEdgeLen Positive number for \code{method="uniform"}:
#'   The target edge length of the remeshed mesh.
#' @param tol Positive number for \code{method="adaptive"}:
#'   Error tolerance. See details
#' @param edgeMin Positive number for \code{method="adaptive"}:
#'   Minimum edge length. See details.
#' @param edgeMax Positive number for \code{method="adaptive"}:
#'   Maximum edge length. See details.
#' @param nIter Positive \code{integer}: Number of iterations.
#' @param nRelaxSteps Positive \code{integer}: Number of relaxation steps.
#' @param out \code{haracter: Indicate output mesh format.
#' @return A \code{CGALmesh} object or a \code{\link[rgl]{mesh3d}} object from package \strong{rgl},
#'   depending on option \code{out}.
#' @details See \url{https://doc.cgal.org/latest/PMP_Remeshing/classCGAL_1_1Polygon__mesh__processing_1_1Adaptive__sizing__field.html},
#'   \url{https://doc.cgal.org/latest/PMP_Remeshing/group__PMP__local__remeshing__grp.html#ga412f696ec3009074bf957f1bba638248}.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' mesh         <- makeMesh(mesh=dataTruncIcosahedron, triangulate=TRUE)
#' mesh_rgl     <- toRGL(mesh)
#' mesh_rem_rgl <- remeshIsotropic(mesh, targetEdgeLen=1, out="rgl")
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' wire3d(mesh_rgl)
#' next3d()
#' wire3d(mesh_rem_rgl)
#'
#' @export
remeshIsotropic <- function(
        x,
        method = c("uniform", "adaptive"),
        tol = 0.001,
        edgeMin = 0.001,
        edgeMax,
        targetEdgeLen,
        nIter = 1,
        nRelaxSteps = 1,
        out = c("CGALmesh", "rgl")) {
    method <- match.arg(method)
    if(!inherits(x, "CGALmesh")) {
        stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
    }
    out <- match.arg(out)
    stopifnot(isStrictPositiveInteger(nIter))
    stopifnot(isStrictPositiveInteger(nRelaxSteps))
    meshCPP <- fromR(x)
    if(method == "uniform") {
      stopifnot(isPositiveNumber(targetEdgeLen))
      meshRe <- remeshIsotropicUniform_cpp(
          meshCPP,
          targetEdgeLen,
          as.integer(nIter),
          as.integer(nRelaxSteps))
    } else if(method == "adaptive") {
      stopifnot(isPositiveNumber(tol))
      stopifnot(isPositiveNumber(edgeMin))
      stopifnot(isPositiveNumber(edgeMax))

      meshRe <- remeshIsotropicAdapt_cpp(
          meshCPP,
          tol,
          edgeMin,
          edgeMax,
          as.integer(nIter),
          as.integer(nRelaxSteps))
    }

    meshReWN <- vcgUpdateNormals(
        tmesh3d(meshRe[["vertices"]],
                meshRe[["faces"]],
                normals=NULL))
    meshOut <- if(out == "rgl") {
        meshReWN
    } else {
        makeMesh(mesh=meshReWN)
    }
    meshOut
}

#' @title Catmull-Clark subdivision and deformation
#' @description Performs the Catmull-Clark subdivision and deformation.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#'   The mesh must be triangle or be able to be made triangle.
#' @param nIter Positive \code{integer}: Number of iterations.
#' @param triangulate Boolean: Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @returns A \code{CGALmesh} object.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' mesh        <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' mesh_rgl    <- toRGL(mesh)
#' mesh_sd     <- subdivideCatmullClark(mesh, nIter=2)
#' mesh_sd_rgl <- toRGL(mesh_sd)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_rgl)
#' next3d()
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_sd_rgl)
#'
#' @export
subdivideCatmullClark <- function(x, nIter = 1, triangulate = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			     " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isBoolean(triangulate))
  meshCPP <- fromR(x)
  meshOut <- subdivideCatmullClark_cpp(meshCPP, as.integer(nIter), triangulate)
  fromCPP(meshOut)
}

#' @title Doo-Sabin subdivision and deformation
#' @description Performs the Doo-Sabin subdivision and deformation.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#'   The mesh must be triangle or be able to be made triangle.
#' @param nIter Positive \code{integer}: Number of iterations.
#' @param triangulate Boolean: Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @returns A \code{CGALmesh} object.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' mesh        <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' mesh_rgl    <- toRGL(mesh)
#' mesh_ds     <- subdivideDooSabin(mesh, nIter=2, triangulate=TRUE)
#' mesh_ds_rgl <- toRGL(mesh_ds)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_rgl)
#' next3d()
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_ds_rgl)
#'
#' @export
subdivideDooSabin <- function(x, nIter = 1, triangulate = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			     " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isBoolean(triangulate))
  meshCPP <- fromR(x)
  meshOut <- subdivideDooSabin_cpp(meshCPP, as.integer(nIter), triangulate)
  fromCPP(meshOut)
}

#' @title Sqrt3 subdivision and deformation
#' @description Performs the 'Sqrt3' subdivision and deformation.
#' @param x A \code{CGALmesh} object, i.e., the output of \code{\link[MeshUtils]{makeMesh}}.
#'   The mesh must be triangle or be able to be made triangle.
#' @param nIter Positive \code{integer}: Number of iterations.
#' @param triangulate Boolean: Whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @returns A \code{CGALmesh} object.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#'
#' mesh        <- makeMesh(mesh=dataPentaPrism, triangulate=TRUE)
#' mesh_rgl    <- toRGL(mesh)
#' mesh_s3     <- subdivideSqrt3(mesh, nIter=2)
#' mesh_s3_rgl <- toRGL(mesh_s3)
#'
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_rgl)
#' next3d()
#' view3d(0, 0, zoom=0.9)
#' wire3d(mesh_s3_rgl)
#'
#' @export
subdivideSqrt3 <- function(x, nIter = 1, triangulate = FALSE) {
  if(!inherits(x, "CGALmesh")) {
      stop("The `x` argument must be of class 'CGALmesh'",
			       " (i.e., the output of the `makeMesh()` function).")
  }
  stopifnot(isStrictPositiveInteger(nIter))
  stopifnot(isBoolean(triangulate))
  meshCPP <- fromR(x)
  meshOut <- subdivideSqrt3_cpp(meshCPP, as.integer(nIter), triangulate)
  fromCPP(meshOut)
}
