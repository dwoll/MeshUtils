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
#' @param triangulate Boolean, whether to triangulate the faces. Ignored if faces
#'   are already triangle.
#' @param normals Boolean, whether to compute the normals.
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
#' mesh1_rgl_a <- tmesh3d(
#'   vertices   =t(vertices),
#'   indices    =t(faces),
#'   homogeneous=FALSE)
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
#' TruncatedIcosahedron[["faces"]]
#' # so we triangulate them:
#' mesh2 <- makeMesh(
#'   mesh       =TruncatedIcosahedron,
#'   triangulate=TRUE,
#'   normals    =FALSE)
#'
#' ## now we can plot the truncated icosahedron
#' mesh2_rgl <- toRGL(mesh2)
#' open3d(windowRect=c(50, 50, 562, 562), zoom=0.9)
#' shade3d(mesh2_rgl, color="orange")
#'
#' @export
makeMesh <- function(vertices,
                     faces,
                     mesh       =NULL,
                     triangulate=FALSE,
                     normals    =FALSE,
                     clean      =TRUE) {
	if(!is.null(mesh)) {
		if(inherits(mesh, "mesh3d")) {
			vft  <- getVFT(mesh, beforeCheck = TRUE)
			mesh <- vft[["rmesh"]]
		}
		vertices <- mesh[["vertices"]]
		faces    <- mesh[["faces"]]
	}
	checkedMesh      <- checkMesh(vertices, faces, aslist = TRUE)
	vertices         <- checkedMesh[["vertices"]]
	faces            <- checkedMesh[["faces"]]
	homogeneousFaces <- checkedMesh[["homogeneousFaces"]]
	isTriangle       <- checkedMesh[["isTriangle"]]
	if(triangulate && isTriangle) {
		message(
				"Ignored option `triangulate`, since the mesh is already triangulated.")
		triangulate <- FALSE
	}

	mesh_r   <- list("vertices"=vertices, "faces"=faces)
	mesh_cpp <- SurfEMesh_cpp(mesh_r, isTriangle, triangulate, clean, normals)
	fromCPP(mesh_cpp)
}

#' @title Conversion to 'rgl' mesh
#' @description Converts a \code{CGALmesh} object (the output of the \code{\link{makeMesh}}
#'   function) to a \strong{rgl} mesh object.
#'
#' @param mesh A \code{CGALmesh} object, i.e., a specific list as produced
#'   by the \code{\link{makeMesh}} function). In order to be
#'   convertible to a \strong{rgl} mesh, its faces must have at most four sides
#' @param ... Arguments passed to \code{\link[rgl]{mesh3d}}
#'
#' @returns A \strong{rgl} mesh object, i.e., a list of class \code{\link[rgl]{mesh3d}}.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' library(rgl)
#' mesh <- makeMesh(
#'   TruncatedIcosahedron[["vertices"]],
#'   TruncatedIcosahedron[["faces"]],
#'   triangulate=TRUE)
#'
#' mesh_rgl <- toRGL(mesh, segments=t(mesh[["edges"]]))
#' open3d(windowRect=c(50, 50, 562, 562), zoom=0.9)
#' shade3d(mesh_rgl, color="darkred")
#'
#' @export
#' @importFrom rgl mesh3d
toRGL <- function(mesh, ...) {
	if(!inherits(mesh, "CGALmesh")) {
		stop("The `mesh` argument must be of class 'CGALmesh'",
				 " (e.g. an output of the `Mesh` function).")
	}
	rgl <- attr(mesh, "toRGL")
	if(isFALSE(rgl)) {
		stop("Impossible to convert this mesh to a 'rgl' mesh ",
				 "(the faces must have at most four sides).")
	}
	if(rgl == 3L) {
		mesh3d(x        =mesh[["vertices"]],
				   normals  =mesh[["normals"]],
				   triangles=t(mesh[["faces"]]),
				   ...)
	} else if(rgl == 4L) {
		mesh3d(x      =mesh[["vertices"]],
				   normals=mesh[["normals"]],
				   quads  =t(mesh[["faces"]]),
				   ...)
	} else {
		faces <- split(mesh[["faces"]], lengths(mesh[["faces"]]))
		mesh3d(x        =mesh[["vertices"]],
				   normals  =mesh[["normals"]],
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
#' @param edgesAsTubes Boolean, whether to draw the edges as tubes.
#' @param tubesRadius The radius of the tubes when \code{edgesAsTubes=TRUE}.
#' @param verticesAsSpheres Boolean, whether to draw the vertices as spheres.
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
#' # we triangulate the truncated icosahedron mesh
#' mesh1 <- makeMesh(
#'   mesh       =TruncatedIcosahedron,
#'   triangulate=TRUE,
#'   normals    =FALSE)
#'
#' # now we can plot the truncated icosahedron
#' mesh1_rgl <- toRGL(mesh1)
#' open3d(windowRect=c(50, 50, 562, 562), zoom=0.9)
#' shade3d(mesh1_rgl, color="gold")
#'
#' # we triangulate the pentagrammic prism mesh
#' mesh2 <- makeMesh(
#'   mesh       =PentagrammicPrism,
#'   triangulate=TRUE,
#'    normals   = FALSE)
#'
#' # now we can plot the pentagrammic prism
#' mesh2_rgl <- toRGL(mesh2)
#' open3d(windowRect=c(50, 50, 562, 562), zoom=0.9)
#' shade3d(mesh2_rgl, color="navy")
#' # we plot the exterior edges only, given in `mesh2[["exteriorEdges"]]`
#' plotEdges(
#'   mesh2[["vertices"]],
#'   mesh2[["exteriorEdges"]],
#'   color        ="gold",
#'   tubesRadius  =0.02,
#'   spheresRadius=0.02)
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
		spheresColor = color
) {
  for(i in seq_len(nrow(edges))) {
		edge <- edges[i, ]
		if(edgesAsTubes) {
			tube <- cylinder3d(
					vertices[edge, ], radius = tubesRadius, sides = 90
			)
			shade3d(tube, color = color)
		}else{
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

#' @title Does mesh bound a volume?
#' @description Does mesh bound a volume?
#'
#' @param x A list with components \code{vertices} and \code{faces}, e.g., a \code{CGALmesh}
#'     object.
#' @returns TRUE or FALSE.
#' @export
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=PentagrammicPrism, triangulate=TRUE)
#' doesBoundVolume(mesh)
#'
#' @export
doesBoundVolume <- function(x) {
    cppMesh <- fromR(x)
    doesBoundVolume_cpp(cppMesh)
}

#' @title Does mesh self intersect?
#' @description Does mesh self intersect?
#'
#' @param x A list with components \code{vertices} and \code{faces}, e.g., a \code{CGALmesh}
#'     object.
#' @returns TRUE or FALSE.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=PentagrammicPrism, triangulate=TRUE)
#' doesSelfIntersect(mesh)
#'
#' @export
doesSelfIntersect <- function(x) {
  cppMesh <- fromR(x)
  doesSelfIntersect_cpp(cppMesh)
}

#' @title Is mesh closed?
#' @description Is mesh closed?
#'
#' @param x A list with components \code{vertices} and \code{faces}, e.g., a \code{CGALmesh}
#'     object.
#' @returns TRUE or FALSE.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=PentagrammicPrism, triangulate=TRUE)
#' isClosed(mesh)
#'
#' @export
isClosed <- function(x) {
  cppMesh <- fromR(x)
  isClosed_cpp(cppMesh)
}

#' @title Orient mesh to bound a volume
#' @description Orient mesh to bound a volume
#'
#' @param x A list with components \code{vertices} and \code{faces}, e.g., a \code{CGALmesh}
#'     object.
#' @returns \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh   <- makeMesh(mesh=PentagrammicPrism, triangulate=TRUE)
#' mesh_o <- orientToBoundVolume(mesh)
#' getVolume(mesh_o)
#'
#' @export
orientToBoundVolume <- function(x) {
  cppMesh <- fromR(x)
  mesh    <- orientToBoundVolume_cpp(cppMesh)
  fromCPP(mesh)
}

#' @title Remove self intersections
#' @description Remove self intersections
#'
#' @param x A list with components \code{vertices} and \code{faces}, e.g., a \code{CGALmesh}
#'     object.
#' @returns \code{CGALmesh} object.
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=PentagrammicPrism, triangulate=TRUE)
#' mesh_nsi <- removeSelfIntersections(mesh)
#' getVolume(mesh_nsi)
#'
#' @export
removeSelfIntersections <- function(x, method=c("auto", "auto_snap")) {
  method_choices <- c("auto", "auto_snap")
  method     <- match.arg(method, choices=method_choices)
  method_int <- match(method, method_choices)
  cppMesh    <- fromR(x)
  mesh       <- removeSelfIntersections_cpp(cppMesh, method_int)
  fromCPP(mesh)
}

#' @title Get mesh volume
#' @description Get the volume of a 3D mesh.
#'
#' @param x A list with components \code{vertices} and \code{faces}, e.g., a \code{CGALmesh}
#'     object.
#' @returns \code{numeric}: The mesh volume - if mesh bounds a volume.
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=PentagrammicPrism, triangulate=TRUE)
#' getVolume(mesh)
#'
#' @export
getVolume <- function(x) {
  cppMesh <- fromR(x)
  getVolume_cpp(cppMesh)
}

#' @title Get mesh centroid
#' @description Get mesh centroid
#'
#' @param x A list with components \code{vertices} and \code{faces}, e.g., a \code{CGALmesh}
#'     object.
#' @returns numeric 3-vector with mesh centroid.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh <- makeMesh(mesh=PentagrammicPrism, triangulate=TRUE)
#' getCentroid(mesh)
#'
#' @export
getCentroid <- function(x) {
  cppMesh <- fromR(x)
  getCentroid_cpp(cppMesh)
}

#' @title Get optimal bounding box
#' @description Get oriented bounding box
#'
#' @param x A list with components \code{vertices} and \code{faces}, e.g., a \code{CGALmesh}
#'     object.
#' @returns A \code{CGALmesh} object.
#'
#' @examples
#' library(MeshUtils)
#' mesh     <- makeMesh(mesh=PentagrammicPrism, triangulate=TRUE)
#' mesh_rgl <- toRGL(mesh)
#' obb      <- getOptimalBoundingBox(mesh)
#' obb_rgl  <- toRGL(obb[["mesh"]])
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' wire3d(mesh_rgl)
#' wire3d(obb_rgl)
#'
#' @export
getOptimalBoundingBox <- function(x) {
  cppMesh <- fromR(x)
  outL    <- optimalBoundingBox_cpp(cppMesh)
  outL[["mesh"]] <- fromCPP(outL[["mesh"]])
  outL
}

#' @title Get axis-parallel bounding box
#' @description Get axis-parallel bounding box
#'
#' @param x A list with components \code{vertices} and \code{faces}, e.g., a \code{CGALmesh}
#'     object.
#' @param out Character to indicate output mesh format.
#'
#' @returns A \code{CGALmesh} object or a \code{\link[rgl]{mesh3d}} object from package \strong{rgl}.
#'
#' @author Originally developed by Stephane Laurent, adapted by Daniel Wollschlaeger.
#'
#' @examples
#' library(MeshUtils)
#' mesh     <- makeMesh(mesh=PentagrammicPrism, triangulate=TRUE)
#' mesh_rgl <- toRGL(mesh)
#' bb_rgl   <- getBoundingBox(mesh, out="rgl")
#' open3d(windowRect=50 + c(0, 0, 800, 400))
#' wire3d(mesh_rgl)
#' wire3d(bb_rgl)
#'
#' @export
#' @importFrom rgl translate3d scale3d cube3d
getBoundingBox <- function(x, out=c("CGALmesh", "rgl")) {
  out     <- match.arg(out)
  cppMesh <- fromR(x)
  outL    <- boundingBox_cpp(cppMesh)
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
    makeMesh(mesh=m_rgl)
  }
}
